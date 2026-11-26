subroutine mprove(h,ch,oldrhs,dtime) 
!     --------------                                                    
  USE header, only : NI,NJ,rc_kind
      integer iter,i,j,l,im1,ip1,ci,cj,nci,ib,ie,is 
      REAL(kind=rc_kind) :: ch(9,NI,NJ),rhs(NI,NJ),oldrhs(NI,NJ),            &
     &     he(0:NI+1,0:NJ+1),h(0:NI+1,0:NJ+1),res,                      &
     &     dtime,tol,maxres,hstar,rlx                     
             
                                                                        
      tol= 1.d-10 
      do j=0,NJ+1 
         do i=0,NI+1 
            he(i,j)= 0.d0 
         end do 
      end do 
                                                                        
!     Form the new RHS : A*(h +dh) - b                                  
      do j=1,NJ 
         i=1 
         rhs(i,j) =(ch(1,i,j)*h(i,j)+ch(2,i,j)*h(i+1,j)+ch(3,i,j)*h(NI,j)+ch(4,i,j)*h(i,j+1)+ch(5,i,j)*h(i,j-1) &
                &  +ch(6,i,j)*h(i+1,j+1)+ch(7,i,j)*h(NI,j+1)+ch(8,i,j)*h(i+1,j-1)+ch(9,i,j)*h(NI,j-1))-oldrhs(i,j)
 
         do i=2,NI-1 
            rhs(i,j) =                                                  &
     &           (ch(1,i,j)*h(i,j)                                      &
     &           +ch(2,i,j)*h(i+1,j)                                    &
     &           +ch(3,i,j)*h(i-1,j)                                    &
     &           +ch(4,i,j)*h(i,j+1)                                    &
     &           +ch(5,i,j)*h(i,j-1)                                    &
     &           +ch(6,i,j)*h(i+1,j+1)                                  &
     &           +ch(7,i,j)*h(i-1,j+1)                                  &
     &           +ch(8,i,j)*h(i+1,j-1)                                  &
     &           +ch(9,i,j)*h(i-1,j-1))                                 &
     &           -oldrhs(i,j)                                           
         end do 
         i=NI 
         rhs(i,j) =                                                     &
     &        (ch(1,i,j)*h(i,j)                                         &
     &        +ch(2,i,j)*h(1,j)                                         &
     &        +ch(3,i,j)*h(i-1,j)                                       &
     &        +ch(4,i,j)*h(i,j+1)                                       &
     &        +ch(5,i,j)*h(i,j-1)                                       &
     &        +ch(6,i,j)*h(1,j+1)                                       &
     &        +ch(7,i,j)*h(i-1,j+1)                                     &
     &        +ch(8,i,j)*h(1,j-1)                                       &
     &        +ch(9,i,j)*h(i-1,j-1))                                    &
     &        -oldrhs(i,j)                                              
      end do 

                                                                        
!     Now solve for the error : he                                      
!     -------------------------                                         
                                                                        

  ! colours in i for the multi-colour SOR sweep (x periodic: odd NI needs a third colour)
  nci= 2
  if (mod(NI,2).ne.0) nci= 3
  if (NI.eq.1) nci= 1

 rlx= 1.72 
 do iter=1,1000 

   ! Multi-colour SOR sweep: the 9-point stencil couples (i,j) to (i+-1,j+-1), so the points are
   ! split into colours by the parity of i and j (a third i-colour for i=NI if NI is odd, as x is
   ! periodic). Points of one colour are independent and are updated in parallel; the result does
   ! not depend on the number of threads (the earlier lexicographic sweep was a data race with OpenMP).
   do cj=1,2
   do ci=1,nci
     if (ci.le.2) then
       ib= ci; ie= NI; is= 2
       if (nci.eq.3) ie= NI-1
     else
       ib= NI; ie= NI; is= 1
     endif
!$OMP PARALLEL DO PRIVATE(j,i,im1,ip1,hstar) SHARED(ch,he,rhs,rlx,ib,ie,is,cj)
     do j=cj,NJ,2
       do i=ib,ie,is
         im1= i-1
         ip1= i+1
         if (i.eq.1) im1= NI
         if (i.eq.NI) ip1= 1
         hstar= ( ch(2,i,j)*he(ip1,j)     &
             &   +ch(3,i,j)*he(im1,j)     &
             &   +ch(4,i,j)*he(i,j+1)     &
             &   +ch(5,i,j)*he(i,j-1)     &
             &   +ch(6,i,j)*he(ip1,j+1)   &
             &   +ch(7,i,j)*he(im1,j+1)   &
             &   +ch(8,i,j)*he(ip1,j-1)   &
             &   +ch(9,i,j)*he(im1,j-1)   &
             &   - rhs(i,j) )/(-ch(1,i,j))
         he(i,j)= (1.d0 -rlx)*he(i,j) +rlx*hstar
       enddo ! i
     enddo ! j
!$OMP END PARALLEL DO
   enddo ! ci
   enddo ! cj


! Boundary condition :                                              
! Since the old bc dh/dn is done perfectly, d(herror)/dn = 0        

  do i=1,NI 
    he(i,0) = he(i,1) 
    he(i,NJ+1) = he(i,NJ) 
  end do 

! Periodic bcs                                                      
  do j=0,NJ+1 
    he(0,j)= he(NI,j) 
    he(NI+1,j)= he(1,j) 
  end do 
                                                                        
  maxres= 0.d0 

!$OMP PARALLEL DO PRIVATE(j,i,res) SHARED(rhs,ch,he) REDUCTION(max:maxres) 
  do j=1,NJ 
    i=1 
    res=   rhs(i,j)-             &
        & (ch(1,i,j)*he(i,j)     &
        & +ch(2,i,j)*he(i+1,j)   &
        & +ch(3,i,j)*he(NI,j)    &
        & +ch(4,i,j)*he(i,j+1)   &
        & +ch(5,i,j)*he(i,j-1)   &
        & +ch(6,i,j)*he(i+1,j+1) &
        & +ch(7,i,j)*he(NI,j+1)  &
        & +ch(8,i,j)*he(i+1,j-1) &
        & +ch(9,i,j)*he(NI,j-1)) 
    if (abs(res).gt.maxres) then 
      maxres= abs(res) 
    end if 

    do i=2,NI-1  !     res = b - A x'                                                    
      res=   rhs(i,j)-               &
          & (ch(1,i,j)*he(i,j)       &
          & +ch(2,i,j)*he(i+1,j)     &
          & +ch(3,i,j)*he(i-1,j)     &
          & +ch(4,i,j)*he(i,j+1)     &
          & +ch(5,i,j)*he(i,j-1)     &
          & +ch(6,i,j)*he(i+1,j+1)   &
          & +ch(7,i,j)*he(i-1,j+1)   &
          & +ch(8,i,j)*he(i+1,j-1)   &
          & +ch(9,i,j)*he(i-1,j-1))  
      if (abs(res).gt.maxres) then 
        maxres= abs(res) 
      end if 
    enddo !i  

    i=NI 
    res=   rhs(i,j)-              &
        & (ch(1,i,j)*he(i,j)      &
        & +ch(2,i,j)*he(1,j)      &
        & +ch(3,i,j)*he(i-1,j)    &
        & +ch(4,i,j)*he(i,j+1)    &
        & +ch(5,i,j)*he(i,j-1)    &
        & +ch(6,i,j)*he(1,j+1)    &
        & +ch(7,i,j)*he(i-1,j+1)  &
        & +ch(8,i,j)*he(1,j-1)    &
        & +ch(9,i,j)*he(i-1,j-1)) 
    if (abs(res).gt.maxres) then 
       maxres= abs(res) 
    end if 
  enddo !j 
!$OMP END PARALLEL DO

   if (maxres.gt.3000.d0) then 
     print*, 'mprove.f90: STOP. res too large, i,j,maxres=',i,j,maxres                                             
     stop 
   end if 
   if (maxres.lt.tol) goto 13 
 enddo ! iter 
 13 continue 
                                                                       
 do j=0,NJ+1 
   do i=0,NI+1 
     h(i,j)= h(i,j) - he(i,j) 
   end do 
 end do 

return 
END                                           
