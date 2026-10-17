subroutine hsolve(h,oldh,hdt,dtime) 
!     --------------                                                    
  USE header, only : NI,NJ,rc_kind
!     modified for periodicew boundaries                                
!     ----------------------------------                                
      integer iter,i,j,kount,l,im1,ip1,ci,cj,nci,ib,ie,is 
      REAL(kind=rc_kind) :: ch(9,NI,NJ),rhs(NI,NJ),                          &
     &     oldh(0:NI+1,0:NJ+1),h(0:NI+1,0:NJ+1),res,                    &
     &     hdt(0:NI+1,0:NJ+1),dti,dtime,tol,maxres,hstar,rlx            
!     &     dti,dtime,tol,maxres,hstar,rlx                              
                                                                        
  tol= 1.d-12 
  dti= 1.d0/dtime 
  kount=0 

!     compute fine grid coefficients and rhs on fine grid               
  call chfine(dtime,ch,rhs) 


  ! colours in i for the multi-colour SOR sweep (x periodic: odd NI needs a third colour)
  nci= 2
  if (mod(NI,2).ne.0) nci= 3
  if (NI.eq.1) nci= 1

  rlx= 1.72 
  do iter=1,100 

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
!$OMP PARALLEL DO PRIVATE(j,i,im1,ip1,hstar) SHARED(h,ch,rhs,rlx,ib,ie,is,cj)
      do j=cj,NJ,2
        do i=ib,ie,is
          im1= i-1
          ip1= i+1
          if (i.eq.1) im1= NI
          if (i.eq.NI) ip1= 1
          hstar= ( ch(2,i,j)*h(ip1,j)     &
              &   +ch(3,i,j)*h(im1,j)     &
              &   +ch(4,i,j)*h(i,j+1)     &
              &   +ch(5,i,j)*h(i,j-1)     &
              &   +ch(6,i,j)*h(ip1,j+1)   &
              &   +ch(7,i,j)*h(im1,j+1)   &
              &   +ch(8,i,j)*h(ip1,j-1)   &
              &   +ch(9,i,j)*h(im1,j-1)   &
              &   - rhs(i,j) )/(-ch(1,i,j))
          h(i,j)= (1.d0 -rlx)*h(i,j) +rlx*hstar
        enddo ! i
      enddo ! j
!$OMP END PARALLEL DO
    enddo ! ci
    enddo ! cj

    do l=1,3 
      call hfill(dtime,h) 
    enddo
                                                                          
    maxres= 0.d0 

!$OMP PARALLEL DO PRIVATE(j,i,res) REDUCTION(max:maxres) 
    do j=1,NJ 
      i=1 
      res= rhs(i,j)-              &
        & (ch(1,i,j)*h(i,j)       &
        & +ch(2,i,j)*h(i+1,j)     &
        & +ch(3,i,j)*h(NI,j)      &
        & +ch(4,i,j)*h(i,j+1)     &
        & +ch(5,i,j)*h(i,j-1)     &
        & +ch(6,i,j)*h(i+1,j+1)   &
        & +ch(7,i,j)*h(NI,j+1)    &
        & +ch(8,i,j)*h(i+1,j-1)   &
        & +ch(9,i,j)*h(NI,j-1))   
      if (abs(res).gt.maxres) then 
        maxres= abs(res) 
      end if 
      do i=2,NI-1    !   res = b - A x'                                                    
        res= rhs(i,j)-              &
          & (ch(1,i,j)*h(i,j)       &
          & +ch(2,i,j)*h(i+1,j)     &
          & +ch(3,i,j)*h(i-1,j)     &
          & +ch(4,i,j)*h(i,j+1)     &
          & +ch(5,i,j)*h(i,j-1)     &
          & +ch(6,i,j)*h(i+1,j+1)   &
          & +ch(7,i,j)*h(i-1,j+1)   &
          & +ch(8,i,j)*h(i+1,j-1)   &
          & +ch(9,i,j)*h(i-1,j-1))  
        if (abs(res).gt.maxres) then 
          maxres= abs(res) 
        end if 
      enddo !i
      i=NI 
      res= rhs(i,j)-              &
        & (ch(1,i,j)*h(i,j)       &
        & +ch(2,i,j)*h(1,j)       &
        & +ch(3,i,j)*h(i-1,j)     &
        & +ch(4,i,j)*h(i,j+1)     &
        & +ch(5,i,j)*h(i,j-1)     &
        & +ch(6,i,j)*h(1,j+1)     &
        & +ch(7,i,j)*h(i-1,j+1)   &
        & +ch(8,i,j)*h(1,j-1)     &
        & +ch(9,i,j)*h(i-1,j-1))  
      if (abs(res).gt.maxres) then 
         maxres= abs(res) 
      end if 
    enddo !j
!$OMP END PARALLEL DO

    !  write(6,*) 'iter,maxres',iter,maxres                          
    if (maxres.gt.3000.d0) then 
      print*, 'hsolve.f90: STOP. res too large, i,j,maxres=',i,j,maxres                                             
      stop 
    end if 
    if (maxres.lt.tol) goto 13 
  enddo !iter 
 
13 continue
                                                                       
  do l=1,1 
     call mprove(h,ch,rhs,dtime) 
  end do 
                                                                        
return 

END                                           
