subroutine mprove(h,ch,oldrhs,dtime) 
!     --------------                                                    
  USE header, only : NI,NJ,rc_kind
      integer iter,i,j,l 
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
                                                                        
 rlx= 1.72 
 do iter=1,1000 

   do j=1,NJ 
     i=1 
     hstar= ( ch(2,i,j)*he(i+1,j)     &
           & +ch(3,i,j)*he(NI,j)      &
           & +ch(4,i,j)*he(i,j+1)     &
           & +ch(5,i,j)*he(i,j-1)     &
           & +ch(6,i,j)*he(i+1,j+1)   &
           & +ch(7,i,j)*he(NI,j+1)    &
           & +ch(8,i,j)*he(i+1,j-1)   &
           & +ch(9,i,j)*he(NI,j-1)    &
           & - rhs(i,j) )/(-ch(1,i,j))                              
     he(i,j)= (1.d0 -rlx)*he(i,j) +rlx*hstar 
     do i=2,NI-1 
       hstar= ( ch(2,i,j)*he(i+1,j)     &
             & +ch(3,i,j)*he(i-1,j)     &
             & +ch(4,i,j)*he(i,j+1)     &
             & +ch(5,i,j)*he(i,j-1)     &
             & +ch(6,i,j)*he(i+1,j+1)   &
             & +ch(7,i,j)*he(i-1,j+1)   &
             & +ch(8,i,j)*he(i+1,j-1)   &
             & +ch(9,i,j)*he(i-1,j-1)   &
             & - rhs(i,j) )/(-ch(1,i,j))                           
       he(i,j)= (1.d0 -rlx)*he(i,j) +rlx*hstar 
     enddo !i
     i=NI 
     hstar= ( ch(2,i,j)*he(1,j)      &
           & +ch(3,i,j)*he(i-1,j)    &
           & +ch(4,i,j)*he(i,j+1)    &
           & +ch(5,i,j)*he(i,j-1)    &
           & +ch(6,i,j)*he(1,j+1)    &
           & +ch(7,i,j)*he(i-1,j+1)  &
           & +ch(8,i,j)*he(1,j-1)    &
           & +ch(9,i,j)*he(i-1,j-1)  &
           & - rhs(i,j) )/(-ch(1,i,j))    
     he(i,j)= (1.d0 -rlx)*he(i,j) +rlx*hstar 
   enddo !j 


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
