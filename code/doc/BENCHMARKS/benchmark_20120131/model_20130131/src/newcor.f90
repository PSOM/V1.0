subroutine newcor(dtimel,n) 
!----------------------------------------------------                   
  USE header
!     modified for periodicew bcs                                       
!     Evaluate the Coriolis terms si,sj                                 
   integer i,j,k,n 
   REAL(kind=rc_kind) :: fac2,ainv,dtimel,con 
!     We are using the values at the ghost points.                      
  fac2= EPS*lambda 
  ainv= 1.d0/apr 
                                                                    
  do k=1,NK 
    do j=1,NJ 
      do i=1,NI 
        si(i,j,k)= 0.d0 
        sj(i,j,k)= 0.d0 
        sk(i,j,k)= 0.d0 
      enddo
    enddo
  enddo
                                                                       
 return 
 END                                           
