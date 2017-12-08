subroutine evalrho_rho(rhonew,n) 
!     ---------------------------------------------                     
  USE header
!     s is the pot. density, so it is just copied into rho              
      REAL(kind=rc_kind) :: rhonew(0:NI+1,0:NJ+1,0:NK+1) 
      integer i,j,k,n 
      do k=0,NK+1 
         do j=0,NJ+1 
            do i=0,NI+1 
               rhonew(i,j,k)= s(i,j,k,n) 
            end do 
         end do 
      end do 
      return 
      END                                           
