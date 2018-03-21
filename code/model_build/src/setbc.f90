subroutine setbc(step) 
!     ------------------------------------------------                  
  USE header
!     sets the initial face velocities  and fills the boundary arrays.  
!                                                                       
      integer i,j,k,step 
      REAL(kind=rc_kind) :: fac,utemp 

      do k=1,NK 
         do i=1,NI 
!     increase by a factor of 4 to speed spin up                        
!              vsouth(i,k)= 1.d-3                                       
            vsouth(i,k)= 0.d0 
            vfbcs(i,k)= vsouth(i,k)*vy(i,1)*Jjfc(i,0,k) 
         end do 
      end do 

!c     assign values to entrance values of s and T :arrays sent & Tent  
      do k=1,NK 
         do i=1,NI 
            ssouth(i,k)= 33d0 - S0 
!=            ssouth(i,k)= 5.7d0 - S0                                   
!     ssouth is used to specify bc in sTbc_periodicew and advecn        
         end do 
      end do 
!-         sbackgrnd= svert(np) - S0                                    
!=      sbackgrnd= 33.50 - S0                                           
      sbackgrnd= 34.0 - S0 
!=      sbackgrnd= 0.d0                                                 
      ! write(6,*)  'background s on shelf ', sbackgrnd +S0 
      return 
      END                                           
