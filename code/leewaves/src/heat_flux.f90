!subroutine heat_flux(Tdif,step)
!subroutine heat_flux(Kdfluxdzt,Kdfluxdzb,sw_source,step)
subroutine heat_flux(step)

  !     ---------------------------------------------                     
  USE header

  !     use level m                                                       
  !     computes d2s/dz2 at the cell centers.                             

  implicit none 

  INTEGER :: i,j,k,step
  REAL(kind=rc_kind),  PARAMETER :: Cp= 4187.0
  REAL(kind=rc_kind) :: swrtemp, sw_source
  REAL(kind=rc_kind) :: fac, Kdfluxdzt, Kdfluxdzb, Kdfluxdzz(0:NK), swrd(NK)
  REAL(kind=rc_kind) :: Tdif(NI,NJ,NK)
  REAL :: wt1,wt2
! ----------------------------------

! This part deals with the water type coefficients: Paulson and Simpson
!  J_lambda1 = 0.6d0     
!  J_lambda2 = 20.d0     
!  J_A       = 0.62d0      
!Melissa's WINTER -SOUTH Bay  - had a mistake
!  J_lambda1 = 3.5d0     
!  J_lambda2 = 13.5d0     
!  J_A       = 0.52d0      
!Melissa's WINTER -NORTH Bay  - corrected April 2016
  J_lambda1 = 0.7d0     
  J_lambda2 = 16.d0     
  J_A       = 0.4d0      
  
! ----------------------------------

! This part deals with short wave radiation and heat loss (which includes longwave, sensible and latent heat)
  swr = 0.d0
  qloss = 0.d0




222 return

END
