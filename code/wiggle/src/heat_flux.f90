
subroutine heat_flux(Tdif,step)
  !     ---------------------------------------------                     
  USE header

  !     use level m                                                       
  !     computes d2s/dz2 at the cell centers.                             

  implicit none 

  INTEGER :: i,j,k,step
  REAL(kind=rc_kind) :: swrtemp
  REAL(kind=rc_kind) :: fac, Kdfluxdzt, Kdfluxdzb, Kdfluxdzz(NK), swrd(NK)
  REAL(kind=rc_kind) :: Tdif(NI,NJ,NK)

! ----------------------------------

! This part deals with the water type coefficients
  J_lambda1 = 0.6d0     
  J_lambda2 = 20.d0     
  J_A       = 0.62d0      
! ----------------------------------

! This part deals with short wave radiation
  swr = 0.d0
  qloss = 0.d0

! Add sinusoidal heat flux
  swrtemp = 0d0*sin( (2*3.14159d0/(24*3600))*step*dtf*(1.d05)    )
  if (swrtemp.lt.(0.d0)) then
    swrtemp = 0.d0
  endif

  !swr  (:) = swrtemp
  do j=1,NJ
    swr  (j) = 900.d0*sin(dble(j-1)/dble(NJ-1)*pi/2)
  enddo
  qloss(:) = 0d0/3.14159
  swr  (:) = 0.d0

! ----------------------------------

  fac= 1.d0/(UL*DL*delta)

  do j=1,NJ 
    do i=1,NI 

      Kdfluxdzt = DL*( swr(j) - qloss(j) )/(R0*4187.d0)
      do k=1,NK
        swrd(k) = swr(j)*( J_A*exp(zf(NI/2,NJ/2,k)*DL/J_lambda1) + (1 - J_A)*exp(zf(NI/2,NJ/2,k)*DL/J_lambda2)   )
        Kdfluxdzz(k) = DL*( swrd(k) )/(R0*4187.d0)
      end do
      Kdfluxdzb = 0.d0

      Tdif(i,j,NK) =fac*Jac(i,j,NK)*wz(i,j,NK)*Kdfluxdzt

      do k=2,NK-1
        Tdif(i,j,k)=fac*Jac(i,j, k)*wz(i,j, k)*(Kdfluxdzz(k) - Kdfluxdzz(k-1))
      end do
      Tdif(i,j, 1) =fac*Jac(i,j, 1)*wz(i,j, 1)*Kdfluxdzb

    enddo
  enddo

return

END
