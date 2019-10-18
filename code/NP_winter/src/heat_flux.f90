subroutine heat_flux(Tdif,step,n)
    !     ---------------------------------------------                     
    USE header

    !     use level m                                                       
    !     computes d2s/dz2 at the cell centers.                             

    implicit none 

    INTEGER :: i,j,k,step,n,counter,io
    character(len=5) :: x1
    REAL(kind=rc_kind) :: fac, Kdfluxdzt, Kdfluxdzb, Kdfluxdzz(NK), swrd(NK)
    REAL(kind=rc_kind) :: Tdif(NI,NJ,NK), ycenter, qlossTS(NJ+2,nsteps), myalpha, temp
    ! Entire time series of Heat flux must be imported to allow pick up. 
    ! The do-loop therefore includes a test for end-of-file
    REAL(kind=rc_kind) :: swrTS(0:56010)
    logical :: exist
    ! ----------------------------------

    ! This part deals with the water type coefficients
    J_lambda1 = 0.6d0
    J_lambda2 = 20.d0
    J_A       = 0.62d0
    ! ----------------------------------

    ! This part deals with short wave radiation and heat loss (which includes longwave, sensible and latent heat)

    ! Set the heat fluxes to ZERO
    ! swr(:) = 0.d0
    ! qloss(:) = 0.d0
    
    ! Open the HF forcing file
    open(unit=1900, file='./input/ini_heatfluxes_108s.in')
    ! Import the heat fluxes time series for model forcing
    counter = 1
    do 
        read(1900,fmt=*,IOSTAT=io) temp
        if (io > 0) then
            WRITE(*,*) 'Check HF input.  Something was wrong'
            STOP
        else if (io < 0) then
            PRINT*,"Read Heat Flux"
            EXIT
        else
            swrTS(counter) = temp
            counter = counter + 1
        end if
    end do
    close (1900)
        
    ! Apply a meridional gradient based on NARR data.
    ! gradient of 1/24 W/m2/km
     ycenter = 0.5*(yc(NJ)+yc(1))

    do j = 0,NJ+1
        ! Compute the meridional HF, accounting for the Runge-Kutta scheme (hence the *3-2 term)
        swr(j)= (-15.0d0/NJ*(yc(j)-ycenter)+ swrTS(step-8000))
        !swr(j) = swrTS(j+1)*myalpha
        !qloss(j)= (12.0d0/NJ*(yc(j)-ycenter)+ qlossTS(step))*myalpha
        qloss = 0.d0
    end do
        
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
