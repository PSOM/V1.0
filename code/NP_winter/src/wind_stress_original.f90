subroutine wind_stress(udif,vdif,step) 
  !     ---------------------------------------------                     
  USE header

  implicit none 
  integer i,j,k,m,step,counter,io
  logical :: exist
  REAL(kind=rc_kind) :: udif(NI,NJ,NK),vdif(NI,NJ,NK)
  REAL(kind=rc_kind) :: stressxTS(1:73201),stressyTS(1:73201),stressprofile(NJ)
  REAL(kind=rc_kind) :: ycenter,ywindmin,ywindmax,edge
  REAL(kind=rc_kind) :: Kdudzt,Kdvdzt,fact,fac,rhoinv
  real(kind=rc_kind) :: yw , yset, y0
  integer            :: iyw, iy0
  real :: stressmax

  udif=0.;vdif=0.;

!*********************************
! COMPUTATION OF THE WIND STRESS

  ! IF NO WIND STRESS , then return with zero values
    stress_top_x=0.;stress_top_y=0;! stress_top= 0
!    stressmax= 1.0d0


    ! Open the wind stress forcing file
    open(unit=17, file='./input/ini_windstress.in')
    ! Import the wind stress time series for model forcing
    counter = 1
    do 
        read(17,*,IOSTAT=io) stressxTS(counter),stressyTS(counter)
        if (io > 0) then
            WRITE(*,*) 'Check wind stress input.  Something was wrong'
            STOP
            ! Reach end-of-file
        else if (io < 0) then
            PRINT*,"Read wind stress"
            EXIT
        else
            counter = counter + 1
        end if
    end do
    close (17)
      
    ! Ramping up heat flux (0 from day 0 to 5, linear from day 5 to 10, full after day 10)
    if (step.ge.10*400) then

    else if (step.le.5*400) then
        stressxTS(step) = 0.d0
        stressyTS(step) = 0.d0
    else 
        stressxTS(step) = stressxTS(step)*(step-(5*400))/(5.d0*400)
        stressyTS(step) = stressyTS(step)*(step-(5*400))/(5.d0*400)
    end if

    ! No meridional wind stress
    stressyTS(step) = 0.d0
    
    
    ! Print the wind stresss into a file windstress.out (every 1 day)
!    if (mod(step,100).eq.0) then
!    
!        inquire(file=TRIM(dirout)//'windstress.out', exist=exist)
!        if (exist) then
!            open(12, file=TRIM(dirout)//'windstress.out', status="old", position="append", action="write")
!        else
!            open(12, file=TRIM(dirout)//'windstress.out', status="new", action="write")
!            WRITE(12,*) "step, ","zonal windstress, ", "meridional windstress"
!        end if
!        WRITE(12,*) step, ",", stressxTS(step), ",", stressyTS(step)
!        close(12)
!    end if
        
    ! Apply a tanh profile in the meridional direction
    ycenter = 0.5*(yc(NJ)+yc(1)) 
    edge = 0.03    ! tightnness of the padding in the wind stress
    ywindmin = abs(atanh(2.d-2 - 1)/edge/PI) ! tapering region from southern bourndary
    ywindmax = yc(NJ)-ywindmin ! tapering region from northern bourndary

    do j=1,NJ
       if (yc(j).lt.yc(NJ/2)) then
           stressprofile(j)= 0.5*(tanh(edge*(yc(j)-ywindmin)*PI)+1.d0)
       else
           stressprofile(j)= -0.5*(tanh(edge*(yc(j)-ywindmax)*PI)-1.d0)
       end if
    end do

    do j=1,NJ
        do i=1,NI
            stress_top_x(i,j)=stressxTS(step)*stressprofile(j)
            stress_top_y(i,j)=stressyTS(step)*stressprofile(j)
        end do
    end do

  !************************************
! COMPUTATION OF THE SOURCE TERM

  fac= 1.d0/(UL*DL*delta) 
  fact = DL/UL 

  do j=1,NJ 
    do i=1,NI 
      stress_top(i,j) = sqrt(stress_top_x(i,j)*stress_top_x(i,j)+ stress_top_y(i,j)*stress_top_y(i,j))
      rhoinv = 1.d0/rho(i,j,NK) 
      !rhoinv = 1.d0/R0 
      Kdudzt= stress_top_x(i,j)*rhoinv*fact 
      Kdvdzt= stress_top_y(i,j)*rhoinv*fact 

      udif(i,j,NK)= fac*Jac(i,j,NK)*wz(i,j,NK)*Kdudzt   
      vdif(i,j,NK)= fac*Jac(i,j,NK)*wz(i,j,NK)*Kdvdzt   
    
    end do ! i
  end do ! j
                                                                        
return 
END                                           
                                                                        
                                                                        
