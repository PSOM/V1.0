subroutine wind_stress(step) 
  !     ---------------------------------------------                     
  USE header
  ! updates  stress_top_x, stress_top_y and stress_top (which is the magnitude of the wind stress)

  integer step
!   implicit none
   
!   integer i,j,k,m,step                                                                      
!   REAL(kind=rc_kind) :: udif(NI,NJ,NK),vdif(NI,NJ,NK)
!   REAL(kind=rc_kind) :: Kdudzt,Kdvdzt,fact,fac,rhoinv
!   real(kind=rc_kind) :: U10,rhoair,CD,ycenter,ywindmin,ywindmax,edge,shapewind(NJ)
!   real :: stressmax_x,stressmax_y
!   REAL :: wt1,wt2
! ----------------------------------

 ! udif=0.;vdif=0.;


!  return

!*********************************
! COMPUTATION OF THE WIND STRESS

! IF NO WIND STRESS , then return with zero values
  stress_top_x=0.;stress_top_y=0; stress_top= 0
!   stressmax_x= 0.0d0
!   stressmax_y= 0.d0
!
! !  lenarry= 352 (in header_expe.h)
!   if (step.eq.1) then
!      rhoair= 1.2   ! kg/m3
!      ! uwind(lenarray), vwind(lenarray) already read in in heat_flux.f90.  daycount(lenarray) also defined.
!      do i=1,lenarray
!         ! daycount, stressx and stressy are saved as global arrays
!         ! Trenberth,Large and Olson (1989)  formula for wind stress
!         U10= sqrt( uwind(i)*uwind(i) +vwind(i)*vwind(i))
!         if (U10.le.3.0) then
!            CD= (0.62 + 1.56/U10)*1d-3
!         else if (U10.le.10.) then
!            CD = 1.14*1d-3
!         else
!            CD = (0.49 +0.065*U10)*1d-3
!         end if
!         stressx(i) = rhoair*CD*U10*uwind(i)
!         stressy(i) = rhoair*CD*U10*vwind(i)
!      end do
!      write(6,*) 'Calc wind stress'
!   end if
!
! ! Find nearest time for interpolation
!   do i=2,lenarray
!      if (daycount(i).gt.time_days) then
!         wt1=(daycount(i)-time_days)/(daycount(i)-daycount(i-1))
!         wt2= 1.-wt1
!         stressmax_x= wt2*stressx(i) + wt1*stressx(i-1)
!         stressmax_y= wt2*stressy(i) + wt1*stressy(i-1)
!         goto 222
!      end if
!   end do
!
! 222 continue
!   ycenter = 0.5*(yc(NJ)+yc(1))
! !  ywindmin= 10.0 ! 10 km from coast
! !  ywindmax= yc(NJ)-10.d0    ! 25 km from boundary
!
! !This caused strong w at the coast
! !  ywindmin= 3.0 ! 10 km from coast
! !  ywindmax= yc(NJ)-3.d0    ! 25 km from boundary
! !  ywindmin= 7.0 ! 10 km from coast   Nidhi sections
! !  ywindmax= yc(NJ)-7.d0    ! Nidhi sections
!   ywindmin= 20.0 ! 20 km from coast   Leg2 sections
!   ywindmax= yc(NJ)-20.d0    ! Leg2 sections
!
! ! tightnness
!   edge=0.06  ! this spreads it out  much more
! !  edge=0.15   ! for Nidhi6 sections
!
!
!   do j=1,NJ
!      if (yc(j).lt.yc(NJ/2)) then
!         shapewind(j)= 0.5*(tanh(edge*(yc(j)-ywindmin)*PI)+1.d0)
!      else
!         shapewind(j)= -0.5*(tanh(edge*(yc(j)-ywindmax)*PI)-1.d0)
!      end if
!   end do
!
!   do j=1,NJ
!      stress_top_x(1:NI,j)=shapewind(j)*stressmax_x
!      stress_top_y(1:NI,j)=shapewind(j)*stressmax_y
!   end do
!
!
!   do j=1,NJ
!      do i=1,NI
!         stress_top(i,j)= sqrt( stress_top_x(i,j)*stress_top_x(i,j) + stress_top_y(i,j)*stress_top_y(i,j) )
!      end do
!   end do

                                                                        
return 
END                                           
                                                                        
                                                                        
