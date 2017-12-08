subroutine inith 
  !     ------------------------------------------------                  
  USE header
  ! Initialize h using thermal wind balance
  ! Any density distribution is fine. 
  implicit none
  integer i,j,k,jlim
  REAL(kind=rc_kind) :: fu(0:NJ),const,hsum,dz
  REAL :: potdens,coef

  REAL :: zl, sigup
  INTEGER :: sig
  INTEGER :: optionl

  ! In case the bottom if flat, h is defined in such a way that it exactly cancels speed at the bottom, 
  !                             supposing thermal wind balace.
  ! otherwise, velocity will be canceled at a given geopotential depth.


  if(lv_flat_bottom) then
    ! Two strategies are possible.
    ! optionl=1 works for solid boundaries in y=0 and y=ymax
    ! optionl=2 is based on the computation of the baroclinic pressure term. It cancels velocity at the CENTER
    !                    of the lowermost cell.
    optionl=1

    if(optionl==1) then
      const=DL/R0 !fu is actually fu / g and is at cell faces
      do i=1,NI
        do j=1,NJ-1
          fu(j)= 0.d0
          do k=1,NK
            ! findzall has been called, so zf can be used. sigma not yet called, so wz not usable
            dz   = zf(i,j,k)-zf(i,j,k-1)
            drho = (rho(i,j+1,k)- rho(i,j,k))*vy(i,j)/LEN
            fu(j)= fu(j) + const*drho*dz  !fu is at cell faces,is actually fu/g
          end do
        end do
        fu(0)= fu(1)
        fu(NJ)= fu(NJ-1)
        ! at k=NK, fu = g*hy
        h(i,0)= 0.d0
        do j=1,NJ+1
          h(i,j)= h(i,j-1) - (LEN/vy(i,j))*fu(j-1)
        end do
      end do
    endif

    if(optionl==2) then    
      do i=1,NI+1
        h(i,0)= 0.d0
        do j=1,NJ
          h(i,j)= h(i,j-1) - HL/(gpr*gj(i,j-1,1,2))*grpjfc(i,j-1,1)
        end do
      end do
    endif

 else

  ! In the general case of non-flat bottom, h cannot cancel bottom speed.
  ! Here, it simply cancels speed at a given depth, zl, that has to be within the domain at every (x,y) location.
  ! It only works for solid boundaries in y=0 and y=ymax and in the case where h is flat at y=0.

  CALL staticsigma
  CALL rpevalgrad_Sche(0)

  optionl=1

  if(optionl==1) then
    zl=-97.;
    do i=1,NI
      h(i,0)= 0.d0; h(i,1)= 0.d0;
      do j=2,NJ
        CALL findsigma(i,j-1,zl,sig,sigup);! sig=1;sigup=1. 
        h(i,j)=h(i,j-2)-1./(0.5*gpr*vy(i,j-1))*(rv4_Sche(i,j-1,sig)*sigup+rv4_Sche(i,j-1,sig-1)*(1.-sigup))
      end do
    end do
   else 
    PRINT*,"optionl=1 is the only option in inith!"
  endif

 endif 
  
  !------------------------------------------------

 ! Additionally, a barotropic jet can be added:
!jlim=118.
!do i=0,NI+1
!  do j=1,NJ
!    coef=0-(tanh((REAL((j-jlim))*5.)/NJ*PI)-1.)/2.   
!    h(i,j)=h(i,j)+1.*(coef-1.)*0.05/0.24*0.1
!  enddo
!enddo




  ! periodicity
  h(:,NJ+1)=h(:,NJ)
  do j=0,NJ+1
    h(0,j)= h(NI,j)
    h(NI+1,j)= h(1,j)
  end do

  ! h should be zero-average.
  hsum=0.d0 
  do i=1,NI 
    do j=1,NJ 
      hsum= hsum +h(i,j)
    end do
  end do
  hmean= hsum/dble(NI*NJ) 
  h=h-hmean
  h=h/HL  

 ! h=0.

  return 
END subroutine inith
