subroutine ini_h 
  !     ------------------------------------------------                  
  USE header
  ! Initialize h using thermal wind balance
  ! Assume the level of no motion at z=-D  (bottom)
  ! Any density distribution is fine. 
  ! findzall has been called, so zf can be used. sigma not yet called, so wz 
  ! not usable
  implicit none
  integer i,j,k,jlim
  REAL(kind=rc_kind) :: fu(0:NJ),const,hsum,dz
  REAL :: potdens,coef

  REAL :: zl, sigup
  INTEGER :: sig
  INTEGER :: optionl


  ! in case the bottom if flat, h is defined in such a way that it exactly cancels speed at the bottom, 
  ! supposing thermal wind balace.
  if(lv_flat_bottom) then

    optionl=1

    if(optionl==1) then
      const=DL/R0
      !fu is actually fu / g and is at cell faces
      do i=1,NI
         do j=1,NJ-1
            fu(j)= 0.d0
            do k=1,NK
               !fu is at cell faces,is actually fu/g
               dz   = zf(i,j,k)-zf(i,j,k-1)
               drho = (rho(i,j+1,k)- rho(i,j,k))*vy(i,j)/LEN
               fu(j)= fu(j) + const*drho*dz
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
        ! at k=NK, fu = g*hy
        h(i,0)= 0.d0
        do j=1,NJ
          h(i,j)= h(i,j-1) - HL/(gpr*gj(i,j-1,1,2))*grpjfc(i,j-1,1)
        end do
      end do
    endif

 ! In the general case of non-flat bottom, h cannot cancel bottom speed.
 ! Here, it simply cancels speed at a given depth, zl, that should be within the domain at every (x,y) location.
 else

  CALL staticsigma
  CALL rpevalgrad_Sche(0)

  optionl=1

  if(optionl==1) then

  zl=-97.;
  do i=1,NI
    ! at k=NK, fu = g*hy
    h(i,0)= 0.d0
    h(i,1)= 0.d0
    do j=2,NJ
      CALL findsigma(i,j-1,zl,sig,sigup);! sig=1;sigup=1. 
      h(i,j)=h(i,j-2)-1./(0.5*gpr*vy(i,j-1))*(rv4_Sche(i,j-1,sig)*sigup+rv4_Sche(i,j-1,sig-1)*(1.-sigup))
      ! h(i,j)=h(i,j-2)-1./(0.5*gpr*vy(i,j-1))*(rv4_Sche(i,j-1,1))
    end do
  end do

! do j=1,NJ
!   do i=3,NI+1
!     CALL findsigma(i-1,j,zl,sig,sigup);! sig=1;sigup=1.
!     h(i,j)=h(i-2,j)-1./(0.5*gpr*ux(i-1,j))*(ru4_Sche(i-1,j,sig)*sigup+ru4_Sche(i-1,j,sig-1)*(1.-sigup))
!     ! h(i,j)=h(i-2,j)-1./(0.5*gpr*ux(i-1,j))*(ru4_Sche(i-1,j,1))
!    end do
! end do

  !h=0.


  else

  

  endif



  ! Addition of a barotropic jet.
  jlim=117.+user4
  do i=0,NI+1
    do j=1,NJ
      coef=0-(tanh((REAL((j-jlim))*3)/NJ*PI)-1.)/2.*user2 
!@
     h(i,j)=h(i,j)+(coef-1.)*0.05/0.24*0.1
    enddo
  enddo


 endif 


  ! periodicity
  h(:,NJ+1)=h(:,NJ)
  do j=0,NJ+1
     h(0,j)= h(NI,j)
     h(NI+1,j)= h(1,j)
  end do

  ! h should be zero-average.
  hsum=0.d0 
  do  i=1,NI 
     do  j=1,NJ 
        hsum= hsum +h(i,j)
     end do
  end do
  hmean= hsum/dble(NI*NJ) 
  h=h-hmean
  h=h/HL  

 ! h=0.


  return 
END subroutine
