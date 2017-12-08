subroutine ini_st

  USE header

!  initializes TS profiles from the Gulf Stream
!  regrided for the model
  implicit none
  integer n,i,j,k

! maximum and minimum potential density values
  real(kind=rc_kind), parameter :: rhoS=1025.8d0

! front parameters
!    + Af :length of the front in the y and z directions
!    + PI: inflection points position of the tanh functions:
!      + yPI:horizontal position of the center of the front.
!      + zPI:same as yPI. Allows to vary the mld along y
!    + mld:mixed layer depth

  real(kind=rc_kind) :: Afy, mld, Afz 
  real(kind=rc_kind) :: Arho, yPI, zPI, z, Fz, Fy(NJ)

  yPI=0.5*(yc(NJ+1)+yc(0))

  ! SUBMESOSCALE front

  Afy=15.d0/4.d0 
  Afy=40.d0/4.d0 
  mld=50.d0
  mld=50.d0
  Afz=3.d2/4.d0

  zPI=mld+Afz

  do j=1, NJ
    !Fy(j)=(tanh( (yc(j)-yPI)/Afy  ))*1.25
    Fy(j)=(tanh( (yc(j)-yPI)/Afy  ))/(4.d0/3.d0)
    !print*, j, Fy(j)
  end do

  do k=0, NK+1
  do j=1, NJ
    z= DL*zc(1,j,k)
    !Fz=(tanh( (zPI*Fy(j)-z)/Afz )+1.d0)/4.
    Fz=(tanh( (zPI*Fy(j)-z)/Afz )+1.d0)/2.
    !s(0,j,k,0)=s(0,j,k,0)+Fz
    s(0,j,k,0)=rhoS+Fz
    !if (j==1 .or. j==NJ) print*, k, z, Fz
  end do
  end do

  do k=1, NK
    s(0,0   ,k,0)=s(0,1 ,k,0)
    s(0,NJ+1,k,0)=s(0,NJ,k,0)
!    do j=350,NJ+1
!      s(0,j,k,0)=s(0,349,k,0)
!    enddo
  end do

  do k=0, NK+1
  do j=0, NJ+1
  do i=0, NI+1
    s(i,j,k,0)= s(0,j,k,0)
  end do
  end do
  end do

  ! reference density values for sponge layers
  do k=0, NK+1
    rho_refS(k)=s(0, 1,k,0)
    rho_refN(k)=s(0,NJ,k,0)
  enddo
  do j=0, NJ+1
    rho_refB(j)=s(0,j,0,0)
  enddo

  return

END


  
