PROGRAM main

  USE header
  use netcdf

#include "cppdefs.h"
#include "netcdf.inc"

  !Variables that will be read from the strain file.
  real(kind=4) ::  dudx(NI,NJ,NK),dudy(NI,NJ,NK),dudz(NI,NJ,NK), &
                    dvdx(NI,NJ,NK),dvdy(NI,NJ,NK),dvdz(NI,NJ,NK), &
                    drdx(NI,NJ,NK),drdy(NI,NJ,NK)

  !Variables that will be used for the computation of the Vertical Velocity (wgeo)
  real(kind=rc_kind), dimension(0:NI+1, 0:NJ+1, 0:NK+1)  :: divQ,wqgeo,Q1,Q2
  real(kind=rc_kind) ::  rcode,Q1dx,Q1dy,Q2dx,Q2dy

  ! Loop indexes, and omega-equation options.
  integer :: f
  integer, parameter :: OMEGA_FROM_TOTAL_VELOCITY=1, OMEGA_FROM_HYDROGRAPHY=0

  ! Some additional declarations
#include "main_declarations.h"
#include "ini_param.h"


   CALL ini_setup(pcorr)

   CALL meanh(NI,NJ,h,hmean)

   CALL sigma           ! calculates metric terms in vertical for the moving part of the grid
   CALL staticsigma     ! calculates metric terms in vertical for the fixed part of the grid
                        ! It is important to call sigma before staticsigma and ini_st
                        ! sigma needs to be called in momentum after each time advancement

  !  The initialization if over.
  !----------------------------------------------------------------------------------------

  ! Main loop:

  ! The user7 and user8 define how many times the vertical velocity should be computed
  ! the user6 decide which type of the equation should be used

  do f=user7,user8

    !-----(1) Read the velocities and the buyoancy gradients from the strain file----
    call read_cdf_strain(f*out3d_int,dudx,dudy,dudz,dvdx,dvdy,dvdz,drdx,drdy)

    !-----(2) Compute the two component of the Q-vector----
    do j=1,NJ
      do i=1,NI
        do k=1,NK
          if(user6.eq.OMEGA_FROM_HYDROGRAPHY) then
            Q1(i,j,k)= FPAR*ffc(i,j)*( dvdx(i,j,k)*dudz(i,j,k) +dvdy(i,j,k)*dvdz(i,j,k) )
            Q2(i,j,k)= -FPAR*ffc(i,j)*( dudx(i,j,k)*dudz(i,j,k) +dudy(i,j,k)*dvdz(i,j,k) )
          elseif(user6.eq.OMEGA_FROM_TOTAL_VELOCITY) then
            Q1(i,j,k)= ( dudx(i,j,k)*drdx(i,j,k) +dvdx(i,j,k)*drdy(i,j,k) )
            Q2(i,j,k)= ( dudy(i,j,k)*drdx(i,j,k) +dvdy(i,j,k)*drdy(i,j,k) )
          endif
        end do
      end do
    end do

    !-----(3) Final settings of Q-vector components----

      !(3a) boundary points
      do i=1,NI
        do k=1,NK
          Q1(i,0,k)= 2.0*Q1(i,1,k) -Q1(i,2,k)
          Q2(i,0,k)= 2.0*Q2(i,1,k) -Q2(i,2,k)
          Q1(i,NJ+1,k)= 2.0*Q1(i,NJ,k) -Q1(i,NJ-1,k)
          Q2(i,NJ+1,k)= 2.0*Q2(i,NJ,k) -Q2(i,NJ-1,k)
        end do
      end do

      do j=0,NJ+1
        do i=1,NI
          Q1(i,j,0)= 2.0*Q1(i,j,1) -Q1(i,j,2)
          Q2(i,j,0)= 2.0*Q2(i,j,1) -Q2(i,j,2)
          Q1(i,j,NK+1)= 2.0*Q1(i,j,NK) -Q1(i,j,NK-1)
          Q2(i,j,NK+1)= 2.0*Q2(i,j,NK) -Q2(i,j,NK-1)
        end do
      end do

      ! (3b) periodic-ew boundaries
      do k=0,NK+1
        do j=0,NJ+1
          Q1(0,j,k)= Q1(NI,j,k)
          Q2(0,j,k)= Q2(NI,j,k)
          Q1(NI+1,j,k)= Q1(1,j,k)
          Q2(NI+1,j,k)= Q2(1,j,k)
        end do
      end do

      !-----(4) Compute the divergence of the Q-vector ----
      do j=1,NJ
        do i=1,NI
          do k=1,NK
            Q1dx = 0.5*((Q1(i+1,j,k)-Q1(i-1,j,k))*ux(i,j) + &
                       (Q1(i,j+1,k) -Q1(i,j-1,k))*vx(i,j) + &
                       (Q1(i,j,k+1) -Q1(i,j,k-1))*wx(i,j,k) )
            Q2dx = 0.5*((Q2(i+1,j,k)-Q2(i-1,j,k))*ux(i,j) + &
                       (Q2(i,j+1,k) -Q2(i,j-1,k))*vx(i,j) + &
                       (Q2(i,j,k+1) -Q2(i,j,k-1))*wx(i,j,k) )
            Q1dy = 0.5*((Q1(i+1,j,k)-Q1(i-1,j,k))*uy(i,j) + &
                       (Q1(i,j+1,k) -Q1(i,j-1,k))*vy(i,j) + &
                       (Q1(i,j,k+1) -Q1(i,j,k-1))*wy(i,j,k) )
            Q2dy = 0.5*((Q2(i+1,j,k)-Q2(i-1,j,k))*uy(i,j) + &
                       (Q2(i,j+1,k) -Q2(i,j-1,k))*vy(i,j) + &
                       (Q2(i,j,k+1) -Q2(i,j,k-1))*wy(i,j,k) )
            divQ(i,j,k) = (Q1dx + Q2dy)/LEN
          end do
        end do
      end do

    !-----(5) Solve the omega equation----
    CALL psolve(wqgeo,divQ)

    !-----(6) Write the resulting vertical velocities----
    call write_cdf_omega(f*out3d_int,divQ,wqgeo)

  end do

END PROGRAM main
