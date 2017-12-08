subroutine initial_tke()
! This initializes turbulence routines in GOTM

!!! Access relevant subroutines from init_tke
use header
use turbulence,           only : eps1d => eps
use turbulence,           only : tke, L, num, nuh, nus
use turbulence,           only : init_turbulence, do_turbulence, clean_turbulence
use mtridiagonal,         only: init_tridiagonal

implicit none

!!! Some local variables
integer                 :: i, j, k, kdp

REAL(kind=rc_kind), dimension(    1:NI  ,1:NJ  , 1:NK  ) :: dum1
REAL(kind=rc_kind), dimension(    1:NI  ,1:NJ  , 1:NK  ) :: dum2

!!!!!!!!!!! Initialize Turbulence Modules !!!!!!!!!!!!!!
call init_turbulence(10, 'gotmturb.nml', NK)
call init_tridiagonal(NK)
write(*,*) 'Turbulence Modules initialized'

call wind_stress(dum1,dum2,1)

do i=1,NI
  do j=1,NJ
    u_fric(i,j) = (stress_top(i,j)/1000.d0)**0.5
  end do
end do

!===== Initialize lengthscale profile =============
do i=1,NI
   do j=1,NJ

!===========
     do k=NK-4, NK
       psom_l(i,j,k) = abs(zf(i,j,k)-zf(i,j,k-1))*1000;
     end do

!      psom_l(i,j,NK) = 1.4794d-2
!      write(6,*) psom_l(i,j,NK)
      do k=0,NK-5
!      psom_l(i,j,k) = 1.45918E-2
!      psom_l(i,j,k) = 1.45918E-1
       psom_l(i,j,k) = 1.45918E00
      end do
!===========
if (i.eq.NI/2 .AND. j.eq.NJ/2) then
write(6,*) 'Length-scale values at near-surface'
do k=NK-10,NK
  write(6,*) psom_l(i,j,k)
end do
end if

   end do
end do
!==================================================


!======= Initialize tke, eps based on lengthscale
do i=0,NI+1
   do j=0,NJ+1
      psom_tke(i,j,NK) = (u_fric(i,j)/cmu0_psom)**2
!      write(6,*) psom_tke(i,j,NK)
       psom_eps(i,j,NK) = (psom_tke(i,j,NK)**1.5)*(2**1.5)/(B1_MY*psom_l(i,j,NK))
!       psom_eps(i,j,NK) = (psom_tke(i,j,NK)**1.5)*(2**1.5)/(psom_l(i,j,NK))
!      write(6,*) psom_eps(i,j,NK)

!===== This part is for finding mixed-layer=========
!===== not necessary ===============================
!do k=0,NK
!if ((zc(i,j,k)*(-1000)).le.mldepth) then
!  ! break out of the loop
!  kdp = k+1
!  exit
!else
!  continue
!end if
!end do
!! Now you have kdp. Set a linear tke profile
!do k=0,NK-1
!if (k.le.kdp) then
!  psom_tke(i,j,k) = 0
!else
!  psom_tke(i,j,k) = (psom_tke(i,j,NK)/mldepth)*(mldepth + (zc(i,j,k)*1000))
!  psom_eps(i,j,k) = (psom_tke(i,j,k)**1.5)*(2**1.5)/(B1_MY*psom_l(i,j,k))
!end if
!end do
!===================================================

      do k = 0,NK-1
         psom_tke(i,j,k) = 1.d-10
!! Change the dissipation rate instead of this
!         psom_eps(i,j,k) = 1.d-14
         psom_eps(i,j,k) = (psom_tke(i,j,k)**1.5)*(2**1.5)/(B1_MY*psom_l(i,j,k))
      end do
!===================================================

   end do
end do

write(6,*) NI, NJ, NK

write(6,*) 'Wind stress at mid point is', stx(NI/2,NJ/2)

write(6,*) 'Friction velocity at mid point is', u_fric(NI/2,NJ/2)

write(6,*) 'Turbulent Kinetic Energy at mid point is', psom_tke(NI/2, NJ/2, NK)

write(6,*) 'Turbulent Kinetic Energy below the surface', psom_tke(NI/2, NJ/2, NK-2)

write(6,*) 'Macro Lengthscale', psom_l(NI/2, NJ/2, NK-1)

write(6,*) 'Macro Lengthscale below the surface', psom_l(NI/2, NJ/2, NK-4)

write(6,*) 'Dissipation Rate of TKEat mid point', psom_eps(NI/2, NJ/2, NK)

write(6,*) 'Dissipation Rate of TKE below the surface', psom_eps(NI/2, NJ/2, NK-4)


return
end subroutine initial_tke


