!!! This subroutine interfaces with gotm and calculates diffusion coefficients
subroutine couple_gotm(I,J)
!! Downloads copy
!!! Subroutines accessed from within couple_gotm
USE header
use turbulence,           only : eps1d => eps
use turbulence,           only : tke, L, num, nuh, nus
!use turbulence,           only : P1d => P
!use turbulence,           only : B1d => B
use turbulence,           only : init_turbulence, do_turbulence, clean_turbulence
use mtridiagonal,         only : init_tridiagonal

! dtf (PSOM) => dt (GOTM) ->> The time step
! NK (PSOM) => nlev (GOTM) ->> Number of levels along vertical
!!! Velocity variables
!!! u, v, w declared at an array size (0:NI+1,0:NJ+1,0:NK+1,0:1)
!!! freqN2(i,j,k) is the buoyancy frequency (mymodules.f90)
!!! u_fric(I,J) is the surface friction velocity at I,J (mymodules.f90)
!!! Note - This variable is only used for setting initial tke, eps, L

implicit none
!!! Local variables
integer, intent(in)     ::  I,J
integer                 ::  iii, ic
integer                 ::  k, MaxItz0b 
REAL(kind=rc_kind)      ::  charnockval, rrb, clip 
!!! Molecular Viscosity, heat and salt diffusivities (used for gotm calculations)
REAL(kind=rc_kind)      :: avmolu , avmolt , avmols  
!!! Variables for interfacing with gotm
REAL(kind=rc_kind)      :: dpth, gotm_dt
!!! Bottom friction velocity, surface and bottom roughness parameters
REAL(kind=rc_kind)      :: u_taub, z0s_gotm, z0b_gotm
!!! Bottom roughness = h0b
!!!  Note: z0b=0.03*h0b+0.1*(molecular viscosity)/ustar,
!!!  Von-Karman constant, bottom friction coefficient,
REAL(kind=rc_kind)      :: h0b   

!-- Computation part -----------------------------
charnockval = 1400.d0; clip = 1.0E0; MaxItz0b = 1;
avmolu = 1.3d-6; avmolt = 1.4d-7; avmols = 1.1d-9;
h0b = 0.05;

!!! Calculate friction velocities, bottom and surface roughness lengths, shear and
!!! buoyancy frequencies, set sea grass production to 0E0


!!! Time Step for PSOM/GOTM is set in momentum equation
!!! according to Runge-Kutta-3rd order scheme...

!!! Calculating layer thickness - added to cell centres
DO k = 1, NK
   thk(k) = (zf(I,J,k) - zf(I,J,k-1))*DL
   lthk(I,J,k) = thk(k)
END DO
   thk(0) = thk(1)

!!! Depth
dpth = (zf(I,J,NK) - zf(I,J,0))*DL
!! Set a time step for gotm

!!! Buoyancy frequency freqN2(I,J,K) and shear squared
do k = 0, NK
   NN1d(k) = nnface(I,J,k)
   SS2(k)  = shq(I,J,k)
!!! Change the values to cell face, use subroutine nnssface
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SeaGrass Turbulent Kinetic Energy Production
TKEP(:) = 0E0

!!! Calculate surface roughness length. Surface friction velocity is known already.
!!! Initialize the wind stress profile
u_fric(I,J) = (stress_top(i,j)/1000.d0)**0.5

!IF (I.eq.(NI/2) .AND. J.eq.(NJ/2))
!  write(6,*) 'Validation- friction velocity'
!  write(6,*) u_fric(I,J), u_fric(I-1,J-1), u_fric(I+1,J+1)
!ENDIF

z0s_gotm = (u_fric(I,J)**2)*charnockval/9.81

!!!!!!!!!!!! Bottom Friction and Roughness !!!!!!!!!!!!!
u_taub = 6.7e-6  ! Limiting case

do iii=1,MaxItz0b
   z0b_gotm = 0.1*avmolu/max(avmolu,u_taub)+0.03*h0b
   
!compute the factor r (version 1, with log-law)   
   rrb=k_von/(log((z0b_gotm+thk(1)/2.d0)/z0b_gotm))

!  compute the factor r (version 2, with meanvalue log-law)
!  frac=(z0b+h(1))/z0b
!  rrb=kappa/((z0b+h(1))/h(1)*log(frac)-1.)

!  compute the friction velocity at the bottom
!   u_taub = rrb * sqrt(u(I,J,1,0)**2 + v(I,J,1,0)**2)

end do
u_taub = 0.d0
gotm_dt = tsp

!!! Keeping a parameter open for gotm - putting bottom roughness 0
 z0b_gotm = 0

!! ---- write(6,*) 'gradupdate'
!!! Update GOTM variables for tke, eps, l from PSOM
do k=0, NK
   tke(k)   = psom_tke(I,J,k) 
   eps1d(k) = psom_eps(I,J,k) 
   L(k) = psom_l(I,J,k) 

end do

u_taub = 0d0
z0b_gotm = 0d0

!!! Now, call do_turbulence subroutine, to update diffusivities
call do_turbulence(NK, gotm_dt, dpth, u_fric(I,J), u_taub, z0s_gotm, z0b_gotm, thk, NN1d, SS2, TKEP)

!!! Eliminate very low values/ Confine then within the limit 1.0d-05
do k=0, NK
!! Update diffusivities in PSOM
if (num(k).ge.(1.0d-05) .AND. num(k).lt.clip) then
      KzMom(I,J,k) = num(k)
    else if (num(k).ge.clip) then
      KzMom(I,J,k) = clip
    else
      KzMom(I,J,k) = 1.0d-05
endif

if (nuh(k).ge.(1.0d-005) .AND. nuh(k).le.clip) then  
     KzTr(I,J,k) = nuh(k)  
   else if (nuh(k).ge.clip) then  
     KzTr(I,J,k) = clip  
   else  
     KzTr(I,J,k)  = 1.0d-05  
endif  

!! Change tke, eps, l only when ivb=3/3rd RK step
if (ivb==3 ) then
!! Update tke, eps, L in PSOM
   psom_tke(I,J,k) = tke(k) 
   psom_eps(I,J,k) = eps1d(k) 
    psom_l(I,J,k) = L(k) 
end if

end do

return
end subroutine couple_gotm

