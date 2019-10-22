subroutine advection_and_mixing(m,n,dtimel,step) 
#include "cppdefs.h"
!     ---------------------------------------------                     
USE header

INTEGER :: i,j,k
INTEGER :: n,m,step

REAL(kind=rc_kind) :: dtimel

REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1) :: var,var0              !var=(s,T,u,v,w,Tr(it,:,:,:)                    
REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1) :: uvarx, uvarx_advec    !uvarx  is the source term, divergence of the 
                                                                                  !advective fluxes 
REAL(kind=rc_kind), dimension(    1:NI  ,1:NJ  , 1:NK  ) :: vardif                !vardif is the source term from diabatic processes 
REAL(kind=rc_kind), dimension(    1:NI  ,1:NJ  , 1:NK  ) :: vardif2               !vardif is the source term from diabatic processes 
REAL(kind=rc_kind), dimension(    1:NI  ,1:NJ  , 1:NK  ) :: T_h,u_w, v_w       

INTEGER :: iv_compute_kz
INTEGER :: av_comp(5+ntr)

iv_compute_kz=1

!Hcontent = 0.d0
!Hmixing  = 0.d0

! By default, all variables will be subject to the following loop
av_comp=1

! If the simulation is a "rhoonly" simulation (rho is stored in s), T=0 and does not need to be updated.
#ifdef rhoonly
  av_comp(1)=0
#endif

u_w=0.;v_w=0.
T_h=0.;

call heat_flux  (T_h,step)
call wind_stress(u_w,v_w,step)

do selvar=1,5+ntr

#ifdef implicit
  mat_A = 0.d0; mat_B = 0.d0; mat_C = 0.d0; mat_D = 0.d0; mat_test = 0.d0;
#endif

  ! Flux terms initialized to 0
  uvarx = 0.d0; uvarx_advec = 0.d0; vardif = 0.d0

  if(av_comp(selvar)==1) then

    if(selvar==1) then; var=T(:,:,:,m); var0(:,:,:) = T(:,:,:,0); endif
    if(selvar==2) then; var=s(:,:,:,m); var0(:,:,:) = s(:,:,:,0); endif
    if(selvar==3) then; var=u(:,:,:,m); var0(:,:,:) = u(:,:,:,0); endif
    if(selvar==4) then; var=v(:,:,:,m); var0(:,:,:) = v(:,:,:,0); endif
    if(selvar==5) then; var=w(:,:,:,m); var0(:,:,:) = w(:,:,:,0); endif
    do it=1,ntr
      if(selvar==5+it) var=Tr(it,:,:,:,m)
    enddo

    ! ---------------------------------------------------------------
    ! computation of the advective fluxes, using QUICK scheme
    CALL advect(var,uvarx)
                         
    ! ---------------------------------------------------------------
    ! computation of the horizontal diabatic fluxes
    if(selvar<5.5) then

      vardif=0.; vardif2=0.
      call mixing_horizontal(var,vardif);
     !call mixing_isopycnal(var,vardif,10.);PRINT*,"VISCOUS REDI";
     !call mixing_isopycnal_biharmonic(var,vardif,1.);PRINT*,"BIHARMONIC REDI"

      uvarx(1:NI,1:NJ,1:NK)=uvarx(1:NI,1:NJ,1:NK)-vardif(1:NI,1:NJ,1:NK)
    endif

    uvarx_advec(1:NI,1:NJ,1:NK) = uvarx(1:NI,1:NJ,1:NK)
                       
    ! ---------------------------------------------------------------
    ! computation of the vertical diabatic fluxes
    if(selvar<4.5) then
      vardif=0.;
      call mixing_vertical(var,vardif,m,step,iv_compute_kz); 

      uvarx(1:NI,1:NJ,1:NK)=uvarx(1:NI,1:NJ,1:NK)-vardif(1:NI,1:NJ,1:NK)
      iv_compute_kz=0;
    endif

    ! ---------------------------------------------------------------
    !* addition of the source term for heat and momentum flux at top surface
    if(selvar==1) then

#ifdef implicit
      mat_D(1:NI,1:NJ,1:NK,selvar)=T_h(1:NI,1:NJ,1:NK)
#else
      uvarx(1:NI,1:NJ,1:NK)=uvarx(1:NI,1:NJ,1:NK)-T_h(1:NI,1:NJ,1:NK)
#endif

    endif

    if(selvar==3) then

#ifdef implicit
      mat_D(1:NI,1:NJ,1:NK,selvar)=u_w(1:NI,1:NJ,1:NK)
#else
      uvarx(1:NI,1:NJ,1:NK)=uvarx(1:NI,1:NJ,1:NK)-u_w(1:NI,1:NJ,1:NK)
#endif

    endif
  
    if(selvar==4) then

#ifdef implicit
      mat_D(1:NI,1:NJ,1:NK,selvar)=v_w(1:NI,1:NJ,1:NK)
#else
      uvarx(1:NI,1:NJ,1:NK)=uvarx(1:NI,1:NJ,1:NK)-v_w(1:NI,1:NJ,1:NK)
#endif

    endif

    ! ---------------------------------------------------------------
#ifdef implicit

#else
    if(selvar==3 .OR. selvar==4) then
      uvarx(1:NI,1:NJ,1)=uvarx(1:NI,1:NJ,1) + RR*(1.d0/(UL*delta) ) * ( Jac(1:NI,1:NJ,1)*wz(1:NI,1:NJ,1) ) * var(1:NI,1:NJ,1)
    end if
#endif

    ! ---------------------------------------------------------------
    ! final summation  

#ifdef implicit
      if (selvar.lt.4.5) then
        do i=1,NI
        do j=1,NJ

          mat_A(i,j,1:NK,selvar)   = dtimel*Jacinv(i,j,1:NK)*mat_A(i,j,1:NK,selvar)
          mat_B(i,j,1:NK,selvar)   = ( (dtimel*Jacinv(i,j,1:NK)*mat_B(i,j,1:NK,selvar)) - 1)
          mat_C(i,j,1:NK-1,selvar) = dtimel*Jacinv(i,j,1:NK-1)*mat_C(i,j,1:NK-1,selvar)
          mat_D(i,j,1:NK,selvar)   = 0.d0 - var0(i,j,1:NK ) - (dtimel*Jacinv(i,j,1:NK)*mat_D(i,j,1:NK,selvar)) + (dtimel*Jacinv(i,j,1:NK)*uvarx_advec(    i,j,1:NK))

          call solve_tridiag(mat_A(i,j,1:NK,selvar), mat_B(i,j,1:NK,selvar), mat_C(i,j,1:NK,selvar), mat_D(i,j,1:NK,selvar), mat_test(i,j,1:NK,selvar), NK )

        end do
        end do
      end if

      if(selvar==1) T(1:NI,1:NJ,1:NK,n)  = mat_test(1:NI,1:NJ,1:NK,selvar)
      if(selvar==2) s(1:NI,1:NJ,1:NK,n)  = mat_test(1:NI,1:NJ,1:NK,selvar)
      if(selvar==3) cx(1:NI,1:NJ,1:NK)   = mat_test(1:NI,1:NJ,1:NK,selvar)
      if(selvar==4) cy(1:NI,1:NJ,1:NK)   = mat_test(1:NI,1:NJ,1:NK,selvar)
      !if(selvar==5) cz(1:NI,1:NJ,1:NK)   = mat_test(1:NI,1:NJ,1:NK,selvar) 
#else
      if(selvar==1) T(1:NI,1:NJ,1:NK,n)=T(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
      if(selvar==2) s(1:NI,1:NJ,1:NK,n)=s(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
      if(selvar==3) cx(1:NI,1:NJ,1:NK) =u(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
      if(selvar==4) cy(1:NI,1:NJ,1:NK) =v(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
      !if(selvar==5) cz(1:NI,1:NJ,1:NK) =w(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
#endif

    if(selvar==5) cz(1:NI,1:NJ,1:NK) =w(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)

    do it=1,ntr
      if(selvar==5+it) Tr(it,1:NI,1:NJ,1:NK,n) =Tr(it,1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
    enddo

  endif
enddo !selvar

#ifdef allow_particle
   PRINT*,"ALLOW PARTICLE"
#endif


return 
                                                                        
END

                                           
