PROGRAM main 

   USE header                   ! Declaration of variables

#ifdef relaxation
   USE relaxation
#endif

#ifdef allow_particle
   USE particles
#endif

#include "cppdefs.h"
#include "main_declarations.h"  ! Some additional declarations
#include "ini_param.h"          ! Namelist is read and major parameters are defined

!==================================================================================================
!                                        Initialization 
!==================================================================================================

#ifdef relaxation
   call set_coef()
#endif

   ! call ini_grids()       ! Initialization of the grid (xc,yc).

#ifdef allow_particle
   ALLOCATE(parti(NPR))
#endif

   CALL ini_setup(pcorr) 
   CALL meanh(NI,NJ,h,hmean) 
   CALL sigma               ! calculates metric terms in vertical for the moving part of the grid
   CALL staticsigma         ! calculates metric terms in vertical for the fixed part of the grid
                            ! It is important to call sigma before staticsigma and ini_st          
                            ! sigma needs to be called in momentum after each time advancement  
!  Tr_p => Tr(1,:,:,:,0)    ! define pointer, map the 5D Tr to the 3D Tr_p for saving
!  CALL tracerinit(0)       ! Initializes tracer
   CALL hsave               ! Saves initial free surface h, hsave must be called before geostroph
   CALL ini_uv(0)           ! find geostrophically balanced velocities from T,s,h
   CALL facediv(EPS,fdiv)   ! checks the divergence using face fluxes
   CALL cdiv(EPS,ctrdiv,0)  ! checks the divergence using center fluxes
   CALL vort(0)             ! calculates vorticity - diagnostic

   step= 0 
   time_nondim= dtf*step
   time_seconds= time_nondim*TL
   time_days= time_seconds/86400.d0
   
#ifdef gotm_call
   CALL initial_tke()
#endif

! #ifdef implicit
!    write(6,*) 'mixing implicit'
! #else
!    write(6,*) 'mixing explicit'
! #endif

   CALL diag_n2 
   CALL diag_n2budget(step) 
   CALL checks              ! display in the terminal
 
   ! 20120406. At this time, dtime is not defined !!
   tim= dtime

   ! BCs set just once in the beginning - Not required to call these for the periodic case
   ! CALL setbc(step) 
   ! CALL correctbc 

! ---------------------------  Pickup  ---------------------------
   lv_test_output_bin=.FALSE.
#ifdef file_output
#ifdef file_output_bin
   lv_test_output_bin=.TRUE.
#endif
#endif

   IF(pickup_step<-0.5) then ! the simulation will not start from a previous file
      initial_step=0
      PRINT*,"  NO PICKUP"
   ELSE ! There will be a pickup
      IF(.NOT.(lv_test_output_bin)) then
         PRINT*,"Error: A pickup is required but the binary I/O is not activated"
         STOP
      ELSE
         PRINT*,"The simulation will start from step ",pickup_step
         initial_step=pickup_step
         step=initial_step
      ENDIF
   ENDIF

! ---------------------------  Output  --------------------------- 
#ifdef file_output
#ifdef file_output_cdf
   call write_cdf(step,0)
#endif
#ifdef file_output_bin
   call write_bin(step) ! write_bin deals with both the binary output and the pickup options
#endif
#endif

   CALL diag_energy(step) 

! ---------------------  Initial conditions  ---------------------
   do j=0,NJ+1
   do k=0,NK+1
      u_ini(j,k)= 0.d0
      do i=0,NI+1
         u_ini(j,k)= u_ini(j,k) + u(i,j,k,0)
      enddo
      u_ini(j,k)= u_ini(j,k)/(dble(NI)+2)
   enddo
   enddo
   
   if (restore_step == 0) then
      if (restore_sT) then
         T_restore= T_ini
         s_restore= s_ini
      endif
      if (restore_u) then
         u_restore= u_ini
      endif
      restore_rate=1./restore_time
   endif

! =================================================================================================
!                                          Main loop 
! =================================================================================================
   DO step = initial_step+1,initial_step+nsteps
!       CALL tracerinit(step)    !initializes tracer

      time_nondim= dtf*step
      time_seconds= time_nondim*TL
      time_days= time_seconds/86400.d0

#ifdef allow_particle
   if (step == ini_particle_time) then
      PRINT*,"INITIALIZATION OF PARTICLES ",NPR
      CALL ini_particles(step)   ! get the particle velocity for the first step
      CALL get_parti_vel(step)
      DO i = 1, NP
         parti(i)%i = parti(i)%i + dtf * parti(i)%u
         parti(i)%j = parti(i)%j + dtf * parti(i)%v
         parti(i)%k = parti(i)%k + dtf * parti(i)%w
         parti(i)%u0 = parti(i)%u
         parti(i)%v0 = parti(i)%v
         parti(i)%w0 = parti(i)%w
      ENDDO
   endif
#endif
   
      CALL diag_n2
      CALL momentum(pcorr,step) 
      CALL diag_n2budget(step)    ! Calcn2 must be called before n2budget. Called at end of momentum  
      CALL diag_energy(step)      ! Compute the KE. Display values every 10 steps             
      CALL meanh(NI,NJ,h,hmean) 

! ---------------------------  Particle  --------------------------- 
#ifdef allow_particle
   if (step>ini_particle_time-1) then
      CALL get_parti_vel(step)
      CALL parti_forward
      PRINT*,step,parti_outfreq,mod(step,parti_outfreq)
      if(mod(step,parti_outfreq)==0) then
         CALL save_parti
      endif
   endif
#endif

! ---------------------------  Restore  ---------------------------
!       if (step == initial_step+restore_step) then
!
!          if (restore_st) then
!             do j=0,NJ+1
!                do k=0,NK+1
!                   T_restore(j,k)= 0.d0
!                   s_restore(j,k)= 0.d0
!                   do i=0,NI+1
!                      T_restore(j,k)= T_restore(j,k) + T(i,j,k,0)
!                      s_restore(j,k)= s_restore(j,k) + s(i,j,k,0)
!                   enddo
!                   T_restore(j,k)= T_restore(j,k)/(dble(NI)+2)
!                   s_restore(j,k)= s_restore(j,k)/(dble(NI)+2)
!          !          rho_restore(j,k)= potdens(s_restore(j,k),T_restore(j,k))
!                enddo
!             enddo
!             restore_rate=1./(3600.*24.*restore_time)
!          end if
!
!          if (restore_u) then
!             do j=0,NJ+1
!                do k=0,NK+1
!                   u_restore(j,k)= 0.d0
!                   do i=0,NI+1
!                      u_restore(j,k)= u_restore(j,k) + u(i,j,k,0)
!                   enddo
!                   u_restore(j,k)= u_restore(j,k)/(dble(NI)+2)
!                enddo
!             enddo
!             restore_rate=1./(3600.*24.*restore_time)
!          end if
!
!       endif
! 
! ---------------------------  Output  ---------------------------
! #ifdef file_output
! #include "write_op.f90"
! #endif

#ifdef file_output

#ifdef file_output_cdf
   call write_cdf(step,0)
#endif

#ifdef file_output_bin
   call write_bin(step)
#endif

#endif

   ENDDO ! end main loop

!==================================================================================================

! close(91)
END PROGRAM main
