PROGRAM main 

! Definition of the parts of the code used
#include "cppdefs.h"

USE header                                 ! Declaration of variables

#ifdef relaxation
  use relaxation
#endif

#ifdef allow_particle
  USE particles
#endif

! Some additional declarations
#include "main_declarations.h"

!=
!= End of the declarations.
!=
!= ------------------------------------------------------------------------
!=
!= Beginning of the code by itself.
!=


! Namelist is read and major parameters (like qpr, kappah, lambda...) are defined:
#include "ini_param.h"


#ifdef relaxation
 call set_coef()
#endif

! call ini_grids()                          ! Initialization of the grid (xc,yc).


#ifdef allow_particle
 ALLOCATE(parti(NPR))
#endif

! Tr_p => Tr(1,:,:,:,0)                    ! define pointer, map the 5-dimensional Tr  to the 3-dimensional Tr_p for saving

 CALL ini_setup(pcorr) 

 CALL meanh(NI,NJ,h,hmean) 

 CALL sigma           ! calculates metric terms in vertical for the moving part of the grid
 CALL staticsigma     ! calculates metric terms in vertical for the fixed part of the grid
                      ! It is important to call sigma before staticsigma and ini_st          
                      ! sigma needs to be called in momentum after each time advancement  

 CALL tracerinit(0)   ! Initializes tracer

 CALL hsave           ! Saves initial free surface h, hsave must be called before geostroph
                             
 CALL ini_uv(0)       ! find geostrophically balanced velocities from T,s,h

 CALL facediv(EPS,fdiv)          ! checks the divergence using face fluxes
 CALL cdiv(EPS,ctrdiv,0)         ! checks the divergence using face fluxes

 CALL vort(0)                 ! calculates vorticity  -diagnostic

 step= 0 
 time_nondim= dtf*step
 time_seconds= time_nondim*TL

#ifdef gotm_call
   CALL initial_tke()
#endif

#ifdef implicit
   write(6,*) 'mixing implicit'        
#else
   write(6,*) 'mixing explicit'    
#endif

 CALL diag_n2 
 CALL diag_n2budget(step) 
 
 ! 20120406. At this time, dtime is not defined !!
 tim= dtime

 !     BCs set just once in the beginning                                
 CALL setbc(step) 
 CALL correctbc 

 CALL checks 

 ! Management of pickup options before the model starts.
 lv_test_output_bin=.FALSE.
#ifdef file_output
#ifdef file_output_bin
 lv_test_output_bin=.TRUE.
#endif
#endif

 IF(pickup_step<-0.5) then ! which means that the simulation will not start from a previous file,
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


 ! Output
#ifdef file_output

#ifdef file_output_cdf
 call write_cdf(step,0)
#endif

#ifdef file_output_bin
 call write_bin(step)     ! write_bin deals with both the binary output and the pickup options (read and write)
                         
#endif

#endif

 CALL diag_energy(step) 


!  The initialization if over.
!----------------------------------------------------------------------------------------




!----------------------------------------------------------------------------------------
! Main loop:

 DO step = initial_step+1,initial_step+nsteps
 !CALL tracerinit(step)    !initializes tracer

 time_nondim= dtf*step
 time_seconds= time_nondim*TL
#ifdef allow_particle
  ! initialize particles
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

   
  ! PRINT*, 'steps='//stepchar(step)

  CALL diag_n2

  CALL momentum(pcorr,step) 

  CALL diag_n2budget(step)                      !     calcn2 must be called before n2budget. Called at end of momentum  
  CALL diag_energy(step)                   !     compute the KE                                                    
  CALL meanh(NI,NJ,h,hmean) 

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

!#ifdef file_output
!#include "write_op.f90"
!#endif

#ifdef file_output

#ifdef file_output_cdf
 call write_cdf(step,0)
#endif

#ifdef file_output_bin
 call write_bin(step)
#endif

#endif



 ENDDO

!================================================================================================================




! close(91)
END PROGRAM main
