PROGRAM main 

#include "cppdefs.h"                       ! Definition of the parts of the code used

USE header                                 ! Declaration of variables

#ifdef relaxation
  use relaxation      ! Use of relaxation
#endif

use grids                                  ! Declarations dealing with the grid. 

#ifdef particle
  USE particles       ! Use of particles
#endif

#include "main_declarations.f90"           ! Some additional declarations


!=
!= End of the declarations.
!=
!= ------------------------------------------------------------------------
!=
!= Beginning of the code by itself.
!=


#include "ini_param.f90"                   ! The namelist is read and major "variables" (like qpr, kappah, lambda...) are defined.

#ifdef relaxation
 call set_coef()
#endif

 call ini_grids()                          ! Initialization of the grid (xc,yc).


#ifdef particle
 ALLOCATE(parti(NPR))
#endif

 Tr_p => Tr(1,:,:,:,0)                    ! define pointer, map the 5-dimensional Tr  to the 3-dimensional Tr_p for saving

 CALL init(pcorr) 

 CALL meanh(NI,NJ,h,hmean) 
 WRITE(*,*) 'initial  hmean',hmean


 CALL sigma           ! calculates metric terms in vertical for the moving part of the grid
 CALL staticsigma     ! calculates metric terms in vertical for the fixed part of the grid
                      ! It is important to call sigma before staticsigma and stprofile          
                      ! sigma needs to be called in momentum after each time advancement  

 CALL tracerinit()    !initializes tracer

 CALL hsave            ! saves initial free surface h, hsave must be called before geostroph
                             
 CALL geostroph(0)     ! find geostrophically balanced velocities from T,s,h

 CALL facediv(EPS,fdiv)          ! checks the divergence using face fluxes
 CALL cdiv(EPS,ctrdiv,0)         ! checks the divergence using face fluxes

 CALL vort(0)                 ! calculates vorticity  -diagnostic

 step= 0 

 CALL calcn2 
 CALL n2budget(step) 
 
 ! 20120406. At this time, dtime is not defined !!
 PRINT*, "# dtime =",dtime
 tim= dtime


 !     BCs set just once in the beginning                                
 CALL setbc(step) 
 CALL correctbc 

 CALL checks 

#ifdef file_output
#ifdef file_output_cdf
   ksurf=NK;call writeksurf(frame_int,step,ksurf,h,consump,Tr,s,T,rho,u,v,w,p,vor,strain,freqN2,xc,yc,zc,DL,LEN,Jac,dtf*TL)
   ksurf= INT(NK/3); call writeksurf(frame_int,step,ksurf,h,consump,Tr,s,T,rho,u,v,w,p,vor,strain,freqN2,xc,yc,zc,DL,LEN,Jac,dtf*TL)
#endif
#endif

 CALL energy(step) 


!  The initialization if over.
!----------------------------------------------------------------------------------------
! Main loop:

 DO step = pickup_step+1,pickup_step+nsteps

#ifdef particle
  ! initialize particles
  if (step == ini_particle_time) then
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

   
  PRINT*, 'steps='//stepchar(step)

  CALL momentum(pcorr,step) 

  CALL n2budget(step)                      !     calcn2 must be called before n2budget. Called at end of momentum  
  CALL energy(step)                        !     compute the KE                                                    
  CALL meanh(NI,NJ,h,hmean) 



#ifdef particle
  if (step>ini_particle_time) then
   !particles
   CALL get_parti_vel(step)
   CALL parti_forward
   call save_parti
  endif
#endif

#ifdef file_output
#include "write_op.f90"
#endif



 ENDDO

!================================================================================================================




! close(91)
END PROGRAM main
