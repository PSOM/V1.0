subroutine checks 
!     Checks and writes a few things 
#include "cppdefs.h"
  !     ------------------                                                
  USE header
  integer i,j,k 
  integer n1, n2, n3
  real r1
                                                              
! do j=1,NJ 
!   do i=0,NI 
!     if (abs(D(i+1,j)-D(i,j)).gt.1.d-4) then 
!       write(6,*) i,j,D(i+1,j),D(i,j) 
!       write(6,*) 'need to modify rpevalgrad (rgradients) for slope in x direcn: grpifc needs to be modified'     
!       write(6,*) 'need to modify velbc_periodicew' 
!!       stop 
!     end if 
!   end do 
! end do 
                                                                   
 do j=0,NJ+1 
   do i=0,NI+1 
     if ((uy(i,j).ne.0.d0).or.(vx(i,j).ne.0.d0)) then 
       write(6,*) 'Need to modify biharmonic for curvi grid' 
       stop 
     end if 
   end do 
 end do 
                                                                   
 if (rect.eqv.(.false.)) then  !  need to modify grpifc,grpjfc to contain cross diff terms       
   write(6,*) 'modify grpifc,grpjfc, stop in checks' 
   stop 
 end if 

 ! =======================================================================================
 if(verbose) then
   write(6,*) "                                                                          "
   write(6,*) "##########################################################################"
   write(6,*) "                                                                          "
   write(6,*) "----------------------------  MODEL SETTING  -----------------------------"
   write(6,*) "                                                                          "
   write(6,*) "# time step x number of time steps = length of simulation "
  
   r1= nsteps*dtf*TL/3600.
   if(r1>24) then
     n1=INT(r1/24.);
     r1=r1-n1*24.; 
     write(6,"(A3,F5.0,A13,I6,A11,I3,A10,F5.2,A2)") " # ",dtime_dim," sec x       ",nsteps,"         = ",n1," days and ",r1," h"
   else
     write(6,"(A3,F5.0,A13,I6,A12,F7.3,A2)") " # ",dtime_dim," sec x       ",nsteps,"          = ",nsteps*dtf*TL/3600., " h"
   endif
   
   ! -------------------------------------------------------------------------------------
!    write(6,*) "                                                  "
   
   if (restore_sT) then
     write(6,*) "# Restoration of salinity and temperature  ENABLED" 
   else
     write(6,*) "# Restoration of salinity and temperature  DISABLED" 
   end if
   if (restore_u) then
     write(6,*) "# Restoration of the zonal mean velocities ENABLED" 
   else
     write(6,*) "# Restoration of the zonal mean velocities DISABLED"
   end if
   if (restore_sT .OR. restore_u) then
     write(6,"(A,F5.0,A)") " # Restoring time scale:", restore_time," sec"
     write(6,"(A,I5,A)")   " # Restoration starts at", restore_step," step"
   end if

   ! =====================================================================================
   write(6,*) "                                                                          "
   write(6,*) "-------------------------  GRID CHARACTERISTICS  -------------------------"
   write(6,*) "                                                                          "
   if(rect) then
     write(6,*) "# Rectangular grid" 
    else 
     write(6,*) "# Non-rectangular grid"
   endif
#ifdef sigma_stretch
   write(6,"(A,D10.3)")     " # Vertical grid stretching factor     : ",pfac
#else
   write(6,*)               "# No vertical stretching."
#endif
   write(6,"(A,3(I8))")     " # Number of grids NI, NJ, NK          : ",NI,NJ,NK
   write(6,"(A,2(F8.2),A2)")" # Grid size delta_x, delta_y          : ",dx,dy, " m"
#ifdef fixed_bottom_thickness
   write(6,"(A,2(F8.2),A2)")" # Top & bottom layer thickness        : ",dztop*DL,dzbot*DL, " m" 
#else
   write(6,"(A,F8.2,A2)")   " # Top & bottom layer thickness        : ",dztop*DL, "variable m" 
#endif
   
   ! ------------------------------------------------------------------------------------- 
!    write(6,*) "                                                     "
   if(lv_flat_bottom) then
     write(6,"(A,F8.2,A2)") " # Flat topography, depth              : ",depmean_dim," m"
   else
     write(6,"(A,F8.2,A2)") " # Topographic wavelength              : ",topo_wavelength, " m"
     write(6,"(A,F8.2,A2)") " # Topographic amplitude               : ",topo_amplitude, " m"
     write(6,"(A,F8.2,A2)") " # Averaged depth of the topography    : ",depmean_dim, " m"     
     write(6,"(A,F8.2,A7,F8.2,A2)") " # Sloping topography, varying between : ",&
                                      & MINVAL(zf(:,:,0)*DL)," m and ",MAXVAL(zf(:,:,0)*DL), " m"
     write(6,"(A,F8.2,A2)") " # Total thickness of corrugated layers: ",distance, " m"
     write(6,"(A,I8)")      " # Number of corrugated layers         : ",Nf 
   endif
!    write(6,"(A,D10.3)")     " # Initial density gradient            : ",ini_drho
!    write(6,"(A,D10.3)")     " # Initial density tightness           : ",ini_tight
   
   ! =====================================================================================
   write(6,*) "                                                                          "
   write(6,*) "----------------------  SIMULATION CHARACTERISTICS  ----------------------"
   write(6,*) "                                                                          "
   
   if(periodicew) then
     write(6,*) "# Periodic in the EW direction" 
   else 
     write(6,*) "# Non periodic in the EW direction (not supported!)"
   endif

#ifdef rhoonly
   write(6,*) "# Only rho is used"
#else
   write(6,*) "# s and T are used"
#endif
 
   if(fplane==1) then
     write(6,"(A,F5.2,A4)") " # f-plane at ",phi0deg," deg"
   else
     write(6,"(A,F5.2,A4)") " # f varies with latitude, centered on ",phi0deg," deg"
   endif

   if(fnhhy .ge. 0.9) then
     write(6,*) "# Non-hydrostatic simulation"
   else
     write(6,*) "# Hydrostatic approximation is used"
   endif

   IF(use_Shchepetkin) then
     write(6,*) "# Shchepetkin scheme ENABLED"
   ELSE
!      if(lv_flat_bottom) then
!        write(6,*) "# Shchepetkin scheme DISABLED"
!        write(6,*) "# Song scheme ENABLED"
!      else
       write(6,*) "# Shchepetkin scheme DISABLED"
!        write(6,*) "------"
!        write(6,"(A97)") "Error: You attempt to run a simulation with a non-flat bottom with the Song scheme.              "
!        write(6,"(A97)") "         To use a baroclinic pressure term computation scheme that is accurate on sloping bottom,"
!        write(6,"(A97)") "          use Shchepetkin scheme: in namelist, use_Shchepetkin=.TRUE.                            "
!        STOP
!      endif
   ENDIF
   
   ! -------------------------------------------------------------------------------------
!    write(6,*) "                                                  "
   
#ifdef biharmonic_horizontal
   write(6,*) "# Horizontal biharmonic mixing ENABLED"
   write(6,"(A,  D10.3)")  " # Horizontal mixing coefficient (m4.s-1): ",binu
#else
   write(6,*) "# Horizontal biharmonic mixing DISABLED"
   write(6,"(A,2(D10.3))")  " # Horizontal mixing coefficient (m2.s-1): ",Kx,Ky
#endif
   write(6,"(A,1(D10.3))") " # Vertical mixing coefficient (m2.s-1)  : ",KzMom

#ifdef fixed_bottom_thickness
   IF(bottom_linear_drag) then
     write(6,"(A,D10.3)")  " # Bottom drag is linear with coefficient: ",RR
   ELSE
     write(6,*) " # Bottom drag not linear (not supported!)"
     STOP
     !write(6,"(A47,D10.3)") "# Bottom drag is quadratic, with coefficient RR =",RR
   ENDIF
#else
   IF(ABS(RR)>1e-13 .AND. .NOT.(lv_flat_bottom)) then
     write(6,*) "------"
     write(6,*) "Error: Attempt to use bottom friction on a sloping topography without keeping the height of the lowermost cell."
     write(6,*) "       You have to modify at least one of these parameters:"
     write(6,*) "        - fixed_bottom_thickness in inc/cppdefs.h (currently undef)"
     write(6,*) "        - RR in namelist (currently non zero)"
     write(6,*) "        - lv_flat_bottom in namelist (currently .FALSE.)"
     stop
   ENDIF
#endif   
    
   ! =====================================================================================
!    write(6,*) "                                                                          "
!    write(6,*) "-----------------------  PARTICLES CHARACTERISTICS  ----------------------"
!    write(6,*) "                                                                        "

#ifdef allow_particle
   write(6,*) "# Particles ENABLED."
   write(6,*) "# Number of particles   : ",NPR
   write(6,*) "# Initial time step     : ",ini_particle_time
   write(6,*) "# Frequency of output   : ",parti_outfreq
   write(6,*) "# Number of output files: ",parti_file_num
   write(6,"(A22,4(E14.7,' , '))") " # pcx, pcy, pcz, pcr: ",pcx,pcy,pcz,pcr
#else
   write(6,*) "# Particles DISABLED. "
#endif 

   ! =====================================================================================
   write(6,*) "                                                                          "
   write(6,*) "------------------------  OUTPUT CHARACTERISTICS  ------------------------"
   write(6,*) "                                                                          "

#ifdef file_output
   write(6,"(A22,A101)") " # Output ENABLED, in  ",dirout

#ifdef file_output_cdf
   write(6,*) "# Output in cdf format ENABLED"
#else  
   write(6,*) "# Output in cdf format DISABLED"
#endif

#ifdef file_output_bin
   write(6,*) "# Output in bin format ENABLED"
#else  
   write(6,*) "# Output in bin format DISABLED"
#endif

   write(6,*) "# Output 1D frequency: ",out1d_int,"steps"
   write(6,*) "# Output 2D frequency: ",out2d_int,"steps"
   write(6,*) "# Output 3D frequency: ",out3d_int,"steps"

#else  
   write(6,*) "# output DISABLED."
#endif

   ! =====================================================================================
   write(6,*) "                                                                          "
   write(6,*) "---------------------------  INITIAL SETTING  ----------------------------"
   write(6,*) "                                                                          "
   write(6,"(A21,D15.7,A5,D15.7)") " #   u lies between: ",MINVAL(u(:,:,:,0))*UL," and",MAXVAL(u(:,:,:,0))*UL    
   write(6,"(A21,D15.7,A5,D15.7)") " #   v lies between: ",MINVAL(v(:,:,:,0))*UL," and",MAXVAL(v(:,:,:,0))*UL    
   write(6,"(A21,D15.7,A5,D15.7)") " #   w lies between: ",MINVAL(w(:,:,:,0))*WL," and",MAXVAL(w(:,:,:,0))*WL    

#ifdef rhoonly
   write(6,"(A21,D15.7,A5,D15.7,A11)") " # rho lies between: ",MINVAL(rho(:,:,:))," and",MAXVAL(rho(:,:,:))," (only rho)"    
#else
   write(6,"(A21,D15.7,A5,D15.7)") " #   s lies between: ",MINVAL(s(:,:,:,0))," and",MAXVAL(s(:,:,:,0))    
   write(6,"(A21,D15.7,A5,D15.7)") " #   T lies between: ",MINVAL(T(:,:,:,0))," and",MAXVAL(T(:,:,:,0))    
   write(6,"(A21,D15.7,A5,D15.7,A19)") " # rho lies between: ",MINVAL(rho(:,:,:))," and",MAXVAL(rho(:,:,:))," (rho from s and T)"
#endif
   write(6,"(A21,D15.7,A5,D15.7)") " #   h lies between: ",MINVAL(h(:,:))," and",MAXVAL(h(:,:))    

   write(6,*) "                                                                          "
   write(6,*) "-----------------------  NONDIMENSIONAL NUMBERS  -------------------------"
   write(6,*) "                                                                          "

   WRITE(6,"(A,D15.7)")  " # EPS     = ",EPS
   WRITE(6,"(A,D15.7)")  " # delta   = ",delta
   WRITE(6,"(A,D15.7)")  " # qpr     = ",qpr 
   WRITE(6,"(A,D15.7)")  " # lambda  = ",lambda
   WRITE(6,"(A,D15.7)")  " # beta    = ",beta
   WRITE(6,"(A,D15.7)")  " # R0      = ",R0
   WRITE(6,"(A,D15.7)")  " # P1      = ",P1
   WRITE(6,"(A,D15.7)")  " # dtf     = ",dtf 
   WRITE(6,"(A,D15.7)")  " # LEN     = ",LEN
   WRITE(6,"(A,D15.7)")  " # DL      = ",DL
   WRITE(6,"(A,D15.7)")  " # UL      = ",UL
   WRITE(6,"(A,D15.7)")  " # WL      = ",WL
   WRITE(6,"(A,D15.7)")  " # HL      = ",HL
   WRITE(6,"(A,D15.7)")  " # HDL     = ",HDL 
   WRITE(6,"(A,D15.7)")  " # TL      = ",TL

   write(6,*) "                                                                          "
   write(6,*) "##########################################################################"
   write(6,*) "                                                                          "
   
 endif ! verobse

!      open(unit=88,file='umax.out',status='new')                       
!      open(unit=89,file='umin.out',status='new')                       
!      close(88 89)                                                     
!      open(unit=79,file='vtimeser.out',status='new')                   
!      open(unit=78,file='utimeser.out',status='new')                   
!      close(78 79)                                                     

      return
      END