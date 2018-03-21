subroutine checks 

#include "cppdefs.h"
  !     ------------------                                                
  USE header
  integer i,j,k 
  integer n1, n2, n3
  real r1

 !     Checks and writes a few things                                    
                                                                        
                                                                   
 do j=1,NJ 
   do i=0,NI 
     if (abs(D(i+1,j)-D(i,j)).gt.1.d-4) then 
       write(6,*) i,j,D(i+1,j),D(i,j) 
       write(6,*) 'need to modify rpevalgrad (rgradients) for slope in x direcn: grpifc needs to be modified'     
       write(6,*) 'need to modify velbc_periodicew' 
       stop 
     end if 
   end do 
 end do 
                                                                   
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

 
 if(verbose) then
   write(6,*) "                                            "
   write(6,*) "############################################"
   write(6,*) "#                                           "
   write(6,*) "#----------------  MODEL SETTING     "
   write(6,*) "#                                           "


   write(6,*)             "#   time step   x number of time steps  = length of simulation "
   !write(6,"(A3,F8.2,A13,I6,A12)") " # ",dtime_dim," sec  x      ",nsteps,"           = "
  
   r1= nsteps*dtf*TL/3600.
   if(r1>24) then
     n1=INT(r1/24.);
     r1=r1-n1*24.; 
     write(6,"(A3,F8.2,A13,I6,A12,I5,A10,F7.3,A2)") " # ",dtime_dim," sec  x      ",nsteps,"           = ",n1," days and ",r1," h"
    else
     write(6,"(A3,F8.2,A13,I6,A12,F7.3,A2)") " # ",dtime_dim," sec  x      ",nsteps,"           = ",nsteps*dtf*TL/3600., "h"
   endif

   write(6,*) "#                                           "

#ifdef rhoonly
   write(6,*) "# only rho is used."
#else
   write(6,*) "# s and T are used."
#endif

   write(6,*) "#                                           "
 
   if(fplane==1) then
     write(6,"(A14,F5.2,A5)") " # f-plane at ",phi0deg, " deg "
    else
     write(6,"(A39,F5.2,A5)") " # f varies with latitude, centered on ",phi0deg," deg"
   endif

   if(fnhhy .ge. 0.9) then
     write(6,*) "# non-hydrostatic simulation "
    else
     write(6,*) "# hydrostatic approximation is used "
  
   endif

   IF(use_Shchepetkin) then
     write(6,*) "# Shchepetkin scheme ENABLED"
    ELSE
     if(lv_flat_bottom) then
       write(6,*) "# Shchepetkin scheme DISABLED"
       write(6,*) "# Song scheme ENABLED"
      else
       write(6,*) "# Shchepetkin scheme DISABLED"
       write(6,*) "------"
       write(6,"(A97)") "Error: You attempt to run a simulation with a non-flat bottom with the Song scheme.              "
       write(6,"(A97)") "         To use a baroclinic pressure term computation scheme that is accurate on sloping bottom," 
       write(6,"(A97)") "          use Shchepetkin scheme: in namelist, use_Shchepetkin=.TRUE.                            "
       STOP
     endif
   ENDIF

   write(6,"(A47,2(D11.3))") "# horizontal diffusion coefficients (m2.s-1) :",Kx,Ky

#ifdef fixed_bottom_thickness
   IF(bottom_linear_drag) then
     write(6,"(A46,D11.3)") " # bottom drag is linear, with coefficient RR=",RR
    ELSE
     write(6,*) "# bottom drag not linear: not supported."
     STOP
     !write(6,"(A46,D11.3)") "# bottom drag is quadratic, with coefficient RR=",RR
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
 
  write(6,*) "#                                            "
  write(6,*) "#----------------  NONDIMENSIONAL NUMBERS    "
  write(6,*) "#                                            "

  WRITE(6,"(A,E14.4)")  " # EPS     = ",EPS
  WRITE(6,"(A,E14.4)")  " # delta   = ",delta
  WRITE(6,"(A,E14.4)")  " # qpr     = ", qpr 
  WRITE(6,"(A,E14.4)")  " # lambda  = ",lambda
  WRITE(6,"(A,E14.4)")  " # beta    = ", beta
  WRITE(6,"(A,E14.4)")  " # dtf     = ", dtf 
  WRITE(6,"(A,E14.4)")  " # WL      = ", WL
  WRITE(6,"(A,E14.4)")  " # R0      = ", R0
  WRITE(6,"(A,E14.4)")  " # UL      = ", UL
  WRITE(6,"(A,E14.4)")  " # P1      = ", P1
  WRITE(6,"(A,E14.4)")  " # HL      = ", HL
  WRITE(6,"(A,E14.4)")  " # HDL     = ", HDL 
  

   write(6,*) "#                                           "
   write(6,*) "#----------------  GRID CHARACTERISTICS     "
   write(6,*) "#                                           "
   if(rect) then
     write(6,*) "# rectangular grid" 
    else 
     write(6,*) "# non-rectangular grid !"
   endif
   if(periodicew) then
     write(6,*) "# periodic in the EW direction" 
    else 
     write(6,*) "# non periodic in the EW direction (not supported!)"
   endif


   write(6,*) "# NI, NJ, NK: ",NI,NJ,NK
   write(6,"(A12,F8.2,A2)") " # delta_x: ",dx, " m"
   write(6,"(A12,F8.2,A2)") " # delta_y: ",dy, " m"


   write(6,*) "#                                           "

   write(6,"(A33,F8.4)") " # top cell thickness (m)      : ",dztop
   
#ifdef fixed_bottom_thickness
   write(6,"(A33,F8.4)") " # bottom cell thickness (m)   : ",dzbot
#else
   write(6,*) "# bottom cell height       : variable"
#endif

   write(6,"(A39,F7.3)") " # vertical grid stretching factor   : ",pfac
 
   if(lv_flat_bottom) then
     write(6,"(A27,F8.3,A2)") " # flat topography, depth: ",total_depth," m"
    else
     write(6,"(A40,F9.3,A7,F9.3,A2)") " # sloping topography, varying between: ",&
                                       &MINVAL(zf(:,:,0)*DL)," m and ",MAXVAL(zf(:,:,0)*DL)," m"
  endif

  write(6,*) "#                                           "
  write(6,*) "#-----------  PARTICLES CHARACTERISTICS     "
  write(6,*) "#                                           "

 
#ifdef allow_particle
   write(6,*) "# particles ENABLED."
   write(6,*) "# number of particles   : ",NPR
   write(6,*) "# initial time step     : ",ini_particle_time
   write(6,*) "# frequency of output   : ",parti_outfreq
   write(6,*) "# number of output files: ",parti_file_num
   write(6,"(A22,4(E14.7,' , '))") " # pcx, pcy, pcz, pcr: ",pcx,pcy,pcz,pcr
   
#else
   write(6,*) "# particles DISABLED. "
#endif 


   write(6,*) "#                                           "
   write(6,*) "#----------------  OUTPUT CHARACTERISTICS     "
   write(6,*) "#                                           "


#ifdef file_output

   write(6,"(A22,A101)") " # output ENABLED, in  ",dirout


#ifdef file_output_cdf
   write(6,*) "# output in cdf format ENABLED."
#else  
   write(6,*) "# output in cdf format DISABLED."
#endif

#ifdef file_output_bin
   write(6,*) "# output in bin format ENABLED."
#else  
   write(6,*) "# output in bin format DISABLED."
#endif


   write(6,*) "# output 1D frequency: ",out1d_int,"steps"
   write(6,*) "# output 2D frequency: ",out2d_int,"steps"
   write(6,*) "# output 3D frequency: ",out3d_int,"steps"


#else  

   write(6,*) "# output DISABLED."

#endif

   write(6,*) "#                                           "
   write(6,*) "#----------------  USER VARIABLES           "
   write(6,*) "#                                           "
   write(6,*) "# user1 : ",user1
   write(6,*) "# user2 : ",user2
   write(6,*) "# user3 : ",user3
   write(6,*) "# user4 : ",user4
   write(6,*) "# user5 : ",user5
   write(6,*) "# user6 : ",user6
   write(6,*) "# user7 : ",user7
   write(6,*) "# user8 : ",user8


   write(6,*) "#                                           "
   write(6,*) "############################################"
   write(6,*) "                                            "
   write(6,*) "                                            "



   write(6,*) "                                            "
   write(6,*) "                                            "
   write(6,*) "############################################"
   write(6,*) "#                                           "
   write(6,*) "#----------------  INITIAL SETTING     "
   write(6,*) "#                                           "
   write(6,"(A21,D15.7,A5,D15.7)") " #   u lies between: ",MINVAL(u(:,:,:,0))," and",MAXVAL(u(:,:,:,0))    
   write(6,"(A21,D15.7,A5,D15.7)") " #   v lies between: ",MINVAL(v(:,:,:,0))," and",MAXVAL(v(:,:,:,0))    
   write(6,"(A21,D15.7,A5,D15.7)") " #   w lies between: ",MINVAL(w(:,:,:,0))," and",MAXVAL(w(:,:,:,0))    

#ifdef rhoonly
   write(6,"(A21,D15.7,A5,D15.7,A11)") " # rho lies between: ",MINVAL(rho(:,:,:))," and",MAXVAL(rho(:,:,:))," (only rho)"    

#else
   write(6,"(A21,D15.7,A5,D15.7)") " #   s lies between: ",MINVAL(s(:,:,:,0))," and",MAXVAL(s(:,:,:,0))    
   write(6,"(A21,D15.7,A5,D15.7)") " #   T lies between: ",MINVAL(T(:,:,:,0))," and",MAXVAL(T(:,:,:,0))    
   write(6,"(A21,D15.7,A5,D15.7,A19)") " # rho lies between: ",MINVAL(rho(:,:,:))," and",MAXVAL(rho(:,:,:))," (rho from s and T)"
#endif
   write(6,"(A21,D15.7,A5,D15.7)") " #   h lies between: ",MINVAL(h(:,:))," and",MAXVAL(h(:,:))    
   write(6,*) "#                                           "
   write(6,*) "############################################"











 endif ! verobse






!                                                                       
!      open(unit=88,file='umax.out',status='new')                       
!      open(unit=89,file='umin.out',status='new')                       
!      close(88 89)                                                     
!      open(unit=79,file='vtimeser.out',status='new')                   
!      open(unit=78,file='utimeser.out',status='new')                   
!      close(78 79)                                                     
                                                                        
      return 
      END                                           
