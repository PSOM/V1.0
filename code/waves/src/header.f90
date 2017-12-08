MODULE header

 ! Definition of the kind of reals:
 integer, parameter :: rc_kind = selected_real_kind(12) 

 ! Definition of the parts of the code used:
#include "cppdefs.h"


 !-------------------------------------------------------------------
 !-- grid size
 !----------------------

 
                                                               ! Grid size and array sizes for mgrid.
#include "size.h"
                                                               ! Definition of (NI, NJ, NK, ngrid, maxout, maxint and int1) and (ntr and nconsume)
                                                               ! NI,NJ,NK are the no. of points in each of the 3 co-ordintate directions.



 !-------------------------------------------------------------------
 !-- physical properties
 !----------------------
  
  REAL(kind=rc_kind),  PARAMETER :: S0=35.7d0, T0=15.d0 ,R0=1027.d0  !@ S0: mean salinity, T0: mean temperature
                                                                     !@   : sreal= S0+s,Treal=T0+T
                                                                     !@ R0: mean density


 !-------------------------------------------------------------------
 !-- physical constants
 !---------------------

  REAL(kind=rc_kind),  PARAMETER :: PI=3.14159265358979323846 

  REAL(kind=rc_kind),  PARAMETER :: OMEGA=7.272d-5, FPAR=1.d-4   !@ OMEGA: Earth angular velocity,
                                                                 !@ FPAR: magnitude of Coriolis parameter.
  REAL(kind=rc_kind),  PARAMETER :: gpr= 0.981                   !@ gpr: gravitational acceleration
  REAL(kind=rc_kind),  PARAMETER :: apr=0.6371                   !@ apr: Earth radius

  REAL(kind=rc_kind),  PARAMETER :: EPS= 0.1d0, AL=1.d7          !@ EPS: Rossby number, AL: magnitude of the Earth radius,



 !-------------------------------------------------------------------
 !-- settings of time stepping
 !----------------------------

  INTEGER :: nsteps                                           !@ number of time steps, in namelist

  REAL(kind=rc_kind) ::  dtf                                  !@ dtf: nondimensional time step
  REAL(kind=rc_kind) ::  dtime_dim                            !@ dtf:    dimensional time step (in sec), in namelist
  REAL(kind=rc_kind) ::  dtime                                !@ dtime: time step (obsolete)
                                                            
  REAL(kind=rc_kind) ::  time_nondim,time_seconds             !@ incremental time (clock)


  REAL(kind=rc_kind) :: kappah,kaphinv                        !@ kappah: implcitness parameter in the Crank-Nicolson scheme f
                                                              !@ kaphinv: inverse of kappah     


 !-------------------------------------------------------------------
 !-- settings of the grid
 !------------------------


  REAL(kind=rc_kind), PARAMETER :: LEN= 1.d5, DL=1.d3         !@ LEN: characteristic length scale, DL: characteristic depth scale

  LOGICAL, PARAMETER :: rect = .TRUE., periodicew = .TRUE.
  REAL(kind=rc_kind) :: dx,dy                                 !@ dx, dy: dimensional (in m), in namelist
  REAL(kind=rc_kind) :: dyM(0:NJ+1)


  LOGICAL :: lv_flat_bottom                                !@ choice for topography, in namelist
  REAL(kind=rc_kind) :: total_depth                        !@ depth of the domain (in m) (only if lv_flat_bottom), in namelist
  LOGICAL :: use_Shchepetkin                               !@ switch for the baroclinic pressure computation, in namelist


  REAL(kind=rc_kind),  PARAMETER :: z0= 0.2d0,zm= 0.1d0 
  REAL(kind=rc_kind) :: dztop                              !@ dztop    : depth of the top layer non-dim by DL
  REAL(kind=rc_kind) :: dztop_dim                          !@ dztop_dim: depth of the top layer (in m), in namelist

#ifdef fixed_bottom_thickness
  REAL(kind=rc_kind) :: dzbot                              !@ dzbot    : if fixed_bottom_thickness, thickness of the lowermost cell
  REAL(kind=rc_kind) :: dzbot_dim                          !@ dzbot_dim: dimensional thickness of the top layer (in m), in namelist
#endif



 !-------------------------------------------------------------------
 !-- settings of the model
 !------------------------


  !- Switch for hydrostatic
 
  REAL(kind=rc_kind) :: fnhhy                                 !@ fnhhy: 0 for hydrostatic, 1 for nonhydrostatic, in namelist
                                                            
 
  !- Planetary vorticity                                      

  REAL(kind=rc_kind) :: phi0deg                               !@ phi0deg: central latitude, in namelist
  INTEGER :: fplane                                           !@ fplane: 0 for latitude-dependance of f,
                                                              !          1 for a constant f            , in namelist



 !-------------------------------------------------------------------
 !-- diffusion and friction
 !--------------------------


  REAL(kind=rc_kind) ::  Kx, Ky                              !@ variables used in mixing_horizontal (in m*m/s), in namelist
                                                           
  LOGICAL, PARAMETER :: bottom_linear_drag=.TRUE.            !@ only linear drag is currently supported (see in mixing_vertical routine). 
  REAL(kind=rc_kind) ::  RR                                  !@ bottom friction (in m/s), initialized in namelist



 !-------------------------------------------------------------------
 !-- derived constants
 !---------------------


  REAL(kind=rc_kind) :: lambda                               !@ lambda: ratio of horizontal length scale to earth's radius = L/A 
                                                            
  REAL(kind=rc_kind) :: DLinv                                !@ DLinv: inverse of DL.
                                                            
  REAL(kind=rc_kind) ::  delta,delinv                        !@ delta: ratio of horizontal to vertical length scale = L/D 
                                                             !@ delinv: inverse of delta
                                                            
  REAL(kind=rc_kind) :: qpr                                  !@ qpr: ratio of non-hyd to hyd pressure = Q/P, qpr= fnhhy*delta
                                                            
  REAL(kind=rc_kind) :: beta                                 !@ beta: beta= 1.d0/(EPS*EPS*delta) 
                                                            
  REAL(kind=rc_kind) :: P1                                   !@ P1: P1= R0*UL*UL*EPS**(-1) 
                                                            
  REAL(kind=rc_kind) :: UL                                   !@ UL: UL= FPAR *LEN*EPS
                                                            
  REAL(kind=rc_kind) :: WL                                   !@ WL: WL= EPS*delta*UL
                                                            
  REAL(kind=rc_kind) :: TL                                   !@ TL: TL= LEN/UL
                                                            
  REAL(kind=rc_kind) :: HL                                   !@ HL: HL= P1/(R0*10.d0)
                                                            
  REAL(kind=rc_kind) :: HDL                                  !@ HDL: HDL= HL/DL


 !-------------------------------------------------------------------
 !-- user-defined constants
 !---------------------

   REAL(KIND=rc_kind) :: user1, user2, user3, user4,&       !@ 8 variables that users can modify to make various simulations. 
                       & user5, user6, user7, user8         !@  initialized in namelist.                                                                        
                                                             
                                                             
                                                                         
                                                                         
 !-------------------------------------------------------------------
 !-- settings of the output                                              
 !--------------------------                                             
                                                                         
  CHARACTER(LEN=151) :: dirout                              !@ dirout: directory where the output will go, in namelist
                                                           
  INTEGER :: out1d_int, out2d_int, out3d_int                !@ outXd_int: frequency of output in X dimensions ( X=1, 2, 3 ),
                                                            !@   in namelist  

  !INTEGER :: save_steps                                   
                                                           
                                                           
  LOGICAL, PARAMETER :: verbose=.TRUE.                      !@ if verbose, the model will produce more messages on the screen.
                                                                         
                                                                         
                                                                         
 !-------------------------------------------------------------------
 !-- more or less local variables                                        
 !------------------------                                               
                                                                         

#ifdef allow_particle
  integer :: NP,NPR,ini_particle_time,parti_file_num               !@ variables used for particles
  integer :: parti_outfreq                                        
  REAL(KIND=rc_kind) :: pcx, pcy, pcz, pcr                        
#endif                                                            
                                                                  
  integer ivb,ivs,ivf                                              !@ For use in momentum.f90
                                                                  
  REAL(kind=rc_kind)  :: pfac                                      !@ pfac: grid stretching factor in z, used in findzall and sigma
                                                                  
!  REAL(kind=rc_kind)  :: stressmax                                 !@ stressmax: defines the maximum stress, used in diffusion_wind.
  REAL(kind=rc_kind), dimension(NI,NJ) :: stress_top,stress_top_x,stress_top_y

  REAL(kind=rc_kind)  :: yfront,dyfront
  REAL(kind=rc_kind)  :: sclwidth,tightness,mldepth                !@ variables that define the initial state, used in ini_st.

  REAL(kind=rc_kind)  :: sigrelease(ntr)                           !@ sigrelease: isopycnal of tracer release, used in tracer*
                                                                  
  REAL(kind=rc_kind)  :: drho                                      !@ drho: used in init_tr and inith
                                                                  
  INTEGER :: conv(0:NI+1,0:NJ+1,0:NK+1)                            !@ conv: used in conadjust.
  INTEGER :: con100(0:NI+1,0:NJ+1,0:NK+1)                          !@ con100: used in conadjust

  REAL(kind=rc_kind), dimension(3) :: mldn2,mldn2init,&
                                     &zbtop,zbbot,&
                                     &zbbotinit,zbtopinit,&
                                     &fricbot,frictop,&
                                     &diabot,diatop,&
                                     &advecpv,friction, diabatic   !@ used in n2budget_topbot and writen2budget

  REAL(kind=rc_kind) :: sbackgrnd                                  !@ sbackgrnd: obsolete. Was used in setbc and surfaceflux.
  REAL(kind=rc_kind) :: r_sponge(0:NJ+1), r_T(0:NJ+1,0:NK+1)      ! r_sponge: used in relaxation.
  REAL(kind=rc_kind) :: rho_refS(0:NK+1), rho_refN(0:NK+1), rho_refB(0:NJ+1)
                                                                
 REAL(kind=rc_kind) :: hmean                                       !@ hmean: mean level of the free surface


 INTEGER :: it
 REAL(kind=rc_kind) ::  dum  


                                                                        
 !-------------------------------------------------------------------
 !-- netcdf output                                        
 !------------------------                                               
  
  integer :: iddatfile, idigit, idudx, idudy, idudz, idvbysection, idvby, idvbz, idvcon, idvcy, idvbx
  integer :: idvcz, idvc, idvdivfeddy, idvdivfreyn, idvdx, idvdz, idvd, idvfb, idvh, idvn2bar
  integer :: idvn2, idvnsq100m, idvnsq30m, idvpe,idvpsiv,idvpsiw,idvpv, idvp,idvrhbar,idvrho,idvrnk
  integer :: idvstrain,idvstress,idvstr,idvs,idvtbar,idvtemp,idvtim,idvtr,idvt,idvu,idvvb,idvvc,idvvor,idvv
  integer :: idvwb,idvwc,idvwpv,idvw,idvy,idvzsave,idvz,idwdx,idwdy,idwdz,iimday,ipos
  INTEGER :: frame_int,ngraph2d

                                                                      
 !-------------------------------------------------------------------
 !-- pickup option                                        
 !------------------------                                               

  INTEGER :: pickup_step,pickup_int
  LOGICAL :: lv_test_output_bin

  INTEGER :: initial_step
  


!------------------------------------------------------------------------!
!---                              ARRAYS                             ----!
!------------------------------------------------------------------------!

 
  REAL(kind=rc_kind), dimension(ntr,0:NI+1,0:NJ+1, 0:NK+1,0:1)  :: Tr
  !REAL(kind=rc_kind), target, dimension(ntr,0:NI+1,0:NJ+1, 0:NK+1,0:1)  :: Tr
  !REAL(kind=rc_kind), pointer, dimension(:,:,:)  :: Tr_p

  REAL(kind=rc_kind), dimension(ntr,0:NI+1,0:NJ+1, 0:NK+1    )  :: wtr
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1,0:1)  :: u,v,w,s,T
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1,0:1)  :: u_bar,v_bar,w_bar,s_bar,T_bar
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1,0:1)  :: u_pri,v_pri,w_pri,s_pri,T_pri
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK,  3  )  :: gqk
  REAL(kind=rc_kind), dimension(    0:NI,    NJ,     NK,  2  )  :: gi,gqi
  REAL(kind=rc_kind), dimension(      NI,  0:NJ,     NK,  2  )  :: gj,gqj
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1,-1:NK+1)      :: zf
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1)      :: zc,wx,wy,wz,p,strain,shear,Jac,Jacinv,rho,rho_bar,rho_pri,vor,pv
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1)      :: freqby,freqbz,freqN2,si,sj,sk,cx,cy,cz,rp,T_ref
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK  )      :: wt,wzk,skfc
  REAL(kind=rc_kind), dimension(      NI,    NJ,   0:NK  )      :: czf,Kz,wf
  REAL(kind=rc_kind), dimension(           0:NJ+1,   NK  )      :: ueast,uwest
  REAL(kind=rc_kind), dimension(      NI,  0:NJ,     NK  )      :: hyn,sjfc,cyf,gj3,gqj3,Jjfc,vf
  REAL(kind=rc_kind), dimension(    0:NI,    NJ  ,   NK  )      :: hxn,sifc,cxf,gi3,gqi3,Jifc,uf
  REAL(kind=rc_kind), dimension(      NI,    NJ,     NK  )      :: uvis,vvis,wvis,fricu,fricv,fricw,fricb,rhoadv,rhoprev
  REAL(kind=rc_kind), dimension(             NJ,     NK  )      :: ufbce,ufbcw,trinit,divreyn,divmean,dcdt,prod
  REAL(kind=rc_kind), dimension(      NI,            NK  )      :: vfbcn,vfbcs,vnorth,vsouth,ssouth
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1,   2   )      :: gradhn
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1        )      :: ux,uy,vx,vy,ffc,bbc,oldh,h,h_bar,h_pri,hdt,&
                                                                    &D,J2d,Ddx,Ddy,g11,g22,g12
  REAL(kind=rc_kind), dimension(      NI  ,0:NJ          )      :: bbj,ffj
  REAL(kind=rc_kind), dimension(    0:NI  ,  NJ          )      :: ffi,bbi
  REAL(kind=rc_kind), dimension(      NI,    NJ          )      :: wfbcb
  REAL(kind=rc_kind), dimension(           0:NJ+1        )      :: yc,latrad
  REAL(kind=rc_kind), dimension(    0:NI+1               )      :: xc
  !REAL(kind=rc_kind), dimension(             NJ          )      :: stressx
 
! --- Sonaljit: I am defining wind stress component stressy here, instead of defining it locally in mixing_vertical routine
  !REAL(kind=rc_kind) :: stressy

  ! Array for output
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1,nconsume) :: consump 


  ! Arrays containing rp variations
  REAL(kind=rc_kind) :: drpx(NI,NJ,NK),drpy(NI,NJ,NK),grpifc(0:NI,NJ,NK),grpjfc(NI,0:NJ,NK)
  REAL(kind=rc_kind) :: ru_Sche(0:NI+1,0:NJ+1,NK),rv_Sche(0:NI+1,0:NJ+1,NK)
  REAL(kind=rc_kind) :: ru2_Sche(0:NI+1,0:NJ+1,NK),rv2_Sche(0:NI+1,0:NJ+1,NK)

  REAL(kind=rc_kind) :: ru3_Sche(0:NI+1,0:NJ+1,NK),rv3_Sche(0:NI+1,0:NJ+1,NK)
  REAL(kind=rc_kind) :: ru4_Sche(0:NI+1,0:NJ+1,NK),rv4_Sche(0:NI+1,0:NJ+1,NK)
  REAL(kind=rc_kind) :: diffS(0:NI+1,0:NJ+1,NK),diffS2(0:NI+1,0:NJ+1,NK)
  REAL(kind=rc_kind) :: diffS3(0:NI+1,0:NJ+1,NK),diffS4(0:NI+1,0:NJ+1,NK)

  ! Arrays for pv and vorticity diagnosis
  REAL(kind=rc_kind) :: pvt(NI,NJ,NK),pv1(NI,NJ,NK),pv2(NI,NJ,NK),pv3(NI,NJ,NK),vor1(NI,NJ,NK),vor2(NI,NJ,NK),vor3(NI,NJ,NK)

!------------------------------------------------------------------------!
!---                         GOTM Variables                          ----!
!------------------------------------------------------------------------!

  ! Assigning a integer for sign changing
  INTEGER :: sign, kppcall, gotmlogic, kzlogic, selvar

                                                                           
  ! Assigning Von Karman constant, reference density, cm0, Mellor Yamada parameter B1
  REAL(kind=rc_kind) :: k_von = 0.41, rhoref = 1027, cmu0_psom = 5.27046261454271d-1, B1_MY = 16.6
 
  ! ASSIGNING THE TKE, EPS, MACROLENGTH, AT CELL CENTERS
  REAL(kind=rc_kind), dimension(  0:NI+1,0:NJ+1,0:NK+1  ) :: psom_tke, psom_eps, psom_l
  
  ! ASSIGNING A RESULTANT WIND STRESS AND FRICTION VELOCITY PROFILE. THESE ARE 2-D ARRAYS
  REAL(kind=rc_kind), dimension(  0:NI+1,0:NJ+1  )    :: stx, u_fric
  REAL(kind=rc_kind), dimension(  0:NI+1,0:NJ+1,0:NK+1  ) :: lthk            

  ! MOMENTUM VISCOSITY AND HEAT DIFFUSIVITY
  REAL(kind=rc_kind), dimension(  0:NI+1,0:NJ+1,0:NK  ) :: KzMom, KzTr     
  REAL(kind=rc_kind), dimension(  0:NI+1,0:NJ+1,0:NK  ) :: nnface, ssface, Rig

  ! I NEED THE VALUES OF TURBULENCE PARAMETERS AT THE LAST STEP
  ! TEST
  REAL(kind=rc_kind) :: tsp
  !  REAL(kind=rc_kind) :: tsp1 = 216.d0, tsp2 = 108.d0, tsp3 = 72.d0

  ! Shear, ushear, vshear
  REAL(kind=rc_kind), dimension(  0:NI+1,0:NJ+1,0:NK  ) :: ush, vsh, shq

  ! Just for getting layer thickness values
  REAL(kind=rc_kind), dimension(  1:NK  ) :: thickn

  ! For buoyancy flux
  REAL(kind=rc_kind) :: qflux

  ! For calculating the number of times 'do_turbulence' subroutine is called
  INTEGER :: cnt

  ! dudz and dvdz at cell face
  REAL(kind=rc_kind), dimension(  0:NI+1,0:NJ+1,0:NK  ) :: udzf, vdzf, rdzf

  ! Implicit mixing
  REAL(kind=rc_kind), dimension(    1:NI,1:NJ,1:NK,1:5   )      :: mat_A, mat_B, mat_C, mat_D, mat_test

  ! Jerlov type variables
  REAL(kind=rc_kind) :: J_lambda1, J_lambda2, J_A

  ! Assigning variables for buoyancy/heat flux. Short wave, heatloss, and Jerlov water type parameters
  REAL(kind=rc_kind), dimension(0:NJ+1)  ::  swr, qloss

  REAL(kind=rc_kind), dimension(  0:NK  )  :: thk, SS2, TKEP, NN1d

  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1)      :: rho_old
! ------------------------------------------------------------------------------------



!------------------------------------------------------------------------!
!---                  EXPERIMENT-SPECIFIC VARIABLES                  ----!
!------------------------------------------------------------------------!

! This include directive enables users to define more variables in the header,
!  that are specific to their experiment without modifying the header source file.
!
! By default, header_expe.h is empty.

#include "header_expe.h"





!------------------------------------------------------------------------!
!---                       END OF MODULE HEADER                      ----!
!------------------------------------------------------------------------!
END MODULE header

