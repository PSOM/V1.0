  
 !****************************
 !   READING NAMELIST FILE
 !****************************

      !--------------
      !- DEFINITION - 
      !--------------


  NAMELIST /PARAM/ nsteps,dtime_dim,fnhhy,fplane,phi0deg,dx,dy,dztop_dim,lv_flat_bottom,&
#ifdef fixed_bottom_thickness
                   dzbot_dim,&
#endif
                   total_depth,use_Shchepetkin,Kx,Ky,RR,&
                   pickup_int,pickup_step,out1d_int,out2d_int,out3d_int,dirout

#ifdef allow_particle
  NAMELIST /traj/ ini_particle_time,parti_file_num,parti_outfreq,NPR,pcx,pcy,pcz,pcr
#endif

  NAMELIST /user/ user1,user2,user3,user4,user5,user6,user7,user8



      !-------------
      !-  READING  - 
      !-------------


  READ (*,PARAM)

#ifdef allow_particle
  READ (*,traj)                        
#endif

  READ (*,user)                        


 
 !******************************************
 !   INITIALIZATION OF NONDIM PARAMETERS
 !******************************************


     !--------------------------
     ! Numerical matters:

  kappah= 0.65d0 

  kaphinv= 1.d0/kappah 


     !--------------------------
     ! Geometric constants:

  lambda= LEN/AL 

  delta= DL/LEN 

  delinv= LEN/DL

     !--------------------------
     ! Vertical grid:

  dztop=dztop_dim/DL
#ifdef fixed_bottom_thickness
  dzbot=dzbot_dim/DL
#endif

      !--------------------------
      ! Characteristics of the flow: 

  qpr= fnhhy*delta 

  beta= 1.d0/(EPS*EPS*delta) 

  UL= FPAR *LEN*EPS 

  WL= EPS*delta*UL 

  P1= R0*UL*UL*EPS**(-1) !     P1= RU^2/EPS                                                      

  HL= P1/(R0*10.d0)      ! 10.d0 stands for the normalization of g. 

  HDL= HL/DL 

  TL= LEN/UL 

     !--------------------------
     ! Time step:
 
 dtf=dtime_dim/TL
 
