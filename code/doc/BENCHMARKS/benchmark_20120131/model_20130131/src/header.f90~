MODULE header

 ! Definition of the kind of reals:

 integer, parameter :: rc_kind = selected_real_kind(12) 

 ! Definition of the parts of the code used:
#include "cppdefs.h"




!------------------------------------------------------------------------!
!---                             CONSTANTS                           ----!
!------------------------------------------------------------------------!


 ! Definition of (NI, NJ, NK, ngrid, maxout, maxint and int1) and (ntr and nconsume)
 ! NI,NJ,NK are the no. of points in each of the 3 co-ordintate directions.

#include "size.h"  ! Grid size and figures for mgrid.
 !-------------------------------------------------------------------
 !-- physical properties
 !----------------------
  
      REAL(kind=rc_kind),  PARAMETER :: S0=35.7d0, T0=15.d0 ,R0=1027.d0    !@ S0: mean salinity, T0: mean temperature
                                                                           !@   : sreal= S0+s,Treal=T0+T
                                                                           !@ R0: mean density

 !-------------------------------------------------------------------
 !-- physical constants
 !---------------------

  REAL(kind=rc_kind),  PARAMETER :: PI=3.14159265358979323846 

  REAL(kind=rc_kind),  PARAMETER :: OMEGA=7.272d-5, gpr= 0.981,apr=0.6371  !@ OMEGA: Earth angular velocity,
                                                                           !@ gpr: gravitational acceleration, apr: Earth radius

  REAL(kind=rc_kind),  PARAMETER :: EPS= 0.1d0, AL=1.d7, FPAR=1.d-4        !@ EPS: Rossby number, AL: magnitude of the Earth radius,
                                                                           !@ FPAR: magnitude of Coriolis parameter.



 !-------------------------------------------------------------------
 !-- settings of the model
 !------------------------

  LOGICAL, PARAMETER :: rect = .TRUE., periodicew = .TRUE.

  REAL(kind=rc_kind),  PARAMETER :: dztop=1.d-3,z0= 0.2d0,zm= 0.1d0  !@ dztop: depth of the top layer non-dim by DL

  REAL(kind=rc_kind), PARAMETER :: LEN= 1.d5, DL=1.d3           !@ LEN: characteristic length scale, DL: characteristic depth scale

  REAL(kind=rc_kind) :: DLinv                                   !@ DLinv: inverse of DL.




!------------------------------------------------------------------------!
!---                             NON-CONSTANTS                       ----!
!------------------------------------------------------------------------!


 !-------------------------------------------------------------------
 !-- settings of the code
 !------------------------

  CHARACTER(LEN=101) :: dirout                                       !@ dirout: directory where the output will go, in namelist



 !-------------------------------------------------------------------
 !-- settings of the model
 !------------------------

  REAL(kind=rc_kind) :: fnhhy                                        !@ fnhhy: 0 for hydrostatic, 1 for nonhydrostatic, in namelist

  REAL(kind=rc_kind) :: total_depth                                  !@ total_depth: in namelist

  REAL(kind=rc_kind) ::  dtf                                         !@ dtf: nondimensional time step, in namelist

  REAL(kind=rc_kind) :: dx,dy                                        !@ dx, dy: dimensional, in namelist

  REAL(kind=rc_kind) :: dtime                                        !@ dtime: time step, in namelist

  REAL(kind=rc_kind) :: phi0deg                                      !@ phi0deg: central latitude, in namelist

  REAL(kind=rc_kind) ::  Kx_TS, Ky_TS_bdy,Ky_TS_int, Kx_m, Ky_m            !@ variables used in diffusion, in namelist

  REAL(kind=rc_kind) :: kappah,kaphinv                               !@ kappah: implcitness parameter in the Crank-Nicolson scheme f
                                                                     !@ kaphinv: inverse of kappah     


 !-------------------------------------------------------------------
 !-- derived constants
 !---------------------


  REAL(kind=rc_kind) :: lambda                                  !@ lambda: ratio of horizontal length scale to earth's radius = L/A 

  REAL(kind=rc_kind) ::  delta,delinv                           !@ delta: ratio of horizontal to vertical length scale = L/D 
                                                                !@ delinv: inverse of delta

  REAL(kind=rc_kind) :: qpr                                     !@ qpr: ratio of non-hyd to hyd pressure = Q/P, qpr= fnhhy*delta

  REAL(kind=rc_kind) :: beta                                    !@ beta: beta= 1.d0/(EPS*EPS*delta) 

  REAL(kind=rc_kind) :: P1                                      !@ P1: P1= R0*UL*UL*EPS**(-1) 

  REAL(kind=rc_kind) :: UL                                      !@ UL: UL= FPAR *LEN*EPS

  REAL(kind=rc_kind) :: WL                                      !@ WL: WL= EPS*delta*UL

  REAL(kind=rc_kind) :: TL                                      !@ TL: TL= LEN/UL

  REAL(kind=rc_kind) :: HL                                      !@ HL: HL= P1/(R0*10.d0)

  REAL(kind=rc_kind) :: HDL                                     !@ HDL: HDL= HL/DL



 !-------------------------------------------------------------------
 !-- more or less local variables
 !------------------------


  REAL(kind=rc_kind)  :: pfac                                   !@ pfac: grid stretching factor in z, used in findzall and sigma

  REAL(kind=rc_kind)  :: stressmax                              !@ stressmax: defines the maximum stress, used in diffusion_wind.

  REAL(kind=rc_kind)  :: yfront,dyfront,sclwidth,tightness,mldepth   !@ variables that define the initial state, used in stprofile.

  REAL(kind=rc_kind)  :: sigrelease(ntr)                        !@ sigrelease: isopycnal of tracer release, used in tracer*

  REAL(kind=rc_kind)  :: drho                                   !@ drho: used in init_tr and inith

  INTEGER :: conv(0:NI+1,0:NJ+1,0:NK+1)                         !@ conv: used in conadjust.

  integer :: NP,NPR,ini_particle_time,parti_file_num            !@ variables used for particles

  INTEGER :: con100(0:NI+1,0:NJ+1,0:NK+1)                       !@ con100: used in conadjust

  INTEGER :: ksurf                                              !@ ksurf: used in netcdf output routines.

  REAL(kind=rc_kind), dimension(3) :: mldn2,mldn2init,zbtop,zbbot,frictop,&
     diatop,zbtopinit,zbbotinit,advecpv,friction, diabatic,diabot,fricbot   !@ used in n2budget_topbot and writen2budget

  REAL(kind=rc_kind) :: sbackgrnd                               !@ sbackgrnd: obsolete. Was used in setbc and surfaceflux.
  REAL(kind=rc_kind) ::  r_sponge(0:NJ+1)                       !@ r_sponge: used in relaxation.



  ! added 20120302
 INTEGER :: it
 REAL(kind=rc_kind) :: hmean
 REAL(kind=rc_kind) ::  dum  



 ! added 20120302
 ! netcdf output

  integer :: iddatfile, idigit, idudx, idudy, idudz, idvbysection, idvby, idvbz, idvcon, idvcy, idvbx
  integer :: idvcz, idvc, idvdivfeddy, idvdivfreyn, idvdx, idvdz, idvd, idvfb, idvh, idvn2bar
  integer :: idvn2, idvnsq100m, idvnsq30m, idvpe,idvpsiv,idvpsiw,idvpv, idvp,idvrhbar,idvrho,idvrnk
  integer :: idvstrain,idvstress,idvstr,idvs,idvtbar,idvtemp,idvtim,idvtr,idvt,idvu,idvvb,idvvc,idvvor,idvv
  integer :: idvwb,idvwc,idvwpv,idvw,idvy,idvzsave,idvz,idwdx,idwdy,idwdz,iimday,ipos


  


!------------------------------------------------------------------------!
!---                              ARRAYS                             ----!
!------------------------------------------------------------------------!

 
  REAL(kind=rc_kind), target, dimension(ntr,0:NI+1,0:NJ+1, 0:NK+1,0:1)  :: Tr
  REAL(kind=rc_kind), pointer, dimension(:,:,:)  :: Tr_p

  REAL(kind=rc_kind), dimension(ntr,0:NI+1,0:NJ+1, 0:NK+1    )  :: wtr
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1,0:1)  :: u,v,w,s,T
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK,  3  )  :: gqk
  REAL(kind=rc_kind), dimension(    0:NI,    NJ,     NK,  2  )  :: gi,gqi
  REAL(kind=rc_kind), dimension(      NI,  0:NJ,     NK,  2  )  :: gj,gqj
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1,-1:NK+1)      :: zf
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1)      :: zc,wx,wy,wz,p,strain,shear,Jac,rho,vor,pv,freqN2,&
                                                                   &si,sj,sk,cx,cy,cz,rp,T_ref
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK  )      :: wt,wzk,skfc
  REAL(kind=rc_kind), dimension(      NI,    NJ,   0:NK  )      :: czf,Kz,wf
  REAL(kind=rc_kind), dimension(           0:NJ+1,   NK  )      :: ueast,uwest
  REAL(kind=rc_kind), dimension(      NI,  0:NJ,     NK  )      :: hyn,sjfc,cyf,gj3,gqj3,Jjfc,vf
  REAL(kind=rc_kind), dimension(    0:NI,    NJ  ,   NK  )      :: hxn,sifc,cxf,gi3,gqi3,Jifc,uf
  REAL(kind=rc_kind), dimension(      NI,    NJ,     NK  )      :: uvis,vvis,wvis,fricu,fricv,fricw,fricb,rhoadv,rhoprev
  REAL(kind=rc_kind), dimension(             NJ,     NK  )      :: ufbce,ufbcw,trinit,divreyn,divmean,dcdt,prod
  REAL(kind=rc_kind), dimension(      NI,            NK  )      :: vfbcn,vfbcs,vnorth,vsouth,ssouth
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1,   2   )      :: gradhn
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1        )      :: ux,uy,vx,vy,ffc,bbc,oldh,h,hdt,D,J2d,Ddx,Ddy,g11,g22,g12
  REAL(kind=rc_kind), dimension(      NI  ,0:NJ          )      :: bbj,ffj
  REAL(kind=rc_kind), dimension(    0:NI  ,  NJ          )      :: ffi,bbi
  REAL(kind=rc_kind), dimension(      NI,    NJ          )      :: wfbcb
  REAL(kind=rc_kind), dimension(           0:NJ+1        )      :: yc,latrad
  REAL(kind=rc_kind), dimension(    0:NI+1               )      :: xc
  REAL(kind=rc_kind), dimension(             NJ          )      :: stressx

  ! Array for output
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1,nconsume) :: consump 


  ! Arrays containing rp variations
  REAL(kind=rc_kind) :: drpx(NI,NJ,NK),drpy(NI,NJ,NK),grpifc(0:NI,NJ,NK),grpjfc(NI,0:NJ,NK)



END MODULE header

