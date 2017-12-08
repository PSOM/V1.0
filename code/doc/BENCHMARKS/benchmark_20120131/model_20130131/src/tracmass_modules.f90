MODULE mod_param
 INTEGER                                   :: IMT, JMT, KM
! INTEGER, POINTER                          :: IMT, JMT, KM
 INTEGER                                   :: JMAX, LBT
 INTEGER, PARAMETER                        :: MR=1001
 INTEGER                                   :: NEND
 INTEGER, PARAMETER                        :: NST=2,NNRJ=8,NTRJ=7
 INTEGER, PARAMETER                        :: NTRACmax= 1 
 REAL(kind=rc_kind), ALLOCATABLE, DIMENSION(:,:)       :: trj
 INTEGER, ALLOCATABLE, DIMENSION(:,:)       :: nrj

 INTEGER                                   :: ncoor,kriva,iter,ngcm
 REAL(kind=rc_kind) ::  tseas,tday,tyear,dtmin,voltr
 REAL(kind=rc_kind) ::  tstep,dstep,tss,partQuant
 REAL(kind=rc_kind), PARAMETER                         :: UNDEF=1.d20 
 REAL(kind=rc_kind), PARAMETER                         :: PI =3.14159265358979323846d0
ENDMODULE mod_param

MODULE mod_grid
  REAL(kind=rc_kind), ALLOCATABLE, DIMENSION(:,:)       :: dxv, dyu
  REAL(kind=rc_kind), ALLOCATABLE, DIMENSION(:)         :: dz
  REAL(kind=rc_kind), ALLOCATABLE, DIMENSION(:,:)       :: dxdy
  REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: dvol
  REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: dzt
  
  REAL(kind=rc_kind) ::  rmin ,tmin ,smin
  REAL(kind=rc_kind) ::  dr ,dtemp ,dsalt
  REAL(kind=rc_kind) ::  arcscale
  INTEGER, ALLOCATABLE, DIMENSION(:,:)      :: kmt
  INTEGER                                   :: subGrid     ,subGridID
  INTEGER                                   :: subGridImin ,subGridImax
  INTEGER                                   :: subGridJmin ,subGridJmax
  CHARACTER(LEN=200)                        :: SubGridFile 
ENDMODULE mod_grid

MODULE mod_vel
  REAL(kind=rc_kind), ALLOCATABLE, DIMENSION(:,:,:,:)    :: uflux ,vflux
  REAL(kind=rc_kind), ALLOCATABLE, DIMENSION(:,:,:,:)    :: wflux
  REAL(kind=rc_kind), ALLOCATABLE, DIMENSION(:,:,:)      :: hs  
  REAL(kind=rc_kind) ::  ff
ENDMODULE mod_vel
