REAL(kind=rc_kind), dimension( NK ) :: Tnorth,Tsouth,snorth,ssouth
!  REAL, PARAMETER :: yrday_start=1.
!  REAL :: mld(0:NI+1,0:NJ+1)
!  REAL :: qflux

  INTEGER,  PARAMETER :: lenarray= 1140    ! length of input array - 2 months of fluxes
  REAL(kind=rc_kind) :: yrday(lenarray),swrday(lenarray),qlossday(lenarray),daycount(lenarray)
  REAL(kind=rc_kind) :: uwind(lenarray),vwind(lenarray),precip(lenarray),stressx(lenarray),stressy(lenarray)

