MODULE mod_nc
  USE netcdf
  USE header
!  use mod_integrate
  integer :: ncid_u,varid_u,rec,ncid_v,varid_v,ncid_p,varid_p,ncid_q,varid_q
CONTAINS
  SUBROUTINE create_ncfiles()
    CALL create_nc_u(ncid_u,varid_u)
    PRINT*, "# create nc file for u field with ID:",ncid_u
!    CALL create_nc_v(ncid_v,varid_v)
!    PRINT*, "# create nc file for v field with ID:",ncid_v
!    CALL create_nc_p(ncid_p,varid_p)
!    PRINT*, "# create nc file for p field with ID:",ncid_p
!    CALL create_nc_q(ncid_q,varid_q)
!    PRINT*, "# create nc file for q field with ID:",ncid_q
  END SUBROUTINE create_ncfiles

  SUBROUTINE dump_data()
    CALL dump_u(ncid_u,varid_u,rec)
!    CALL dump_v(ncid_v,varid_v,rec)
!    CALL dump_p(ncid_p,varid_p,rec)
!    CALL dump_q(ncid_q,varid_q,rec)
  END SUBROUTINE dump_data

  SUBROUTINE close_ncfiles()
    CALL close_nc(ncid_u)
!    CALL close_nc(ncid_v)
!    CALL close_nc(ncid_p)
!    CALL close_nc(ncid_q)
  END SUBROUTINE close_ncfiles

  SUBROUTINE create_nc_q(ncid,varid)
    IMPLICIT NONE
    REAL(kind=rc_kind) ::  lats(NJ), lons(NI)
    INTEGER :: ncid,varid
    CHARACTER (len = *), PARAMETER :: FILE_NAME="q.nc",VAR_NAME="q"
    INTEGER, PARAMETER :: NDIMS = 3
    INTEGER :: start(NDIMS), countt(NDIMS),rec
    CHARACTER (len = *), PARAMETER :: LAT_NAME = "Y"
    CHARACTER (len = *), PARAMETER :: LON_NAME = "X"
    CHARACTER (len = *), PARAMETER :: REC_NAME = "time"
    INTEGER :: lon_dimid, lat_dimid, rec_dimid

    INTEGER :: lon_varid, lat_varid
    INTEGER :: dimids(NDIMS)
    CHARACTER (len = *), PARAMETER :: UNITS = "1/s"
    CHARACTER (len = *), PARAMETER :: LAT_UNITS = "meter"
    CHARACTER (len = *), PARAMETER :: LON_UNITS = "meter"

    lats=y_q(:,1)
    lons=x_q(1,:)

    ! Create the file. 
    CALL check( nf90_create(casename//'_'//FILE_NAME, nf90_clobber, ncid) )
    CALL check( nf90_def_dim(ncid, LAT_NAME, Ny, lat_dimid) )
    CALL check( nf90_def_dim(ncid, LON_NAME, Nx, lon_dimid) )
    CALL check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )
    CALL check( nf90_def_var(ncid, LAT_NAME, NF90_REAL, lat_dimid, lat_varid) )
    CALL check( nf90_def_var(ncid, LON_NAME, NF90_REAL, lon_dimid, lon_varid) )
    ! Assign units attributes to coordinate variables.
    CALL check( nf90_put_att(ncid, lat_varid, "units","meter") )
    CALL check( nf90_put_att(ncid, lon_varid, "units", "meter") )
    dimids = (/ lon_dimid, lat_dimid, rec_dimid /)
    ! Define the netCDF variables 
    CALL check( nf90_def_var(ncid, VAR_NAME, NF90_FLOAT, dimids, varid) )
    ! Assign units attributes to the netCDF variables.
    CALL check( nf90_put_att(ncid, varid, "units",UNITS) )
    CALL assign_common_att(ncid,varid)

    CALL check( nf90_enddef(ncid) )
    CALL check( nf90_put_var(ncid, lat_varid, lats) )
    CALL check( nf90_put_var(ncid, lon_varid, lons) )
  END SUBROUTINE create_nc_q

  SUBROUTINE dump_q(ncid,varid,rec)
    INTEGER ::ncid,varid,rec,start(3),countt(3)
    countt = (/ Nx, Ny, 1 /)
    start = (/ 1, 1, 1/)
    start(3) = rec
    CALL check( nf90_put_var(ncid, varid,transpose(qm), start = start, &
         count = countt) )
  END SUBROUTINE dump_q

  SUBROUTINE create_nc_p(ncid,varid)
    IMPLICIT NONE
    REAL(kind=rc_kind) ::  lats(Ny+1), lons(Nx+1)
    INTEGER :: ncid,varid
    CHARACTER (len = *), PARAMETER :: FILE_NAME="p.nc",VAR_NAME="p"
    INTEGER, PARAMETER :: NDIMS = 3
    INTEGER :: start(NDIMS), countt(NDIMS),rec
    CHARACTER (len = *), PARAMETER :: LAT_NAME = "Y"
    CHARACTER (len = *), PARAMETER :: LON_NAME = "X"
    CHARACTER (len = *), PARAMETER :: REC_NAME = "time"
    INTEGER :: lon_dimid, lat_dimid, rec_dimid
    INTEGER :: lon_varid, lat_varid
    INTEGER :: dimids(NDIMS)
    CHARACTER (len = *), PARAMETER :: UNITS = "m^2/s"
    CHARACTER (len = *), PARAMETER :: LAT_UNITS = "meter"
    CHARACTER (len = *), PARAMETER :: LON_UNITS = "meter"

    lats=y_p(:,1)
    lons=x_p(1,:)

    ! Create the file. 
    CALL check( nf90_create(casename//'_'//FILE_NAME, nf90_clobber, ncid) )
    CALL check( nf90_def_dim(ncid, LAT_NAME, Ny+1, lat_dimid) )
    CALL check( nf90_def_dim(ncid, LON_NAME, Nx+1, lon_dimid) )
    CALL check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )
    CALL check( nf90_def_var(ncid, LAT_NAME, NF90_REAL, lat_dimid, lat_varid) )
    CALL check( nf90_def_var(ncid, LON_NAME, NF90_REAL, lon_dimid, lon_varid) )
    ! Assign units attributes to coordinate variables.
    CALL check( nf90_put_att(ncid, lat_varid, "units","meter") )
    CALL check( nf90_put_att(ncid, lon_varid, "units", "meter") )

    ! The dimids array is used to pass the dimids of the dimensions of
    ! the netCDF variables. Both of the netCDF variables we are creating
    ! share the same four dimensions. In Fortran, the unlimited
    ! dimension must come last on the list of dimids.
    dimids = (/ lon_dimid, lat_dimid, rec_dimid /)

    ! Define the netCDF variables 
    CALL check( nf90_def_var(ncid, VAR_NAME, NF90_REAL, dimids, varid) )

    ! Assign units attributes to the netCDF variables.
    CALL check( nf90_put_att(ncid, varid, "units",UNITS) )
    CALL assign_common_att(ncid,varid)

    CALL check( nf90_enddef(ncid) )
    CALL check( nf90_put_var(ncid, lat_varid, lats) )
    CALL check( nf90_put_var(ncid, lon_varid, lons) )
  END SUBROUTINE create_nc_p
  SUBROUTINE dump_p(ncid,varid,rec)
    INTEGER ::ncid,varid,rec,start(3),countt(3)
    countt = (/ Nx+1, Ny+1, 1 /)
    start = (/ 1, 1, 1 /)
    start(3) = rec
    CALL check( nf90_put_var(ncid, varid, TRANSPOSE(pm), start = start, &
         count = countt) )
  END SUBROUTINE dump_p


  SUBROUTINE create_nc_v(ncid,varid)
    IMPLICIT NONE
    REAL(kind=rc_kind) ::  lats(NJ+2), lons(NI+2)
    INTEGER :: ncid,varid
    CHARACTER (len = *), PARAMETER :: FILE_NAME="v.nc",VAR_NAME="v"
    INTEGER, PARAMETER :: NDIMS = 3
    INTEGER :: start(NDIMS), countt(NDIMS),rec
    CHARACTER (len = *), PARAMETER :: LAT_NAME = "Y"
    CHARACTER (len = *), PARAMETER :: LON_NAME = "X"
    CHARACTER (len = *), PARAMETER :: REC_NAME = "time"
    INTEGER :: lon_dimid, lat_dimid, rec_dimid

    INTEGER :: lon_varid, lat_varid
    INTEGER :: dimids(NDIMS)
    CHARACTER (len = *), PARAMETER :: UNITS = "m/s"
    CHARACTER (len = *), PARAMETER :: LAT_UNITS = "meter"
    CHARACTER (len = *), PARAMETER :: LON_UNITS = "meter"

    lats=y_v(:,1)
    lons=x_v(1,:)

    ! Create the file. 
    CALL check( nf90_create(casename//'_'//FILE_NAME, nf90_clobber, ncid) )
    CALL check( nf90_def_dim(ncid, LAT_NAME, Ny+1, lat_dimid) )
    CALL check( nf90_def_dim(ncid, LON_NAME, Nx, lon_dimid) )
    CALL check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )
    CALL check( nf90_def_var(ncid, LAT_NAME, NF90_REAL, lat_dimid, lat_varid) )
    CALL check( nf90_def_var(ncid, LON_NAME, NF90_REAL, lon_dimid, lon_varid) )
    ! Assign units attributes to coordinate variables.
    CALL check( nf90_put_att(ncid, lat_varid, "units","meter") )
    CALL check( nf90_put_att(ncid, lon_varid, "units", "meter") )

    ! The dimids array is used to pass the dimids of the dimensions of
    ! the netCDF variables. Both of the netCDF variables we are creating
    ! share the same four dimensions. In Fortran, the unlimited
    ! dimension must come last on the list of dimids.
    dimids = (/ lon_dimid, lat_dimid, rec_dimid /)

    ! Define the netCDF variables 
    CALL check( nf90_def_var(ncid, VAR_NAME, NF90_REAL, dimids, varid) )

    ! Assign units attributes to the netCDF variables.
    CALL check( nf90_put_att(ncid, varid, "units",UNITS) )
    CALL assign_common_att(ncid,varid)

    CALL check( nf90_enddef(ncid) )
    CALL check( nf90_put_var(ncid, lat_varid, lats) )
    CALL check( nf90_put_var(ncid, lon_varid, lons) )
  END SUBROUTINE create_nc_v
  SUBROUTINE dump_v(ncid,varid,rec)
    INTEGER ::ncid,varid,rec,start(3),countt(3)
    countt = (/ Nx, Ny+1, 1 /)
    start = (/ 1, 1, 1 /)
    start(3) = rec
    CALL check( nf90_put_var(ncid, varid, TRANSPOSE(vm), start = start, &
         count = countt) )
  END SUBROUTINE dump_v

  SUBROUTINE create_nc_u(ncid,varid)
    IMPLICIT NONE
    INTEGER :: ncid,varid
    CHARACTER (len = *), PARAMETER :: FILE_NAME="u.nc",VAR_NAME="u"
    INTEGER, PARAMETER :: NDIMS = 3
    INTEGER :: start(NDIMS), countt(NDIMS)
    CHARACTER (len = *), PARAMETER :: LAT_NAME = "Y"
    CHARACTER (len = *), PARAMETER :: LON_NAME = "X"
    CHARACTER (len = *), PARAMETER :: REC_NAME = "time"
    INTEGER :: lon_dimid, lat_dimid, rec_dimid
    REAL(kind=rc_kind) ::  lats(NJ+2), lons(NI+2)
    INTEGER :: lon_varid, lat_varid
    INTEGER :: dimids(NDIMS)
    CHARACTER (len = *), PARAMETER :: UNITS = "m/s"
    CHARACTER (len = *), PARAMETER :: LAT_UNITS = "meter"
    CHARACTER (len = *), PARAMETER :: LON_UNITS = "meter"
    lats=et_c(:,1)
    lons=xi_c(1,:)
    ! Create the file. 
    CALL check( nf90_create(casename//'_'//FILE_NAME, nf90_clobber, ncid) )
    CALL check( nf90_def_dim(ncid, LAT_NAME, NJ+2, lat_dimid) )
    CALL check( nf90_def_dim(ncid, LON_NAME, NI+2, lon_dimid) )
    CALL check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )
    CALL check( nf90_def_var(ncid, LAT_NAME, NF90_REAL, lat_dimid, lat_varid) )
    CALL check( nf90_def_var(ncid, LON_NAME, NF90_REAL, lon_dimid, lon_varid) )
    ! Assign units attributes to coordinate variables.
    CALL check( nf90_put_att(ncid, lat_varid, "units","meter") )
    CALL check( nf90_put_att(ncid, lon_varid, "units", "meter") )

    ! The dimids array is used to pass the dimids of the dimensions of
    ! the netCDF variables. Both of the netCDF variables we are creating
    ! share the same four dimensions. In Fortran, the unlimited
    ! dimension must come last on the list of dimids.
    dimids = (/ lon_dimid, lat_dimid, rec_dimid /)

    ! Define the netCDF variables 
    CALL check( nf90_def_var(ncid, VAR_NAME, NF90_REAL, dimids, varid) )

    ! Assign units attributes to the netCDF variables.
    CALL check( nf90_put_att(ncid, varid, "units",UNITS) )
    CALL assign_common_att(ncid,varid)

    CALL check( nf90_enddef(ncid) )
    CALL check( nf90_put_var(ncid, lat_varid, lats) )
    CALL check( nf90_put_var(ncid, lon_varid, lons) )
  END SUBROUTINE create_nc_u
  SUBROUTINE dump_u(ncid,varid,rec)
    INTEGER ::ncid,varid,rec,start(3),countt(3)
    countt = (/ Nx+1, Ny, 1 /)
    start = (/ 1, 1, 1 /)
    start(3) = rec
    CALL check( nf90_put_var(ncid, varid, TRANSPOSE(um), start = start, &
         count = countt) )
  END SUBROUTINE dump_u


  SUBROUTINE assign_common_att(ncid,varid)
    INTEGER :: ncid, varid
    CALL check( nf90_put_att(ncid, varid, "delta_m_bdy", delta_m_bdy) )
    CALL check( nf90_put_att(ncid, varid, "delta_m_int", delta_m_int) )
    CALL check( nf90_put_att(ncid, varid, "delta_s0", delta_s0) )
    CALL check( nf90_put_att(ncid, varid, "r_bdy", r_bdy) )
    CALL check( nf90_put_att(ncid, varid, "r_int", r_int) )
    CALL check( nf90_put_att(ncid, varid, "LX", LX) )    
    CALL check( nf90_put_att(ncid, varid, "WY", WY) )    
    CALL check( nf90_put_att(ncid, varid, "beta", beta) )    
    CALL check( nf90_put_att(ncid, varid, "save_dts", dts) )    
    CALL check( nf90_put_att(ncid, varid, "dt", dt) )    
    CALL check( nf90_put_att(ncid, varid, "we_gyre", we_gyre) ) 
    CALL check( nf90_put_att(ncid, varid, "we_jet", we_jet) ) 
    CALL check( nf90_put_att(ncid, varid, "casename", casename) ) 
CALL check( nf90_put_att(ncid, varid, "L_jet", L_jet) ) 
CALL check( nf90_put_att(ncid, varid, "loc_jet", loc_jet) ) 
CALL check( nf90_put_att(ncid, varid, "q_init", q_init) )
  END SUBROUTINE assign_common_att

  SUBROUTINE pickup(q)
    REAL(kind=rc_kind), DIMENSION(Nx,Ny) :: q_tem
    REAL(kind=rc_kind), DIMENSION(Ny,Nx) :: q
    INTEGER :: pres_varid,ncid
    CALL check( nf90_open(casename//"_pickup.nc", nf90_NoWrite, ncid) )
    !if(status /= nf90_NoErr) call handle_err(status)
    CALL check( nf90_inq_varid(ncid, "q", pres_varid) )
    !if(status /= nf90_NoErr) call handle_err(status)
    CALL check( nf90_get_var(ncid, pres_varid, q_tem) )
    !if(status /= nf90_NoErr) call
    !handle_err(status)
    q=TRANSPOSE(q_tem)
  END SUBROUTINE pickup

  SUBROUTINE close_nc(ncid)
    !  end do 
    CALL check( nf90_close(ncid) )
    ! Close the file. This causes netCDF to flush all buffers and make
    ! sure your data are really written to disk.

    PRINT *,"# *** SUCCESS writing file ", ncid, "!"

  END SUBROUTINE close_nc

  SUBROUTINE check(status)
    !use netcdf
    IMPLICIT NONE
    INTEGER :: status
    IF(status /= nf90_noerr) THEN 
       PRINT *, TRIM(nf90_strerror(status))
       STOP "# Stopped"
    END IF
  END SUBROUTINE check

END MODULE mod_nc
