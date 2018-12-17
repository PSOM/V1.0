subroutine write_cdf_omega(nstp,divQ,wqgeo)
  !subroutine readcdfcn(Tr,uf,vf,wf,xc,yc,zc,xf,yf,zf,nbegin)
  ! ------------------------------------------------------
  ! reads in the netCDF files written by the routine outcdf.f

  !use header, only : NI,NJ,NK,ntr,rc_kind
  use header

  implicit none

#include "netcdf.inc"

  ! include 'dims.f90'
  integer nstp
  REAL(kind=rc_kind), dimension(0:NI+1, 0:NJ+1, 0:NK+1)  :: divQ,wqgeo

  integer :: idvwqgeo,idvdivq

  REAL(kind=rc_kind) ::rcode

  character (len = 550) :: outname_data
  integer :: idOutFile

  integer NI2,NJ2,NK2
  parameter ( NI2=NI+2,NJ2=NJ+2,NK2=NK+2)
  integer start(3), count(3), dims(3)

  DATA start /1, 1, 1/
  DATA count /NI2, NJ2, NK2/

  WRITE(outname_data,'("omega_",I5.5,".cdf")') nstp
  idOutFile =nccre(TRIM(dirout)//outname_data,NCCLOB,rcode)

  dims(1) = ncddef(idOutFile,'x',NI+2,rcode);
  dims(2) = ncddef(idOutFile,'y',NJ+2,rcode);
  dims(3) = ncddef(idOutFile,'sigma',NK+2,rcode);

  idvdivq = ncvdef(idOutFile,'divQ',NCDOUBLE,3,dims,rcode)
  idvwqgeo = ncvdef(idOutFile,'wqgeo',NCDOUBLE,3,dims,rcode)
  CALL ncendf(idOutFile,rcode)

  CALL ncvpt(idOutFile, idvdivq, start, count, divQ, rcode)
  CALL ncvpt(idOutFile, idvwqgeo, start, count, wqgeo, rcode)

  CALL ncclos(idOutFile,rcode)

  PRINT*,"Saving output to:",TRIM(outname_data)

  return
end
