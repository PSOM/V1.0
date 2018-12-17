subroutine read_cdf_strain(nstp,dudx,dudy,dudz,dvdx,dvdy,dvdz,drdx,drdy)
  !subroutine readcdfcn(Tr,uf,vf,wf,xc,yc,zc,xf,yf,zf,nbegin)
  !     ------------------------------------------------------
  !     reads in the netCDF files written by the routine outcdf.f

  !use header, only : NI,NJ,NK,ntr,rc_kind
  use header

  implicit none

#include "netcdf.inc"
  ! include 'dims.f90'

  integer :: nstp
  integer :: idInFile

  REAL(kind=4) ::  dudx(NI,NJ,NK),dudy(NI,NJ,NK),dudz(NI,NJ,NK), &
                         dvdx(NI,NJ,NK),dvdy(NI,NJ,NK),dvdz(NI,NJ,NK), &
                         drdx(NI,NJ,NK),drdy(NI,NJ,NK)

  integer :: idvdy,idrdy,idrdx
  REAL(kind=rc_kind) ::  rcode
  !
  !double precision Tr(0:NI+1,0:NJ+1,0:NK+1,ntr), &
  !     uf(0:NI,NJ,NK),vf(NI,0:NJ,NK),wf(NI,NJ,0:NK)

  character (len = 550) :: inname_data

  integer start(3), count(3)

  DATA start /1, 1, 1/
  DATA count /NI, NJ, NK/

  WRITE(inname_data,'("strain_",I5.5,".cdf")') nstp
  idInFile = ncopn(TRIM(dirout)//inname_data, NCNOWRIT,rcode)

  idudx = ncvid(idInFile,'udx',rcode)
  idudy = ncvid(idInFile,'udy',rcode)
  idudz = ncvid(idInFile,'udz',rcode)
  idvdx = ncvid(idInFile,'vdx',rcode)
  idvdy = ncvid(idInFile,'vdy',rcode)
  idvdz = ncvid(idInFile,'udz',rcode)
  idrdx = ncvid(idInFile,'rdx',rcode)
  idrdy = ncvid(idInFile,'rdy',rcode)
  call ncvgt( idInFile, idudx, start, count, dudx, rcode )
  call ncvgt( idInFile, idudx, start, count, dudy, rcode )
  call ncvgt( idInFile, idudx, start, count, dudz, rcode )
  call ncvgt( idInFile, idvdx, start, count, dvdx, rcode )
  call ncvgt( idInFile, idvdy, start, count, dvdy, rcode )
  call ncvgt( idInFile, idvdz, start, count, dvdz, rcode )
  call ncvgt( idInFile, idrdx, start, count, drdx, rcode )
  call ncvgt( idInFile, idrdy, start, count, drdy, rcode )
  call ncclos(idInFile, rcode)

  return
end
