subroutine writestrain(step,n)
!---------------------
! writes out the terms for the strainrate tenstor, namely ux,vx,wx, uy,vy,wy,
! uz,vz,wz
  USE header, ONLY : NI,NJ,NK,xc,yc,zc,u,v,w,uf,vf,wf,Jac,rho,UL,LEN,WL,DL,DLinv,gpr,R0,dirout,rc_kind,&
   & ux,uy,vx,vy,wx,wy,wz,time_seconds

#include "netcdf.inc"
  integer dims(3),start(3),count(3),nstp,step

  integer i,idv,idvdy,j,k,n,kp1,km1
  REAL(kind=4) ::  dudx(NI,NJ,NK),dudy(NI,NJ,NK),dudz(NI,NJ,NK),dvdx(NI,NJ,NK),     &
       dvdy(NI,NJ,NK),dvdz(NI,NJ,NK),dwdx(NI,NJ,NK),dwdy(NI,NJ,NK),       &
       dwdz(NI,NJ,NK),rdx(NI,NJ,NK),rdy(NI,NJ,NK),rdz(NI,NJ,NK)
  REAL(kind=rc_kind) ::  LENinv,const,fac

  REAL(kind=rc_kind) ::  rcode
  integer :: idvx
  REAL ::time_days


  character(LEN=150) outname

  DATA start /1, 1, 1/


 integer :: iddatfile, idigit, idudx, idudy, idudz, idvbysection, idvby, idvbz, idvcon, idvcy, idvbx
  integer :: idvcz, idvc, idvdivfeddy, idvdivfreyn, idvdx, idvdz, idvd, idvfb, idvh, idvn2bar
  integer :: idvn2, idvnsq100m, idvnsq30m, idvpe,idvpsiv,idvpsiw,idvpv, idvp,idvrhbar,idvrho,idvrnk
  integer :: idvstrain,idvstress,idvstr,idvs,idvtbar,idvtemp,idvtim,idvtr,idvt,idvu,idvvb,idvvc,idvvor,idvv
  integer :: idvwb,idvwc,idvwpv,idvw,idvy,idvzsave,idvz,idwdx,idwdy,idwdz,ipos

  REAL(kind=rc_kind) ::  vol


    time_days= time_seconds/86400.d0

!     name the output file
! First caculate ux,vx,wx

  do k=1,NK
     do j=1,NJ
        do i=1,NI
          dudx(i,j,k)= (uf(i,j,k)-uf(i-1,j,k))/Jac(i,j,k)
          dvdy(i,j,k)= (vf(i,j,k)-vf(i,j-1,k))/Jac(i,j,k)
          dwdz(i,j,k)= (wf(i,j,k)-wf(i,j,k-1))/Jac(i,j,k)
       end do
    end do
 end do
! Now convert to dimensional form
 dudx=dudx*UL/LEN
 dvdy=dvdy*UL/LEN
 dwdz=dwdz*WL/DL


 !     periodic_e-w boundaries
 do k=0,NK+1
   do j=0,NJ+1
      u(0,j,k,n)= u(NI,j,k,n)
      v(0,j,k,n)= v(NI,j,k,n)
      w(0,j,k,n)= w(NI,j,k,n)
      rho(0,j,k)= rho(NI,j,k)

      u(NI+1,j,k,n)= u(1,j,k,n)
      v(NI+1,j,k,n)= v(1,j,k,n)
      w(NI+1,j,k,n)= w(1,j,k,n) 
      rho(NI+1,j,k)= rho(1,j,k)
   end do
end do

DLinv= 1.0/DL
LENinv=1.0/LEN
const = 2.d0*gpr*10.d0/R0

!     for diff in j (solid bndry)
do k=1,NK
   do i=1,NI
      u(i,0,k,n)= 2.d0*u(i,1,k,n) -u(i,2,k,n)
      v(i,0,k,n)= 2.d0*v(i,1,k,n) -v(i,2,k,n)
      w(i,0,k,n)= 2.d0*w(i,1,k,n) -w(i,2,k,n)
      rho(i,0,k)= 2.d0*rho(i,1,k) -rho(i,2,k)
            
      u(i,NJ+1,k,n)= 2.d0*u(i,NJ,k,n) -u(i,NJ-1,k,n)
      v(i,NJ+1,k,n)= 2.d0*v(i,NJ,k,n) -v(i,NJ-1,k,n)
      w(i,NJ+1,k,n)= 2.d0*w(i,NJ,k,n) -w(i,NJ-1,k,n)
      rho(i,NJ+1,k)= 2.d0*rho(i,NJ,k) -rho(i,NJ-1,k)
   end do
end do


do k=1,NK
   kp1= k+1
   km1= k-1
   if (k.eq.NK) kp1=NK
   if (k.eq.1) km1= 1
   fac= 1.d0/dble(kp1-km1)
   do j=1,NJ
      do i=1,NI
         dvdx(i,j,k)= LENinv*UL*                                  &
              (0.5*(v(i+1,j,k,n)-v(i-1,j,k,n))*ux(i,j) +             &
              0.5*(v(i,j+1,k,n) -v(i,j-1,k,n))*vx(i,j) +             &
              fac*(v(i,j,kp1,n) -v(i,j,km1,n))*wx(i,j,k) )
         dwdx(i,j,k)= LENinv*WL*                                  &
              (0.5*(w(i+1,j,k,n)-w(i-1,j,k,n))*ux(i,j) +             &
              0.5*(w(i,j+1,k,n) -w(i,j-1,k,n))*vx(i,j) +             &
              fac*(w(i,j,kp1,n) -w(i,j,km1,n))*wx(i,j,k) )
         dudy(i,j,k)= LENinv*UL*                                  &
              (0.5*(u(i+1,j,k,n)-u(i-1,j,k,n))*uy(i,j) +             &
              0.5*(u(i,j+1,k,n) -u(i,j-1,k,n))*vy(i,j) +             &
              fac*(u(i,j,kp1,n) -u(i,j,km1,n))*wy(i,j,k) )           
         dwdy(i,j,k)= LENinv*WL*                                  &
              (0.5*(w(i+1,j,k,n)-w(i-1,j,k,n))*uy(i,j) +             &
              0.5*(w(i,j+1,k,n) -w(i,j-1,k,n))*vy(i,j) +             &
              fac*(w(i,j,kp1,n) -w(i,j,km1,n))*wy(i,j,k) )            

         rdx(i,j,k) = LENinv*                                           &
              (0.5*(rho(i+1,j,k)-rho(i-1,j,k))*ux(i,j) +         &
              0.5*(rho(i,j+1,k) -rho(i,j-1,k))*vx(i,j) +         &
              fac*(rho(i,j,kp1) -rho(i,j,km1))*wx(i,j,k) )

         rdy(i,j,k)= LENinv*                                            &
              (0.5*(rho(i+1,j,k)-rho(i-1,j,k))*uy(i,j) +         &
              0.5*(rho(i,j+1,k) -rho(i,j-1,k))*vy(i,j) +         &
              fac*(rho(i,j,kp1) -rho(i,j,km1))*wy(i,j,k) )

         rdz(i,j,k)= fac*(rho(i,j,kp1) -rho(i,j,km1))*wz(i,j,k)*DLinv

         dudz(i,j,k)= UL*fac*(u(i,j,kp1,n) -u(i,j,km1,n))*wz(i,j,k)   &
              *DLinv
         dvdz(i,j,k)= UL*fac*(v(i,j,kp1,n) -v(i,j,km1,n))*wz(i,j,k)   &
              *DLinv

      end do
   end do
end do

!     name the output file

     WRITE(outname,'("strain_",I5.5,".cdf")') step

     count(1)= NI
      count(2)= NJ
      count(3)= NK



idDatFile =  nccre(TRIM(dirout)//outname,NCCLOB,rcode)      

dims(1) = ncddef(idDatFile,'x',NI,rcode)
dims(2) = ncddef(idDatFile,'y',NJ,rcode)
dims(3) = ncddef(idDatFile,'sigma',NK,rcode)

idvx = ncvdef(idDatFile,'xc',NCFLOAT,1,dims(1),rcode)
idvy = ncvdef(idDatFile,'yc',NCFLOAT,1,dims(2),rcode)
idvz = ncvdef(idDatFile,'zc',NCFLOAT,3,dims,rcode)
idudx = ncvdef(idDatFile,'udx',NCFLOAT,3,dims,rcode)
idudy= ncvdef(idDatFile,'udy',NCFLOAT,3,dims,rcode)
idudz = ncvdef(idDatFile,'udz',NCFLOAT,3,dims,rcode)
idvdx = ncvdef(idDatFile,'vdx',NCFLOAT,3,dims,rcode)
idvdy = ncvdef(idDatFile,'vdy',NCFLOAT,3,dims,rcode)
idvdz = ncvdef(idDatFile,'vdz',NCFLOAT,3,dims,rcode)
idwdx = ncvdef(idDatFile,'wdx',NCFLOAT,3,dims,rcode)
idwdy = ncvdef(idDatFile,'wdy',NCFLOAT,3,dims,rcode)
idwdz = ncvdef(idDatFile,'wdz',NCFLOAT,3,dims,rcode)
!idvn = ncvdef(idDatFile,'Nsq',NCFLOAT,3,dims,rcode)
!      idv1 = ncvdef(idDatFile,'vor1',NCFLOAT,3,dims,rcode)
!      idv2 = ncvdef(idDatFile,'vor2',NCFLOAT,3,dims,rcode)
!      idv3 = ncvdef(idDatFile,'vor3',NCFLOAT,3,dims,rcode)
!      idpv1 = ncvdef(idDatFile,'pv1',NCFLOAT,3,dims,rcode)
!      idpv2 = ncvdef(idDatFile,'pv2',NCFLOAT,3,dims,rcode)
!      idpv3 = ncvdef(idDatFile,'pv3',NCFLOAT,3,dims,rcode)


CALL ncendf(idDatFile,rcode)

CALL ncvpt(idDatFile,idvx, start(1), count(1), xc, rcode)
CALL ncvpt(idDatFile,idvy, start(2), count(2), yc, rcode)
CALL ncvpt(idDatFile,idvz, start, count, zc, rcode)

!CALL ncvpt(idDatFile,idrho, start, count, rho, rcode)
!      CALL ncvpt(idDatFile,idu, start, count, u, rcode)
!      CALL ncvpt(idDatFile,idv, start, count, v, rcode)
!      CALL ncvpt(idDatFile,idw, start, count, w, rcode)
CALL ncvpt(idDatFile,idudx, start, count, dudx, rcode)
CALL ncvpt(idDatFile,idudy, start, count, dudy, rcode)
CALL ncvpt(idDatFile,idudz, start, count, dudz, rcode)
CALL ncvpt(idDatFile,idvdx, start, count, dvdx, rcode)
CALL ncvpt(idDatFile,idvdy, start, count, dvdy, rcode)
CALL ncvpt(idDatFile,idvdz, start, count, dvdz, rcode)
CALL ncvpt(idDatFile,idwdx, start, count, dwdx, rcode)
CALL ncvpt(idDatFile,idwdy, start, count, dwdy, rcode)
CALL ncvpt(idDatFile,idwdz, start, count, dwdz, rcode)
!      CALL ncvpt(idDatFile,idvn, start, count, freqN2, rcode)
!     CALL ncvpt(idDatFile,idv1, start, count, vor1, rcode)
!      CALL ncvpt(idDatFile,idv2, start, count, vor2, rcode)
!      CALL ncvpt(idDatFile,idv3, start, count, vor3, rcode)
!      CALL ncvpt(idDatFile,idpv1, start, count, pv1, rcode)
!      CALL ncvpt(idDatFile,idpv2, start, count, pv2, rcode)
!      CALL ncvpt(idDatFile,idpv3, start, count, pv3, rcode)


CALL ncclos(idDatFile,rcode)

return
end subroutine writestrain
