subroutine write_cdf_2D_isopycnal(sigma,counter_2d,n)
  !     ------------------------------------------------------------------
USE header,ONLY : NI,NJ,NK,ntr,DL,LEN,stress_top_x,s,T,rho,u,v,Tr,consump,xc,yc,zc,zf, & 
     vor,strain,shear,freqN2,nconsume,time_seconds,time_days,dirout,rc_kind
  !     writes the property T on isopycnal sig                            
#include "netcdf.inc"                                                                        
      integer i,j,k,it,n 
                                                                        
      REAL(kind=rc_kind) :: sig,sigup,sigdn,dsig       
                                                                        
      REAL(kind=4) ::  xsurf(NI),ysurf(NJ),zsurf(NI,NJ),vorsurf(NI,NJ),   &     !shearsurf(NI,NJ),
     &     strainsurf(NI,NJ),                  &
     &     n2surf(NI,NJ),Tsurf(NI,NJ),ssurf(NI,NJ),usurf(NI,NJ),vsurf(NI,NJ),Trsurf(NI,NJ,ntr),consurf(NI,NJ),timday     
      character(LEN=150) outname 
                                                                        
      REAL :: sigma
      INTEGER :: counter_2d

 integer :: iddatfile, idigit, idudx, idudy, idudz, idvbysection, idvby, idvbz, idvcon, idvcy, idvbx
  integer :: idvcz, idvc, idvdivfeddy, idvdivfreyn, idvdx, idvdz, idvd, idvfb, idvh, idvn2bar
  integer :: idvn2, idvnsq100m, idvnsq30m, idvpe,idvpsiv,idvpsiw,idvpv, idvp,idvrhbar,idvrho,idvrnk
  integer :: idvstrain,idvstress,idvstr,idvtbar,idvtemp,idvtim,idvtr,idvt,idvs,idvu,idvvb,idvvc,idvvor,idvv
  integer :: idvwb,idvwc,idvwpv,idvw,idvy,idvzsave,idvz,idwdx,idwdy,idwdz,iimday,ipos
  REAL(kind=rc_kind) ::  rcode
  integer :: idvx,idvshear



                                                                        
      integer start(3),count(3),dims(3),dimstr(4),starttr(4),           &
     &     counttr(4),dimstim,counttim,iconsump
  
      timday= time_days
      iconsump=1

     sig=sigma+1000.0
                                                                        
      do j=1,NJ 
         do i=1,NI 
            do k=2,NK 
               if (rho(i,j,k).le.sig) then 
                  dsig= rho(i,j,k-1)-rho(i,j,k) 
                  sigup= (rho(i,j,k-1) - sig)/dsig 
                  sigdn= (sig-rho(i,j,k))/dsig 
                  zsurf(i,j)= sigup*zc(i,j,k) +sigdn*zc(i,j,k-1) 
                  vorsurf(i,j)= sigup*vor(i,j,k) +sigdn*vor(i,j,k-1) 
                  strainsurf(i,j)= sigup*strain(i,j,k)               &
     &                 +sigdn*strain(i,j,k-1)                           
!                  shearsurf(i,j)= sigup*shear(i,j,k)                 &
!     &                 +sigdn*shear(i,j,k-1)                            
                  n2surf(i,j)= sigup*freqN2(i,j,k)                   &
     &                 +sigdn*freqN2(i,j,k-1)                           
                  Tsurf(i,j)= sigup*T(i,j,k,n) +                  &
     &                 sigdn*T(i,j,k-1,n)                            
                  ssurf(i,j)= sigup*s(i,j,k,n) +                  &
     &                 sigdn*s(i,j,k-1,n)                            
                  usurf(i,j)= sigup*u(i,j,k,n) +                  &
     &                 sigdn*u(i,j,k-1,n)                            
                  vsurf(i,j)= sigup*v(i,j,k,n) +                  &
     &                 sigdn*v(i,j,k-1,n)                            
                  consurf(i,j)= sigup*consump(i,j,k,iconsump) +                  &
     &                 sigdn*consump(i,j,k-1,iconsump)                            
                  do it=1,ntr 
                     Trsurf(i,j,it)= sigup*Tr(it,i,j,k,n) +                  &
                          sigdn*Tr(it,i,j,k-1,n) 
                  end do
                  goto 101 
               end if 
            end do 
!     If NK reached and rho above not less than sig, then sig is at surf
            zsurf(i,j)= zc(i,j,NK) 
            vorsurf(i,j)= vor(i,j,NK) 
            strainsurf(i,j)= strain(i,j,NK) 
!            shearsurf(i,j)= shear(i,j,NK) 
            n2surf(i,j)= freqN2(i,j,NK) 
            ssurf(i,j)= s(i,j,NK,n) 
            Tsurf(i,j)= T(i,j,NK,n) 
            usurf(i,j)= u(i,j,NK,n) 
            vsurf(i,j)= v(i,j,NK,n) 
            consurf(i,j)= consump(i,j,NK,iconsump) 
            do it=1,ntr
               Trsurf(i,j,it)= Tr(it,i,j,NK,n) 
            end do

101      end do
      end do
      zsurf= zsurf*1000.
      do i=1,NI 
         xsurf(i)= xc(i) 
      end do 
      do j=1,NJ 
         ysurf(j)= yc(j) 
      end do 

     WRITE(outname,'("isopycnal_",I9.9,".cdf")') INT(sigma*1000.)

      if (counter_2d.eq.1) then 
         idDatFile =  nccre(TRIM(dirout)//outname,NCCLOB,rcode) 
                                                                        
         dims(1) = ncddef(idDatFile,'x',NI,rcode) 
         dims(2) = ncddef(idDatFile,'y',NJ,rcode) 
         dims(3) = ncddef(idDatFile,'time',NCUNLIM,rcode) 
                                                                        
         dimstr(1) = dims(1) 
         dimstr(2) = dims(2) 
         dimstr(3) = ncddef(idDatFile,'ntr',ntr,rcode) 
         dimstr(4) = dims(3) 
                                                                        
         dimstim= dims(3) 
                                                                        
         idvx = ncvdef(idDatFile,'xc',NCFLOAT,1,dims(1),rcode) 
         idvy = ncvdef(idDatFile,'yc',NCFLOAT,1,dims(2),rcode) 
         idvtim = ncvdef(idDatFile,'day',NCFLOAT,1,dimstim,rcode) 
         idvd = ncvdef(idDatFile,'zc',NCFlOAT,3,dims,rcode) 
         idvt = ncvdef(idDatFile,'temp',NCFLOAT,3,dims,rcode) 
         idvs = ncvdef(idDatFile,'s',NCFLOAT,3,dims,rcode) 
         idvu = ncvdef(idDatFile,'u',NCFLOAT,3,dims,rcode) 
         idvv = ncvdef(idDatFile,'v',NCFLOAT,3,dims,rcode) 
         idvcon = ncvdef(idDatFile,'consump',NCFLOAT,3,dims,rcode) 
         idvtr = ncvdef(idDatFile,'Tr',NCFLOAT,4,dimstr,rcode) 
         idvvor = ncvdef(idDatFile,'vor',NCFlOAT,3,dims,rcode) 
         idvstr = ncvdef(idDatFile,'strain',NCFlOAT,3,dims,rcode) 
!         idvshear = ncvdef(idDatFile,'shear',NCFlOAT,4,dimstr,rcode) 
         idvn2 = ncvdef(idDatFile,'n2',NCFlOAT,3,dims,rcode) 
                                                                        
         CALL ncendf(idDatFile,rcode) 
                                                                        
      else 
         idDatFile = ncopn(TRIM(dirout)//outname, NCWRITE,rcode) 
         idvtim = NCVID(idDatFile,'day',RCODE) 
         idvd = NCVID(idDatFile, 'zc', RCODE) 
         idvt = NCVID(idDatFile, 'temp', RCODE) 
         idvs = NCVID(idDatFile, 's', RCODE) 
         idvu = NCVID(idDatFile, 'u', RCODE) 
         idvv = NCVID(idDatFile, 'v', RCODE) 
         idvcon = NCVID(idDatFile, 'consump', RCODE) 
         idvtr = NCVID(idDatFile, 'Tr', RCODE) 
         idvvor = NCVID(idDatFile, 'vor', RCODE) 
         idvstr = NCVID(idDatFile, 'strain', RCODE) 
!         idvshear = NCVID(idDatFile, 'shear', RCODE) 
         idvn2 = NCVID(idDatFile, 'n2', RCODE) 
      endif 
                                                                        
      counttim =1 
                                                                        
      count(1)= NI 
      count(2)= NJ 
      count(3)= 1 
                                                                        
      counttr(1)= NI 
      counttr(2)= NJ 
      counttr(3)= NTR 
      counttr(4)= 1 
                                                                        
                                                                        
      start(1)= 1 
      start(2)= 1 
      start(3)= counter_2d
                                                                        
                                                                        
      starttr(1)= 1 
      starttr(2)= 1 
      starttr(3)= 1 
      starttr(4)= counter_2d 
                                                                        
      if (counter_2d.eq.1) then 
         CALL ncvpt(idDatFile,idvx, start(1), count(1), xsurf, rcode) 
         CALL ncvpt(idDatFile,idvy, start(2), count(2), ysurf, rcode) 
      endif 
      CALL ncvpt(idDatFile,idvtim,start(3),counttim,timday,rcode) 
      CALL ncvpt(idDatFile,idvd, start, count, zsurf, rcode) 
      CALL ncvpt(idDatFile,idvt, start, count, Tsurf, rcode) 
      CALL ncvpt(idDatFile,idvs, start, count, ssurf, rcode) 
      CALL ncvpt(idDatFile,idvu, start, count, usurf, rcode) 
      CALL ncvpt(idDatFile,idvv, start, count, vsurf, rcode) 
      CALL ncvpt(idDatFile,idvcon, start, count, consurf, rcode) 
      CALL ncvpt(idDatFile,idvtr, starttr, counttr, Trsurf, rcode) 
      CALL ncvpt(idDatFile,idvvor, start, count, vorsurf, rcode) 
      CALL ncvpt(idDatFile,idvstr, start, count, strainsurf, rcode) 
!      CALL ncvpt(idDatFile,idvshear, starttr, counttr, shearsurf, rcode) 
      CALL ncvpt(idDatFile,idvn2, start, count, n2surf, rcode) 
                                                                        
      CALL ncclos(idDatFile,rcode) 
                                                                        
      return 
      END                                           
