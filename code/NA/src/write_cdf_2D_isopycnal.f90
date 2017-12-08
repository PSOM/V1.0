subroutine write_cdf_2D_isopycnal(sigma,counter_2d,n)
  !     ------------------------------------------------------------------
USE header,ONLY : NI,NJ,NK,ntr,DL,LEN,stressx,T,Tr,rho,xc,yc,zc,zf,vor,strain,shear,freqN2,nconsume,time_seconds,dirout,rc_kind
  !     writes the property T on isopycnal sig                            
#include "netcdf.inc"                                                                        
      integer i,j,k,it,n 
                                                                        
      REAL(kind=rc_kind) :: sig(ntr),sigup,sigdn,dsig       
                                                                        
      REAL(kind=4) ::  xsurf(NI),ysurf(NJ),zsurf(NI,NJ,ntr),vorsurf(NI,NJ,ntr),   &
     &     strainsurf(NI,NJ,ntr),shearsurf(NI,NJ,ntr),                  &
     &     n2surf(NI,NJ,ntr),Tsurf(NI,NJ,ntr),Trsurf(NI,NJ,ntr),stresssurf(NJ),timday     
      character(LEN=150) outname 
                                                                        
      REAL :: sigma,time_days
      INTEGER :: counter_2d

 integer :: iddatfile, idigit, idudx, idudy, idudz, idvbysection, idvby, idvbz, idvcon, idvcy, idvbx
  integer :: idvcz, idvc, idvdivfeddy, idvdivfreyn, idvdx, idvdz, idvd, idvfb, idvh, idvn2bar
  integer :: idvn2, idvnsq100m, idvnsq30m, idvpe,idvpsiv,idvpsiw,idvpv, idvp,idvrhbar,idvrho,idvrnk
  integer :: idvstrain,idvstress,idvstr,idvs,idvtbar,idvtemp,idvtim,idvtr,idvt,idvu,idvvb,idvvc,idvvor,idvv
  integer :: idvwb,idvwc,idvwpv,idvw,idvy,idvzsave,idvz,idwdx,idwdy,idwdz,iimday,ipos
  REAL(kind=rc_kind) ::  rcode
  integer :: idvx,idvshear



                                                                        
      integer start(3),count(3),dims(3),dimstr(4),starttr(4),           &
     &     counttr(4),dimstim,counttim,startwind(2),countwind(2),       &
     &     dimswind(2)                                                  
  

      time_days=time_seconds/86400.d0

     sig=sigma
                                                                        
!      n=0 
!     it=1 is the tracer on the light side of the front-write surface va
!     it=1 
!     k=NK 
!     do j=1,NJ 
!        do i=1,NI 
!           zsurf(i,j,it)= rho(i,j,k) 
!           vorsurf(i,j,it)= vor(i,j,k) 
!           strainsurf(i,j,it)= strain(i,j,k) 
!           shearsurf(i,j,it)= shear(i,j,k) 
!           n2surf(i,j,it)= freqN2(i,j,k) 
!           Tsurf(i,j,it)= T(it,i,j,k,n) 
!        end do 
!     end do 
                                                                        
      do it=1,ntr 
      do j=1,NJ 
         do i=1,NI 
            do k=2,NK 
               if (rho(i,j,k).le.sig(it)) then 
                  dsig= rho(i,j,k-1)-rho(i,j,k) 
                  sigup= (rho(i,j,k-1) - sig(it))/dsig 
                  sigdn= (sig(it)-rho(i,j,k))/dsig 
                  zsurf(i,j,it)= sigup*zc(i,j,k) +sigdn*zc(i,j,k-1) 
                  vorsurf(i,j,it)= sigup*vor(i,j,k) +sigdn*vor(i,j,k-1) 
                  strainsurf(i,j,it)= sigup*strain(i,j,k)               &
     &                 +sigdn*strain(i,j,k-1)                           
                  shearsurf(i,j,it)= sigup*shear(i,j,k)                 &
     &                 +sigdn*shear(i,j,k-1)                            
                  n2surf(i,j,it)= sigup*freqN2(i,j,k)                   &
     &                 +sigdn*freqN2(i,j,k-1)                           
                  Tsurf(i,j,it)= sigup*T(i,j,k,n) +                  &
     &                 sigdn*T(i,j,k-1,n)                            
                  Trsurf(i,j,it)= sigup*Tr(1,i,j,k,n) +                  &
     &                 sigdn*Tr(1,i,j,k-1,n)      
                   goto 101 
               end if 
            end do 
!     If NK reached and rho above not less than sig, then sig is at surf
            zsurf(i,j,it)= zc(i,j,NK) 
            vorsurf(i,j,it)= vor(i,j,NK) 
            strainsurf(i,j,it)= strain(i,j,NK) 
            shearsurf(i,j,it)= shear(i,j,NK) 
            n2surf(i,j,it)= freqN2(i,j,NK) 
            Tsurf(i,j,it)= T(i,j,NK,n) 
            Trsurf(i,j,it)= Tr(1,i,j,NK,n) 
  101    end do 
      end do 
      end do 
                                                                        
      do i=1,NI 
         xsurf(i)= xc(i) 
      end do 
      do j=1,NJ 
         ysurf(j)= yc(j) 
         stresssurf(j)= stressx(j) 
      end do 

     WRITE(outname,'("isopycnal_",I9.9,".cdf")') INT(sigma*1000.)

      if (counter_2d.eq.1) then 
         idDatFile =  nccre(TRIM(dirout)//outname,NCCLOB,rcode) 
                                                                        
         dims(1) = ncddef(idDatFile,'x',NI,rcode) 
         dims(2) = ncddef(idDatFile,'y',NJ,rcode) 
         dims(3) = ncddef(idDatFile,'time',NCUNLIM,rcode) 
                                                                        
         dimswind(1)= dims(2) 
         dimswind(2)= dims(3) 
                                                                        
         dimstr(1) = dims(1) 
         dimstr(2) = dims(2) 
         dimstr(3) = ncddef(idDatFile,'ntr',ntr,rcode) 
         dimstr(4) = dims(3) 
                                                                        
         dimstim= dims(3) 
                                                                        
         idvx = ncvdef(idDatFile,'xc',NCFLOAT,1,dims(1),rcode) 
         idvy = ncvdef(idDatFile,'yc',NCFLOAT,1,dims(2),rcode) 
         idvtim = ncvdef(idDatFile,'day',NCFLOAT,1,dimstim,rcode) 
         idvd = ncvdef(idDatFile,'zc',NCFlOAT,4,dimstr,rcode) 
         idvt = ncvdef(idDatFile,'T',NCFLOAT,4,dimstr,rcode) 
         idvtr = ncvdef(idDatFile,'Tr',NCFLOAT,4,dimstr,rcode) 
         idvvor = ncvdef(idDatFile,'vor',NCFlOAT,4,dimstr,rcode) 
         idvstr = ncvdef(idDatFile,'strain',NCFlOAT,4,dimstr,rcode) 
         idvshear = ncvdef(idDatFile,'shear',NCFlOAT,4,dimstr,rcode) 
         idvn2 = ncvdef(idDatFile,'n2',NCFlOAT,4,dimstr,rcode) 
         idvstress = ncvdef(idDatFile,'stressx',NCFLOAT,2,dimswind,     &
     &        rcode)                                                    
                                                                        
         CALL ncendf(idDatFile,rcode) 
                                                                        
      else 
         idDatFile = ncopn(TRIM(dirout)//outname, NCWRITE,rcode) 
         idvtim = NCVID(idDatFile,'day',RCODE) 
         idvd = NCVID(idDatFile, 'zc', RCODE) 
         idvt = NCVID(idDatFile, 'T', RCODE) 
         idvtr = NCVID(idDatFile, 'Tr', RCODE) 
         idvvor = NCVID(idDatFile, 'vor', RCODE) 
         idvstr = NCVID(idDatFile, 'strain', RCODE) 
         idvshear = NCVID(idDatFile, 'shear', RCODE) 
         idvn2 = NCVID(idDatFile, 'n2', RCODE) 
         idvstress = NCVID(idDatFile, 'stressx', RCODE) 
      endif 
                                                                        
      counttim =1 
                                                                        
      count(1)= NI 
      count(2)= NJ 
      count(3)= 1 
                                                                        
      countwind(1)= NJ 
      countwind(2)= 1 
                                                                        
      counttr(1)= NI 
      counttr(2)= NJ 
      counttr(3)= NTR 
      counttr(4)= 1 
                                                                        
                                                                        
      start(1)= 1 
      start(2)= 1 
      start(3)= counter_2d
                                                                        
      startwind(1)=1 
      startwind(2)= start(3) 
                                                                        
      starttr(1)= 1 
      starttr(2)= 1 
      starttr(3)= 1 
      starttr(4)= counter_2d 
                                                                        
      if (counter_2d.eq.1) then 
         CALL ncvpt(idDatFile,idvx, start(1), count(1), xsurf, rcode) 
         CALL ncvpt(idDatFile,idvy, start(2), count(2), ysurf, rcode) 
      endif 
      CALL ncvpt(idDatFile,idvtim,start(3),counttim,timday,rcode) 
      CALL ncvpt(idDatFile,idvd, starttr, counttr, zsurf, rcode) 
      CALL ncvpt(idDatFile,idvt, starttr, counttr, Tsurf, rcode) 
      CALL ncvpt(idDatFile,idvtr, starttr, counttr, Trsurf, rcode) 
      CALL ncvpt(idDatFile,idvvor, starttr, counttr, vorsurf, rcode) 
      CALL ncvpt(idDatFile,idvstr, starttr, counttr, strainsurf, rcode) 
      CALL ncvpt(idDatFile,idvshear, starttr, counttr, shearsurf, rcode) 
      CALL ncvpt(idDatFile,idvn2, starttr, counttr, n2surf, rcode) 
      CALL ncvpt(idDatFile,idvstress,startwind,countwind,               &
     &     stresssurf,rcode)                                            
                                                                        
      CALL ncclos(idDatFile,rcode) 
                                                                        
      return 
      END                                           
