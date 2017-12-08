subroutine writegeopotential(zd,counter_2d,n)

  !     ------------------------------------------------------------------
USE header,ONLY : dx,dy,NI,NJ,NK,ntr,h,consump,Tr,s,T,rho,u,v,w,p,pvt,vor,strain,freqN2,   &
 & xc,yc,zc,DL,LEN,Jac,nconsume,time_seconds,dirout,rc_kind
  !     periodically writes out the data on the k-th surface              
#include "netcdf.inc"                                                                        
!      INCLUDE '/sw/lib/netcdf-gfortran/include/netcdf.inc'             
                                                                        
      integer i,j,k,l,n,it,keuphotic,kdum,ksurf,k30,k100  
      REAL(kind=rc_kind) :: zd,zmax,zeuphotic,dep,dsig,sigup,sigdn                      
                                                  
      REAL(kind=4) ::  xsurf(NI),ysurf(NJ),zsurf(NI,NJ),ssurf(NI,NJ),              &
     &     tempsurf(NI,NJ),Trsurf(NI,NJ,ntr),rsurf(NI,NJ),rNK(NI,NJ),     &
     &     usurf(NI,NJ),vsurf(NI,NJ),wsurf(NI,NJ),hsurf(NI,NJ),          &
     &psurf(NI,NJ),pvsurf(NI,NJ),         &
     &     vorsurf(NI,NJ),strainsurf(NI,NJ),n2surf(NI,NJ),              &
     &     nsq30m(NI,NJ),nsq100m(NI,NJ),                                &
     &     integcons(NI,NJ,nconsume)                             
!     ,integconsBBL(NI,NJ,ntr)                                          
!                                          
      REAL :: time_days                             
       INTEGER :: counter_2d

 integer :: iddatfile, idigit, idudx, idudy, idudz, idvbysection, idvby, idvbz, idvcon, idvcy, idvbx
  integer :: idvcz, idvc, idvdivfeddy, idvdivfreyn, idvdx, idvdz, idvd, idvfb, idvh, idvn2bar
  integer :: idvn2, idvnsq100m, idvnsq30m, idvpe,idvpsiv,idvpsiw,idvpv, idvp,idvrhbar,idvrho,idvrnk
  integer :: idvstrain,idvstress,idvstr,idvs,idvtbar,idvtp,idvtim,idvtr,idvt,idvu,idvvb,idvvc,idvvor,idvv
  integer :: idvwb,idvwc,idvwpv,idvw,idvy,idvzsave,idvz,idwdx,idwdy,idwdz,ipos
 
  REAL(kind=rc_kind) ::  rcode,vol
  integer :: idvx

                                                                        
      character(LEN=150) :: outname 
                                                                        
      integer start(3),count(3),dims(3),start2d(2),count2d(2),dims2d(2),&
     &     dimstr(4),dimsconstr(4),starttr(4),counttr(4),countconstr(4) &
     &     ,dimstim,counttim                                            
                 
    time_days= time_seconds/86400.d0
      zeuphotic= -100.d0/DL 
      do k=1,NK 
         if (zc(NI/2,NJ/2,k).gt.zeuphotic) then 
            keuphotic = k 
            goto 101 
         end if 
      end do 
! 101  write(6,*) 'keuphotic = ',keuphotic                              
  101 do j=1,NJ 
         do i=1,NI 
            do it=1,nconsume 
               integcons(i,j,it)= 0.0 
!               integconsBBL(i,j,it)= 0.0                               
            end do 
!     sum over euphotic zone (trinit is zero here)                      
!=            do kdum=keuphotic,NK                                      
            do kdum=keuphotic,NK 
               do it=1,nconsume 
                  integcons(i,j,it)= integcons(i,j,it) +0.5d0*          &
     &                 (consump(i,j,kdum,it)+                           &
     &                 abs(consump(i,j,kdum,it)))*                     &
     &                 Jac(i,j,kdum)*DL*LEN*LEN                         
!     integcons is the consump in milimols/timestep if                  
!     consump is in mili-mols/m^3/timestep                              
               end do 
            end do 
!            do kdum=2,3                                                
!               do it=1,ntr                                             
!                  integconsBBL(i,j,it)= integconsBBL(i,j,it) +0.5d0*   
!     &              (consump(i,j,kdum,it)+dabs(consump(i,j,kdum,it)))* 
!     &              Jac(i,j,kdum)*DL*LEN*LEN                           
!     integcons is the consump in milimols/s if consump is in mili-mols/
!               end do                                                  
!     end do                                                            
         end do 
      end do 
 

zsurf=0.;
vorsurf=0.;
strainsurf=0.;
n2surf=0.;
tempsurf=0.;
ssurf=0.;
rsurf=0.;
usurf=0.;
vsurf=0.;
wsurf=0.;
psurf=0.;
pvsurf=0.;
hsurf=0.;
Trsurf=0.;
vorsurf=0.;
                                                                       
      do j=1,NJ 
         do i=1,NI 
            zsurf(i,j)= zd
            
            if (zc(i,j,1)*DL .le. zd) then 
            do k=2,NK 
              if (zc(i,j,k)*DL .ge. zd) then 
                dsig= (zc(i,j,k)-zc(i,j,k-1) )*DL
                sigup= -(zc(i,j,k-1)*DL - zd)/dsig 
                sigdn=  (zc(i,j,k)*DL   - zd)/dsig 
                zsurf(i,j)      = sigup*zc(i,j,k)     + sigdn*zc(i,j,k-1) 
                vorsurf(i,j)    = sigup*vor(i,j,k)    + sigdn*vor(i,j,k-1) 
                strainsurf(i,j) = sigup*strain(i,j,k) + sigdn*strain(i,j,k-1)    
                n2surf(i,j)     = sigup*freqN2(i,j,k) + sigdn*freqN2(i,j,k-1)    
                tempsurf(i,j)   = sigup*T(i,j,k,n)    + sigdn*T(i,j,k-1,n)    
                ssurf(i,j)      = sigup*s(i,j,k,n)    + sigdn*s(i,j,k-1,n)    
                rsurf(i,j)      = sigup*rho(i,j,k)  + sigdn*rho(i,j,k-1)    
                usurf(i,j)      = sigup*u(i,j,k,n)    + sigdn*u(i,j,k-1,n)    
                vsurf(i,j)      = sigup*v(i,j,k,n)    + sigdn*v(i,j,k-1,n)    
                wsurf(i,j)      = sigup*w(i,j,k,n)    + sigdn*w(i,j,k-1,n)    
                psurf(i,j)      = sigup*p(i,j,k)      + sigdn*p(i,j,k-1)    
                pvsurf(i,j)     = sigup*pvt(i,j,k)    + sigdn*pvt(i,j,k-1)    
                pvsurf(i,j)=zsurf(i,j)
                hsurf(i,j)= h(i,j) 
                rNK(i,j)= rho(i,j,NK) 
                do it=1,ntr 
                   Trsurf(i,j,it) = Tr(it,i,j,k,n)     
                end do 
                goto 102 
              end if 
         end do 
              end if 
  102 l=0 
      end do 
      end do 

     vorsurf(1,:)=0.;vorsurf(NI,:)=0.;
     vorsurf(:,1)=0.;vorsurf(:,NJ)=0.;

     do j=2,NJ-1
       do i=2,NI-1
         vorsurf(i,j)=0.5* ( (vsurf(i+1,j)-vsurf(i-1,j))/dx-(usurf(i,j+1)-usurf(i,j-1))/dy )
       enddo
     enddo


      do i=1,NI 
         xsurf(i)= xc(i) 
      end do 
      do j=1,NJ 
         ysurf(j)= yc(j) 
      end do 
                                                                        
!     write average N2 in upper 0-30m and in 30-130 m                   
!     write to a netcdf file                                            
            
k30=0;k100=0;                                                            
      do  k=NK,1,-1 
         dep=-DL*zc(NI/2,NJ/2,k)
         if (dep.gt.30.0) then
            k30=k+1 
            goto 199
         end if
      end do
199   do k=NK,1,-1    
         dep=-DL*zc(NI/2,NJ/2,k)
         if (dep.gt.100.0) then 
            k100= k+1 
            goto 201 
         end if 
      end do 
!201   write(6,*) 'k100',k100,'k30',k30 
201   do j=1,NJ 
         do i=1,NI 
            vol=0.d0 
            nsq30m(i,j)= 0.d0 
            do k=NK,k30,-1 
               nsq30m(i,j) = nsq30m(i,j)+ freqN2(i,j,k)*Jac(i,j,k) 
               vol= vol+ Jac(i,j,k) 
            end do 
            nsq100m(i,j)= nsq30m(i,j) 
            nsq30m(i,j)= nsq30m(i,j)/vol 
            do k=k30+1,k100 
               nsq100m(i,j) = nsq100m(i,j)+ freqN2(i,j,k)*Jac(i,j,k) 
               vol= vol+ Jac(i,j,k) 
            end do 
            nsq100m(i,j)= nsq100m(i,j)/vol 
         end do 
      end do 
                                                                        

     WRITE(outname,'("gslice_",I3.3,".cdf")') INT(-zd)
                                                      
      if (counter_2d.eq.1) then 
         idDatFile =  nccre(TRIM(dirout)//outname,NCCLOB,rcode) 
                                                                        
!         dims2d(1) = ncddef(idDatFile,'xi',NI,rcode)                   
!         dims2d(2) = ncddef(idDatFile,'eta',NJ,rcode)                  
                                                                        
         dims(1) = ncddef(idDatFile,'x',NI,rcode) 
         dims(2) = ncddef(idDatFile,'y',NJ,rcode) 
         dims(3) = ncddef(idDatFile,'time',NCUNLIM,rcode) 
                                                                        
!         dimstr(1) = ncddef(idDatFile,'x',NI,rcode)                    
!         dimstr(2) = ncddef(idDatFile,'y',NJ,rcode)                    
!         dimstr(3) = ncddef(idDatFile,'ntr',ntr,rcode)                 
!         dimstr(4) = ncddef(idDatFile,'time',NCUNLIM,rcode)            
         dimstr(1) = dims(1) 
         dimstr(2) = dims(2) 
         dimstr(3) = ncddef(idDatFile,'ntr',ntr,rcode) 
         dimstr(4) = dims(3) 
                                                                        
         dimsconstr(1) = dims(1) 
         dimsconstr(2) = dims(2) 
         dimsconstr(3) = ncddef(idDatFile,'nconsume',nconsume,rcode) 
         dimsconstr(4) = dims(3) 
                                                                        
         dimstim= dims(3) 
                                                                        
         idvx = ncvdef(idDatFile,'xc',NCFLOAT,1,dims(1),rcode) 
         idvy = ncvdef(idDatFile,'yc',NCFLOAT,1,dims(2),rcode) 
         idvz = ncvdef(idDatFile,'zc',NCFLOAT,2,dims,rcode) 
         idvtim = ncvdef(idDatFile,'day',NCFLOAT,1,dimstim,rcode) 
         idvh = ncvdef(idDatFile,'h',NCFlOAT,3,dims,rcode) 
         idvt = ncvdef(idDatFile,'Tr',NCFLOAT,4,dimstr,rcode) 
         idvs = ncvdef(idDatFile,'s',NCFLOAT,3,dims,rcode) 
         idvtp = ncvdef(idDatFile,'tp',NCFLOAT,3,dims,rcode) 
         idvrho = ncvdef(idDatFile,'rho',NCFLOAT,3,dims,rcode) 
         idvrnk = ncvdef(idDatFile,'rhosurf',NCFLOAT,3,dims,rcode) 
         idvu = ncvdef(idDatFile,'u',NCFLOAT,3,dims,rcode) 
         idvv = ncvdef(idDatFile,'v',NCFLOAT,3,dims,rcode) 
         idvw = ncvdef(idDatFile,'w',NCFLOAT,3,dims,rcode) 
         idvp = ncvdef(idDatFile,'p',NCFLOAT,3,dims,rcode) 
         idvpv = ncvdef(idDatFile,'pv',NCFLOAT,3,dims,rcode) 
         idvvor = ncvdef(idDatFile,'vor',NCFLOAT,3,dims,rcode) 
         idvstrain = ncvdef(idDatFile,'strain',NCFLOAT,3,dims,rcode) 
         idvn2 = ncvdef(idDatFile,'n2',NCFLOAT,3,dims,rcode) 
         idvnsq30m = ncvdef(idDatFile,'nsq30m',NCFLOAT,3,dims,rcode) 
         idvnsq100m = ncvdef(idDatFile,'nsq100m',NCFLOAT,3,dims,rcode) 
         idvc = ncvdef(idDatFile,'integcons',NCFLOAT,4,dimsconstr,rcode) 
!         idvb =ncvdef(idDatFile,'integconsBBL',NCFLOAT,4,dimstr,rcode) 
                                                                        
         CALL ncendf(idDatFile,rcode) 
                                                                        
      else 
         idDatFile = ncopn(TRIM(dirout)//outname, NCWRITE,rcode) 
!         call ncredf(idDatFile,rcode)                                  
!         idvx = NCVID(idDatFile, 'xc', RCODE)                          
!         idvy = NCVID(idDatFile, 'yc', RCODE)                          
!         idvz = NCVID(idDatFile, 'zc', RCODE)                          
         idvtim = NCVID(idDatFile,'day',RCODE) 
         idvh = NCVID(idDatFile, 'h', RCODE) 
         idvt = NCVID(idDatFile, 'Tr', RCODE) 
         idvs = NCVID(idDatFile, 's', RCODE)                          
         idvtp = NCVID(idDatFile, 'tp', RCODE)                          
         idvrho = NCVID(idDatFile, 'rho', RCODE) 
         idvrnk = NCVID(idDatFile, 'rhosurf', RCODE) 
         idvu = NCVID(idDatFile, 'u', RCODE) 
         idvv = NCVID(idDatFile, 'v', RCODE) 
         idvw = NCVID(idDatFile, 'w', RCODE) 
!=         idvp = NCVID(idDatFile, 'p', RCODE)                          
         idvpv = NCVID(idDatFile, 'pv', RCODE)                          
         idvvor = NCVID(idDatFile, 'vor', RCODE) 
         idvstrain = NCVID(idDatFile, 'strain', RCODE) 
         idvn2 = NCVID(idDatFile, 'n2', RCODE) 
         idvnsq30m = NCVID(idDatFile, 'nsq30m', RCODE) 
         idvnsq100m = NCVID(idDatFile, 'nsq100m', RCODE) 
         idvc = NCVID(idDatFile, 'integcons', RCODE) 
!         idvb = NCVID(idDatFile, 'integconsBBL', RCODE)                
      endif 
                                                                        
      counttim =1 
                                                                        
      count2d(1)= NI 
      count2d(2)= NJ 
                                                                        
      count(1)= NI 
      count(2)= NJ 
      count(3)= 1 
                                                                        
      counttr(1)= NI 
      counttr(2)= NJ 
      counttr(3)= NTR 
      counttr(4)= 1 
                                                                        
      countconstr(1)= NI 
      countconstr(2)= NJ 
      countconstr(3)= nconsume 
      countconstr(4)= 1 
                                                                        
      start2d(1)= 1 
      start2d(2)= 1 
                                                                        
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
         CALL ncvpt(idDatFile,idvz, start2d, count2d, zsurf, rcode) 
      endif 
      CALL ncvpt(idDatFile,idvtim,start(3),counttim,time_days,rcode) 
      CALL ncvpt(idDatFile,idvh, start, count, hsurf, rcode) 
      CALL ncvpt(idDatFile,idvt, starttr, counttr, Trsurf, rcode) 
      CALL ncvpt(idDatFile,idvs, start, count, ssurf, rcode)          
      CALL ncvpt(idDatFile,idvtp, start, count, tempsurf, rcode)          
      CALL ncvpt(idDatFile,idvrho, start, count, rsurf, rcode) 
      CALL ncvpt(idDatFile,idvrnk, start, count, rNK, rcode) 
      CALL ncvpt(idDatFile,idvu, start, count, usurf, rcode) 
      CALL ncvpt(idDatFile,idvv, start, count, vsurf, rcode) 
      CALL ncvpt(idDatFile,idvw, start, count, wsurf, rcode) 
!=      CALL ncvpt(idDatFile,idvp, start, count, psurf, rcode)          
      CALL ncvpt(idDatFile,idvpv, start, count, pvsurf, rcode)          
      CALL ncvpt(idDatFile,idvvor, start, count, vorsurf, rcode) 
      CALL ncvpt(idDatFile,idvstrain, start, count, strainsurf, rcode) 
      CALL ncvpt(idDatFile,idvn2, start, count, n2surf, rcode) 
      CALL ncvpt(idDatFile,idvnsq30m, start, count, nsq30m, rcode) 
      CALL ncvpt(idDatFile,idvnsq100m, start, count, nsq100m, rcode) 
      CALL ncvpt(idDatFile,idvc, starttr,countconstr,integcons,rcode) 
!      CALL ncvpt(idDatFile,idvb, starttr,counttr,integconsBBL,rcode)   
                                                                        
      CALL ncclos(idDatFile,rcode) 
                                                                        
      return 
      END                                           
