subroutine writexslice(islice,counter_2d,n)

  !     ------------------------------------------------------------      
  USE header,ONLY : NI,NJ,NK,ntr,s,T,rho,u,v,w,consump,Tr,shear,vor,freqN2,yc,zc,fricb,stressx,DL,time_seconds,dirout,rc_kind
  !     reads in a netcdf file and interpolates the data from the sigma gr
  !     onto a z level                                                    
#include "netcdf.inc"                                                                        
!      INCLUDE '/usr/local/include/netcdf.inc' 
!      INCLUDE '/sw/lib/netcdf-gfortran/include/netcdf.inc'             
                                                                        
      integer NI2,NJ2,NK2,i,j,k,n,islice,it,itr 
      parameter ( NI2=NI+2,NJ2=NJ+2,NK2=NK+2) 
      REAL(kind=rc_kind) :: trbar,fbmean 
      REAL(kind=rc_kind) :: svert(NK),Tvert(NK),rvert(NK),uvert(NK),         &
     &     vvert(NK),wvert(NK),zvert(NK),seval                          
                                                                        
      REAL(kind=4) ::  yslice(NJ),zslice(NJ,NK),conslice(ntr,NJ,NK),              &
           sslice(NJ,NK),Tslice(NJ,NK),rslice(NJ,NK),                   &
           Trslice(ntr,NJ,NK),shearslice(NJ,NK),                        &
           uslice(NJ,NK),vslice(NJ,NK),wslice(NJ,NK),vorslice(NJ,NK),   &
           n2slice(NJ,NK),n2barslice(NJ,NK),zcave(NK),fricbsl(NJ,NK),   &
           timday,stressxsp(NJ)  
!     Zonally averaged psi and by,bz                                    
                                                                        
      REAL(kind=rc_kind) ::   psiw(NJ,NK),psiv(NJ,NK),bybar(NJ,NK),          &
     &     bzbar(NJ,NK),vbbar(NJ,NK),wbbar(NJ,NK),psieddy(NJ,NK),       &
     &     rhobar(NJ,NK),n2bar(NJ,NK),drhobardz,drhodz,                 &
     &     wpvbar(NJ,NK),pvbar(NJ,NK),by(NJ,NK)                      
      REAL(kind=rc_kind) :: Feddydiv(NJ,NK),Freyndiv(NJ,NK),                 &
     &     cybar(ntr,NJ,NK),czbar(ntr,NJ,NK),                           &
     &     vcbar(ntr,NJ,NK),wcbar(ntr,NJ,NK),csfeddy(ntr,NJ,NK)         
                                                                        
      REAL(kind=rc_kind) ::  byslice(NJ,NK),bysection(NJ,NK),bzslice(NJ,NK),            &
     &      psivslice(NJ,NK),            &
     &     psiwslice(NJ,NK),vbbarsl(NJ,NK),wbbarsl(NJ,NK),              &
     &     psieddysl(NJ,NK),rhobarslice(NJ,NK),                         &
     &     wpvbarsl(NJ,NK),pvbarsl(NJ,NK),trbarsl(ntr,NJ,NK)            
                                                                        
      REAL(kind=rc_kind) ::  Feddydivsl(NJ,NK),Freyndivsl(NJ,NK),cybarsl(ntr,NJ,NK),    &
     &     czbarsl(ntr,NJ,NK),vcbarsl(ntr,NJ,NK),                       &
     &     wcbarsl(ntr,NJ,NK),csfeddysl(ntr,NJ,NK)                      
                                                                        
      character(LEN=150) outname 
     
      REAL :: time_days
       INTEGER :: counter_2d
                                                                   
      integer start(3),count(3),dims(3),start2d(2),dims4(4),start4(4),  &
           count4(4),dimstim,counttim,count1d,dims1d,dimswind(2),       &
           countwind(2),startwind(2)

 integer :: iddatfile, idigit, idudx, idudy, idudz, idvbysection, idvby, idvbz, idvcon, idvcy, idvbx
  integer :: idvcz, idvc, idvdivfeddy, idvdivfreyn, idvdx, idvdz, idvd, idvfb, idvh, idvn2bar
  integer :: idvn2, idvnsq100m, idvnsq30m, idvpe,idvpsiv,idvpsiw,idvpv, idvp,idvrhbar,idvrho,idvrnk
  integer :: idvstrain,idvstress,idvstr,idvs,idvtbar,idvtemp,idvtim,idvtr,idvt,idvu,idvvb,idvvc,idvvor,idvv
  integer :: idvwb,idvwc,idvwpv,idvw,idvy,idvzsave,idvz,idwdx,idwdy,idwdz,iimday,ipos
  integer :: idvzave,idvshear
  REAL(kind=rc_kind) ::  rcode

     time_days=time_seconds/86400.d0

                                                                          
      call streamfunction(rhobar,n2bar,psiv,psiw,by,bybar,bzbar,vbbar,  &
     &     wbbar,psieddy,wpvbar,pvbar,Feddydiv,Freyndiv,                &
     &     cybar,czbar,vcbar,wcbar)                                     
                                                                        
      do k=1,NK 
         do j=1,NJ 
            fbmean=0.d0 
            do i=1,NI 
               fbmean= fricb(i,j,k)+ fbmean 
            end do 
            fbmean= fbmean/dble(NI) 
            fricbsl(j,k)= fbmean 
         end do 
      end do 
                                                                        
      do k=1,NK 
         do j=1,NJ 
            rhobarslice(j,k)= rhobar(j,k) 
            psivslice(j,k)= psiv(j,k) 
            psiwslice(j,k)= psiw(j,k) 
            byslice(j,k)= bybar(j,k) 
            bysection(j,k)= by(j,k) 
            bzslice(j,k)= bzbar(j,k) 
            vbbarsl(j,k)= vbbar(j,k) 
            wbbarsl(j,k)= wbbar(j,k) 
            psieddysl(j,k)= psieddy(j,k) 
            wpvbarsl(j,k)= wpvbar(j,k) 
            pvbarsl(j,k)= pvbar(j,k) 
                                                                        
            Feddydivsl(j,k)= Feddydiv(j,k) 
            Freyndivsl(j,k)= Freyndiv(j,k) 
            do it=1,ntr 
               cybarsl(it,j,k)= cybar(it,j,k) 
               czbarsl(it,j,k)= czbar(it,j,k) 
               vcbarsl(it,j,k)= vcbar(it,j,k) 
               wcbarsl(it,j,k)= wcbar(it,j,k) 
!               csfeddysl(it,j,k)= csfeddy(it,j,k)                      
            end do 
                                                                        
         end do 
      end do 
                                                                        
                                                                        
      do k=1,NK 
         do j=1,NJ 
            do itr=1,ntr 
               trbar= 0.d0 
               do i=1,NI 
                  trbar= trbar + Tr(itr,i,j,k,n) 
               end do 
               trbar= trbar/dble(NI) 
               trbarsl(itr,j,k)= trbar 
            end do 
         end do 
      end do 
      i= islice 
      do k=1,NK 
         zcave(k)= 0.d0 
         do j=1,NJ 
            zcave(k)= zcave(k)+ zc(i,j,k)*DL 
            zslice(j,k)= zc(i,j,k)*DL 
            yslice(j)= yc(j) 
                                                                        
            sslice(j,k)= s(i,j,k,n)        
            Tslice(j,k)= T(i,j,k,n)
            rslice(j,k)= rho(i,j,k) 
            uslice(j,k)= u(i,j,k,n) 
            vslice(j,k)= v(i,j,k,n) 
            shearslice(j,k)=shear(i,j,k)
            wslice(j,k)= w(i,j,k,n) 
            vorslice(j,k)= vor(i,j,k) 
            do itr=1,ntr 
               conslice(itr,j,k)= consump(i,j,k,itr) 
               Trslice(itr,j,k)= Tr(itr,i,j,k,n) 
            end do 
!            if (k.eq.1) then                                           
!               dz = (zc(i,j,k+1)-zc(i,j,k))*DL                         
!               drhodz= (rho(i,j,k+1)-rho(i,j,k))/dz                    
!               drhobardz= (rhobar(j,k+1)-rhobar(j,k))/dz               
!            else if (k.eq.NK) then                                     
!               dz = (zc(i,j,k)-zc(i,j,k-1))*DL                         
!               drhodz= (rho(i,j,k)-rho(i,j,k-1))/dz                    
!               drhobardz= (rhobar(j,k)-rhobar(j,k-1))/dz               
!            else                                                       
!               dz = (zc(i,j,k+1)-zc(i,j,k-1))*DL                       
!               drhodz= (rho(i,j,k+1)-rho(i,j,k-1))/dz                  
!               drhobardz= (rhobar(j,k+1)-rhobar(j,k-1))/dz             
!            end if                                                     
!            n2slice(j,k) = -drhodz*9.81/R0                             
!            n2barslice(j,k) = -drhobardz*9.81/R0                       
            n2slice(j,k) = freqN2(i,j,k) 
            n2barslice(j,k) = n2bar(j,k) 
         end do 
         zcave(k)= zcave(k)/dble(NJ) 
      end do 
      do j=1,NJ
         stressxsp(j)= stressx(j)
      end do
!                                                                       
                                                                        
!     write to a netcdf file                  
      WRITE(outname,'("xslice_",I3.3,".cdf")') islice
                                                                        
                                                                        
      if (counter_2d.eq.1) then 
         idDatFile =  nccre(TRIM(dirout)//outname,NCCLOB,rcode) 
                                                                        
         dims(1) = ncddef(idDatFile,'y',NJ,rcode) 
         dims(2) = ncddef(idDatFile,'z',NK,rcode) 
         dims(3) = ncddef(idDatFile,'time',NCUNLIM,rcode) 
                                                                        
         dims1d= dims(2) 
                                                                        
         dimswind(1)= dims(1)
         dimswind(2)= dims(3)

         dims4(1) = ncddef(idDatFile,'ntr',ntr,rcode) 
         dims4(2) = dims(1) 
         dims4(3) = dims(2) 
         dims4(4) = dims(3) 
                                                                        
         dimstim= dims(3) 
                                                                        
         idvy = ncvdef(idDatFile,'yc',NCFLOAT,1,dims(1),rcode) 
         idvzave = ncvdef(idDatFile,'zcave',NCFLOAT,1,dims1d,rcode) 
         idvz = ncvdef(idDatFile,'zc',NCFLOAT,2,dims,rcode) 
         idvtim = ncvdef(idDatFile,'day',NCFLOAT,1,dimstim,rcode) 
         idvstress = ncvdef(idDatFile,'stressx',NCFLOAT,2,dimswind,rcode) 
         idvcon = ncvdef(idDatFile,'consump',NCFLOAT,4,dims4,rcode) 
         idvtr = ncvdef(idDatFile,'tr',NCFLOAT,4,dims4,rcode) 
         idvtbar = ncvdef(idDatFile,'trbar',NCFLOAT,4,dims4,rcode) 
         idvs = ncvdef(idDatFile,'s',NCFLOAT,3,dims,rcode)          
         idvt = ncvdef(idDatFile,'temp',NCFLOAT,3,dims,rcode)          
         idvrho = ncvdef(idDatFile,'rho',NCFLOAT,3,dims,rcode) 
         idvu = ncvdef(idDatFile,'u',NCFLOAT,3,dims,rcode) 
         idvv = ncvdef(idDatFile,'v',NCFLOAT,3,dims,rcode) 
         idvw = ncvdef(idDatFile,'w',NCFLOAT,3,dims,rcode) 
         idvvor = ncvdef(idDatFile,'vor',NCFLOAT,3,dims,rcode) 
         idvshear = ncvdef(idDatFile,'shear',NCFLOAT,3,dims,rcode) 
         idvn2 = ncvdef(idDatFile,'n2',NCFLOAT,3,dims,rcode) 
         idvfb = ncvdef(idDatFile,'fricb',NCFLOAT,3,dims,rcode) 
                                                                        
         idvrhbar = ncvdef(idDatFile,'rhobar',NCFLOAT,3,dims,rcode) 
         idvn2bar = ncvdef(idDatFile,'n2bar',NCFLOAT,3,dims,rcode) 
         idvpsiv = ncvdef(idDatFile,'psivbar',NCFLOAT,3,dims,rcode) 
         idvpsiw = ncvdef(idDatFile,'psiwbar',NCFLOAT,3,dims,rcode) 
         idvbysection = ncvdef(idDatFile,'by',NCFLOAT,3,dims,rcode) 
         idvby = ncvdef(idDatFile,'bybar',NCFLOAT,3,dims,rcode) 
         idvbz = ncvdef(idDatFile,'bzbar',NCFLOAT,3,dims,rcode) 
         idvvb = ncvdef(idDatFile,'vbbar',NCFLOAT,3,dims,rcode) 
         idvwb = ncvdef(idDatFile,'wbbar',NCFLOAT,3,dims,rcode) 
         idvpe = ncvdef(idDatFile,'psieddy',NCFLOAT,3,dims,rcode) 
         idvpv = ncvdef(idDatFile,'pvbar',NCFLOAT,3,dims,rcode) 
         idvwpv = ncvdef(idDatFile,'wpvbar',NCFLOAT,3,dims,rcode) 
                                                                        
         idvdivfreyn=ncvdef(idDatFile,'divfreyn',NCFLOAT,3,dims,rcode) 
         idvdivfeddy=ncvdef(idDatFile,'divfeddy',NCFLOAT,3,dims,rcode) 
                                                                        
         idvcy = ncvdef(idDatFile,'cybar',NCFLOAT,4,dims4,rcode) 
         idvcz = ncvdef(idDatFile,'czbar',NCFLOAT,4,dims4,rcode) 
         idvvc = ncvdef(idDatFile,'vcbar',NCFLOAT,4,dims4,rcode) 
         idvwc = ncvdef(idDatFile,'wcbar',NCFLOAT,4,dims4,rcode) 
!         idvpec = ncvdef(idDatFile,'psieddytracer',NCFLOAT,4,          
!     &        dims4,rcode)                                             
                                                                        
                                                                        
         CALL ncendf(idDatFile,rcode) 
                                                                        
      else 
         idDatFile = ncopn(TRIM(dirout)//outname, NCWRITE,rcode) 
!         call ncredf(idDatFile,rcode)                                  
         idvtim = NCVID(idDatFile,'day',RCODE) 
         idvstress = NCVID(idDatFile, 'stressx', RCODE) 
         idvcon = NCVID(idDatFile, 'consump', RCODE) 
         idvtr = NCVID(idDatFile, 'tr', RCODE) 
         idvtbar = NCVID(idDatFile, 'trbar', RCODE) 
         idvs = NCVID(idDatFile, 's', RCODE)                       
         idvt = NCVID(idDatFile, 'temp', RCODE)                       
         idvrho = NCVID(idDatFile, 'rho', RCODE) 
         idvu = NCVID(idDatFile, 'u', RCODE) 
         idvv = NCVID(idDatFile, 'v', RCODE) 
         idvw = NCVID(idDatFile, 'w', RCODE) 
         idvvor = NCVID(idDatFile, 'vor', RCODE) 
         idvshear = NCVID(idDatFile, 'shear', RCODE) 
         idvn2 = NCVID(idDatFile, 'n2', RCODE) 
         idvfb = NCVID(idDatFile, 'fricb', RCODE) 
                                                                        
         idvrhbar = NCVID(idDatFile, 'rhobar', RCODE) 
         idvn2bar = NCVID(idDatFile, 'n2bar', RCODE) 
         idvpsiv = NCVID(idDatFile, 'psivbar', RCODE) 
         idvpsiw = NCVID(idDatFile, 'psiwbar', RCODE) 
         idvbysection = NCVID(idDatFile, 'by', RCODE) 
         idvby = NCVID(idDatFile, 'bybar', RCODE) 
         idvbz = NCVID(idDatFile, 'bzbar', RCODE) 
         idvvb = NCVID(idDatFile, 'vbbar', RCODE) 
         idvwb = NCVID(idDatFile, 'wbbar', RCODE) 
         idvpe = NCVID(idDatFile, 'psieddy', RCODE) 
         idvpv = NCVID(idDatFile, 'pvbar', RCODE) 
         idvwpv = NCVID(idDatFile, 'wpvbar', RCODE) 
                                                                        
         idvdivfreyn = ncvid(idDatFile,'divfreyn',rcode) 
         idvdivfeddy = ncvid(idDatFile,'divfeddy',rcode) 
                                                                        
         idvcy = NCVID(idDatFile,'cybar',rcode) 
         idvcz = NCVID(idDatFile,'czbar',rcode) 
         idvvc = NCVID(idDatFile,'vcbar',rcode) 
         idvwc = NCVID(idDatFile,'wcbar',rcode) 
!         idvpec = NCVID(idDatFile,'psieddytracer',rcode)               
                                                                        
      endif 
                                                                        
      counttim =1 
      count1d= NK 
                                                                        
      count(1)= NJ 
      count(2)= NK 
      count(3)= 1 
                                                                        
      countwind(1)= NJ
      countwind(2)= 1

      count4(1)= ntr 
      count4(2)= NJ 
      count4(3)= NK 
      count4(4)= 1 
                                                                        
      start2d(1)= 1 
      start2d(2)= 1 
                                                                        
      start(1)= 1 
      start(2)= 1 
      start(3)= counter_2d 

      startwind(1)= 1
      startwind(2)= counter_2d

                                                                        
      start4(1)= 1 
      start4(2)= 1 
      start4(3)= 1 
      start4(4)= counter_2d
                                                                        
      if (counter_2d.eq.1) then 
         CALL ncvpt(idDatFile,idvy, start(1), count(1), yslice, rcode) 
         CALL ncvpt(idDatFile,idvzave, start(1),count1d,zcave, rcode) 
         CALL ncvpt(idDatFile,idvz, start, count, zslice, rcode) 
      endif 
      CALL ncvpt(idDatFile,idvtim,start(3),counttim,timday,rcode) 
      CALL ncvpt(idDatFile,idvstress, startwind, countwind, stressxsp, rcode) 
      CALL ncvpt(idDatFile,idvcon, start4, count4, conslice, rcode) 
      CALL ncvpt(idDatFile,idvtr, start4, count4, Trslice, rcode) 
      CALL ncvpt(idDatFile,idvtbar, start4, count4, Trbarsl, rcode) 
      CALL ncvpt(idDatFile,idvs, start, count, sslice, rcode) 
      CALL ncvpt(idDatFile,idvt, start, count, Tslice, rcode) 
      CALL ncvpt(idDatFile,idvrho, start, count, rslice, rcode) 
      CALL ncvpt(idDatFile,idvu, start, count, uslice, rcode) 
      CALL ncvpt(idDatFile,idvv, start, count, vslice, rcode) 
      CALL ncvpt(idDatFile,idvw, start, count, wslice, rcode) 
      CALL ncvpt(idDatFile,idvvor, start, count, vorslice, rcode) 
      CALL ncvpt(idDatFile,idvshear, start, count, shearslice, rcode) 
      CALL ncvpt(idDatFile,idvn2, start, count, n2slice, rcode) 
      CALL ncvpt(idDatFile,idvfb, start, count, fricbsl, rcode) 
                                                                        
      CALL ncvpt(idDatFile,idvrhbar, start, count, rhobarslice, rcode) 
      CALL ncvpt(idDatFile,idvn2bar, start, count, n2barslice, rcode) 
      CALL ncvpt(idDatFile,idvpsiv, start, count, psivslice, rcode) 
      CALL ncvpt(idDatFile,idvpsiw, start, count, psiwslice, rcode) 
      CALL ncvpt(idDatFile,idvbysection, start, count, bysection, rcode) 
      CALL ncvpt(idDatFile,idvby, start, count, byslice, rcode) 
      CALL ncvpt(idDatFile,idvbz, start, count, bzslice, rcode) 
      CALL ncvpt(idDatFile,idvvb, start, count, vbbarsl, rcode) 
      CALL ncvpt(idDatFile,idvwb, start, count, wbbarsl, rcode) 
      CALL ncvpt(idDatFile,idvpe, start, count, psieddysl, rcode) 
      CALL ncvpt(idDatFile,idvpv, start, count, pvbarsl, rcode) 
      CALL ncvpt(idDatFile,idvwpv, start, count, wpvbarsl, rcode) 
                                                                        
      CALL ncvpt(idDatFile,idvdivfreyn, start, count, Freyndivsl,rcode) 
      CALL ncvpt(idDatFile,idvdivfeddy, start, count, Feddydivsl,rcode) 
                                                                        
      CALL ncvpt(idDatFile,idvcy, start4, count4, cybarsl, rcode) 
      CALL ncvpt(idDatFile,idvcz, start4, count4, czbarsl, rcode) 
      CALL ncvpt(idDatFile,idvvc, start4, count4, vcbarsl, rcode) 
      CALL ncvpt(idDatFile,idvwc, start4, count4, wcbarsl, rcode) 
!      CALL ncvpt(idDatFile,idvpec, start4, count4, csfeddysl, rcode)   
      CALL ncclos(idDatFile,rcode) 
                                                                        
      return 
      END                                           
