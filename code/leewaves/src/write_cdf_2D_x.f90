subroutine write_cdf_2D_x(islice,counter_2d,n)

   USE header,ONLY : NI,NJ,NK,ntr,nconsume,s,T,rho,u,v,w,consump,Tr,            &
  &    KzMom,KzTr,shear,vor,freqN2,yc,zc,fricb,stress_top_x,DL,                 &
  &    time_seconds,time_days,dirout,rc_kind,pvt,rp,stress_top_y,pv3,swr,qloss, &
  &    T_ini,s_ini,T_restore,s_restore                    !!! ADDED
  !     reads in a netcdf file and interpolates the data from the sigma gr
  !     onto a z level                                                    
#include "netcdf.inc"                                                                        
!    INCLUDE '/usr/local/include/netcdf.inc'
!    INCLUDE '/sw/lib/netcdf-gfortran/include/netcdf.inc'
                                                                        
   integer NI2,NJ2,NK2,i,j,k,n,islice,it,itr 
   parameter(NI2=NI+2,NJ2=NJ+2,NK2=NK+2) 
   REAL(kind=rc_kind) :: trbar,fbmean 
   REAL(kind=rc_kind) :: svert(NK),Tvert(NK),rvert(NK),uvert(NK),vvert(NK),wvert(NK),zvert(NK),seval

   REAL(kind=4) :: yslice(NJ),zslice(NJ,NK),conslice(nconsume,NJ,NK),                                    &
  &     sslice(NJ,NK),Tslice(NJ,NK),rslice(NJ,NK),                                                       &
  &     s_ini_no_ghost(NJ,NK),T_ini_no_ghost(NJ,NK),s_restore_no_ghost(NJ,NK),T_restore_no_ghost(NJ,NK), &
  &     Trslice(ntr,NJ,NK),shearslice(NJ,NK),                                                            &
  &     uslice (NJ,NK),vslice (NJ,NK),wslice(NJ,NK),vorslice(NJ,NK),                                     &
  &     ugslice(NJ,NK),vgslice(NJ,NK),                                                                   &
  &     n2slice(NJ,NK),n2barslice(NJ,NK),zcave(NK),fricbsl(NJ,NK),                                       &
  &     stressxsp(NJ),stressysp(NJ),pvtslice(NJ,NK),pv3slice(NJ,NK),pslice(NJ,NK),swrsl(NJ),qlosssl(NJ)
        
!    Zonally averaged psi and by,bz

   REAL(kind=rc_kind) :: psiw(NJ,NK),psiv(NJ,NK),bybar(NJ,NK),       &
  &     bzbar(NJ,NK),vbbar(NJ,NK),wbbar(NJ,NK),psieddy(NJ,NK),       &
  &     rhobar(NJ,NK),n2bar(NJ,NK),drhobardz,drhodz,                 &
  &     wpvbar(NJ,NK),pvbar(NJ,NK),by(NJ,NK)
  
   REAL(kind=rc_kind) :: Feddydiv(NJ,NK),Freyndiv(NJ,NK),            &
  &     cybar(ntr,NJ,NK),czbar(ntr,NJ,NK),                           &
  &     vcbar(ntr,NJ,NK),wcbar(ntr,NJ,NK),csfeddy(ntr,NJ,NK)

   REAL(kind=4) :: byslice(NJ,NK),bysection(NJ,NK),bzslice(NJ,NK),   &
  &     psivslice(NJ,NK),                                            &
  &     psiwslice(NJ,NK),vbbarsl(NJ,NK),wbbarsl(NJ,NK),              &
  &     psieddysl(NJ,NK),rhobarslice(NJ,NK),                         &
  &     wpvbarsl(NJ,NK),pvbarsl(NJ,NK),trbarsl(ntr,NJ,NK),           &
  &     KzTrsl(NJ,NK),KzMomsl(NJ,NK)

   REAL(kind=rc_kind) :: Feddydivsl(NJ,NK),Freyndivsl(NJ,NK),cybarsl(ntr,NJ,NK),    &
  &     czbarsl(ntr,NJ,NK),vcbarsl(ntr,NJ,NK),                                      &
  &     wcbarsl(ntr,NJ,NK),csfeddysl(ntr,NJ,NK)               

   REAL(kind=rc_kind), dimension(0:NI+1,0:NJ+1,0:NK+1) :: ugeo,vgeo

   character(LEN=150) outname 

   INTEGER :: counter_2d

   integer start(4),count(4),dims(4),start2d(2),dims4(4),dims4c(4),start4(4),      &
           count4(4),count4c(4),dimstim,counttim,count1d,dims1d,dimswind(3),       &
           countwind(3),startwind(3)

   integer :: iddatfile,idigit,idudx,idudy,dudz,idvbysection,idvby,idvbz,idvcon,idvcy,idvbx
   integer :: idvcz,idvc,idvdivfeddy,idvdivfreyn,idvdx,idvdz,idvd,idvfb,idvh,idvn2bar
   integer :: idvn2,idvnsq100m,idvnsq30m,idvpe,idvpsiv,idvpsiw,idvpv,idvp,idvrhbar,idvrho,idvrnk
   integer :: idvstrain,idvstressx,idvstressy,idvstr,idvs,idvtbar,idvtemp,&
  &           idvtim,idvtr,idvt,idvu,idvvb,idvvc,idvvor,idvv,idvug,idvvg
   integer :: idvwb,idvwc,idvwpv,idvw,idvy,idvzsave,idvz
   integer :: idvzave,idvshear,idvpvt,idvpv3
!    integer :: idvtini,idvsini,idvtrestore,idvsrestore   !!! ADDED
   REAL(kind=rc_kind) :: rcode

! =================================================================================================
!                                       Initialization
! =================================================================================================
!    call diag_streamfunction(rhobar,n2bar,psiv,psiw,by,bybar,bzbar,vbbar,  &
!   &     wbbar,psieddy,wpvbar,pvbar,Feddydiv,Freyndiv,                     &
!   &     cybar,czbar,vcbar,wcbar)                                  
!    call uv_geostrophic(ugeo,vgeo)
!
!    do k=1,NK
!       do j=1,NJ
!          fbmean=0.d0
!          do i=1,NI
!             fbmean= fricb(i,j,k)+ fbmean
!          end do
!          fbmean= fbmean/dble(NI)
!          fricbsl(j,k)= fbmean
!       end do
!    end do

! =================================================================================================
!                                         2D slices
! ================================================================================================= 
!    do k=1,NK
!       do j=1,NJ
!          rhobarslice(j,k)= rhobar(j,k)
!          psivslice(j,k)= psiv(j,k)
!          psiwslice(j,k)= psiw(j,k)
!          byslice(j,k)= bybar(j,k)
!          bysection(j,k)= by(j,k)
!          bzslice(j,k)= bzbar(j,k)
!          vbbarsl(j,k)= vbbar(j,k)
!          wbbarsl(j,k)= wbbar(j,k)
!          psieddysl(j,k)= psieddy(j,k)
!          wpvbarsl(j,k)= wpvbar(j,k)
!          pvbarsl(j,k)= pvbar(j,k)
!          Feddydivsl(j,k)= Feddydiv(j,k)
!          Freyndivsl(j,k)= Freyndiv(j,k)
!          do it=1,ntr
!             cybarsl(it,j,k)= cybar(it,j,k)
!             czbarsl(it,j,k)= czbar(it,j,k)
!             vcbarsl(it,j,k)= vcbar(it,j,k)
!             wcbarsl(it,j,k)= wcbar(it,j,k)
!             csfeddysl(it,j,k)= csfeddy(it,j,k)
!          end do
!       end do
!    end do
!
!    do k=1,NK
!       do j=1,NJ
!          do itr=1,ntr
!             trbar= 0.d0
!             do i=1,NI
!                trbar= trbar + Tr(itr,i,j,k,n)
!             end do
!             trbar= trbar/dble(NI)
!             trbarsl(itr,j,k)= trbar
!          end do
!       end do
!    end do

   i= islice 
   do k=1,NK 
      zcave(k)= 0.d0 
      do j=1,NJ 
         zcave(k)= zcave(k)+ zc(i,j,k)*DL 
         zslice(j,k)= zc(i,j,k)*DL 
         yslice(j)= yc(j) 
         
!          s_ini_no_ghost(j,k) = s_ini(j,k)
!          T_ini_no_ghost(j,k) = T_ini(j,k)
!          s_restore_no_ghost(j,k) = s_restore(j,k)
!          T_restore_no_ghost(j,k) = T_restore(j,k)

         sslice(j,k)= s(i,j,k,n)        
         Tslice(j,k)= T(i,j,k,n)
         rslice(j,k)= rho(i,j,k) 
         uslice(j,k)= u(i,j,k,n) 
         vslice(j,k)= v(i,j,k,n) 
!          ugslice(j,k)=ugeo(i,j,k)
!          vgslice(j,k)=vgeo(i,j,k)
!          shearslice(j,k)=shear(i,j,k)
         wslice(j,k)= w(i,j,k,n) 
!          vorslice(j,k)= vor(i,j,k)
!          pvtslice(j,k)= pvt(i,j,k)
!          pv3slice(j,k)= pv3(i,j,k)
!          pslice  (j,k)= rp(i,j,k)
!          KzMomsl (j,k)= KzMom(i,j,k)
!          KzTrsl  (j,k)= KzTr(i,j,k)
!          do itr=1,ntr
!             Trslice(itr,j,k)= Tr(itr,i,j,k,n)
!          end do
!          do itr=1,nconsume
!             conslice(itr,j,k)= consump(i,j,k,itr)
!          end do
!
!          if (k.eq.1) then
!             dz = (zc(i,j,k+1)-zc(i,j,k))*DL
!             drhodz= (rho(i,j,k+1)-rho(i,j,k))/dz
!             drhobardz= (rhobar(j,k+1)-rhobar(j,k))/dz
!          else if (k.eq.NK) then
!             dz = (zc(i,j,k)-zc(i,j,k-1))*DL
!             drhodz= (rho(i,j,k)-rho(i,j,k-1))/dz
!             drhobardz= (rhobar(j,k)-rhobar(j,k-1))/dz
!          else
!             dz = (zc(i,j,k+1)-zc(i,j,k-1))*DL
!             drhodz= (rho(i,j,k+1)-rho(i,j,k-1))/dz
!             drhobardz= (rhobar(j,k+1)-rhobar(j,k-1))/dz
!          end if
!          n2slice(j,k) = -drhodz*9.81/R0
!          n2barslice(j,k) = -drhobardz*9.81/R0
!          n2slice(j,k) = freqN2(i,j,k)
!          n2barslice(j,k) = n2bar(j,k)
      end do 
      zcave(k)= zcave(k)/dble(NJ) 
   end do 
      
!    do j=1,NJ
!       stressxsp(j)= stress_top_x(NI/2,j)
!       stressysp(j)= stress_top_y(NI/2,j)
!       swrsl(j)= swr(j)
!       qlosssl(j)= qloss(j)
!    end do
                                                      
! =================================================================================================
!                                         Start writing
! =================================================================================================               
   WRITE(outname,'("xslice_",I3.3,".cdf")') islice

   if (counter_2d.eq.1) then 
      idDatFile =  nccre(TRIM(dirout)//outname,NCCLOB,rcode) 

      !---------------------------------------
      !              Dimension
      !---------------------------------------
      dims(1) = ncddef(idDatFile,'x',1,rcode) 
      dims(2) = ncddef(idDatFile,'y',NJ,rcode) 
      dims(3) = ncddef(idDatFile,'z',NK,rcode) 
      dims(4) = ncddef(idDatFile,'time',NCUNLIM,rcode) 

!       dims1d= dims(2)
!
!       dimswind(1)= dims(1)
!       dimswind(2)= dims(2)
!       dimswind(3)= dims(4)
!
!       dims4(1) = ncddef(idDatFile,'ntr',ntr,rcode)
!       dims4(2) = dims(2)
!       dims4(3) = dims(3)
!       dims4(4) = dims(4)
!
!       dims4c(1) = ncddef(idDatFile,'nconsume',nconsume,rcode)
!       dims4c(2) = dims(2)
!       dims4c(3) = dims(3)
!       dims4c(4) = dims(4)
!
!       dimstim(1) = ncddef(idDatFile,'x',1,rcode)
!       dimstim(2) = ncddef(idDatFile,'y',1,rcode)
!       dimstim(3) = ncddef(idDatFile,'z',1,rcode)
!       dimstim = ncddef(idDatFile,'time',NCUNLIM,rcode)
      dimstim = dims(4)  ! dimension for time day

      !---------------------------------------
      !                  ID
      !---------------------------------------
      idvy      = ncvdef(idDatFile,'yc',       NCFLOAT,1,dims(2), rcode) 
!       idvzave   = ncvdef(idDatFile, 'zcave',   NCFLOAT,1,dims1d,  rcode)
      idvz      = ncvdef(idDatFile,'zc',       NCFLOAT,3,dims,    rcode) 
      idvtim    = ncvdef(idDatFile,'day',      NCFLOAT,1,dimstim, rcode) 
         
!      idvstressx = ncvdef(idDatFile, 'stressx', NCFLOAT,3,dimswind,rcode)
!      idvstressy = ncvdef(idDatFile, 'stressy', NCFLOAT,3,dimswind,rcode)
!      idvswr     = ncvdef(idDatFile, 'swr',     NCFLOAT,3,dimswind,rcode) ! short wave radiation
!      idvqloss   = ncvdef(idDatFile, 'qloss',   NCFLOAT,3,dimswind,rcode) ! heat loss
!      idvcon     = ncvdef(idDatFile, 'consump', NCFLOAT,4,dims4c,  rcode) ! consumption nutrition
!      idvtr      = ncvdef(idDatFile, 'tr',      NCFLOAT,4,dims4,   rcode) ! tracer
!      idvtbar    = ncvdef(idDatFile, 'trbar',   NCFLOAT,4,dims4,   rcode)

      idvs       = ncvdef(idDatFile,'s',       NCFLOAT,4,dims,    rcode) ! salinity
      idvt       = ncvdef(idDatFile,'temp',    NCFLOAT,4,dims,    rcode)      
!       idvsini    = ncvdef(idDatFile,'sini',    NCFLOAT,4,dims,    rcode) !!! ADDED
!       idvtini    = ncvdef(idDatFile,'tini',    NCFLOAT,4,dims,    rcode) !!! ADDED
!       idvsrestore= ncvdef(idDatFile,'srestore',NCFLOAT,4,dims,    rcode) !!! ADDED
!       idvtrestore= ncvdef(idDatFile,'trestore',NCFLOAT,4,dims,    rcode) !!! ADDED

      idvrho     = ncvdef(idDatFile,'rho',     NCFLOAT,4,dims,    rcode) 
      idvu       = ncvdef(idDatFile,'u' ,      NCFLOAT,4,dims,    rcode) 
!       idvug      = ncvdef(idDatFile,'ug',      NCFLOAT,4,dims,    rcode)
      idvv       = ncvdef(idDatFile,'v' ,      NCFLOAT,4,dims,    rcode) 
!       idvvg      = ncvdef(idDatFile,'vg',      NCFLOAT,4,dims,    rcode) 
      idvw       = ncvdef(idDatFile,'w' ,      NCFLOAT,4,dims,    rcode) 
!       idvvor     = ncvdef(idDatFile,'vor',     NCFLOAT,4,dims,    rcode)
!       idvshear   = ncvdef(idDatFile,'shear',   NCFLOAT,4,dims,    rcode)
!       idvn2      = ncvdef(idDatFile,'n2',      NCFLOAT,4,dims,    rcode)
!       idvKzMom   = ncvdef(idDatFile,'KzMom',   NCFLOAT,4,dims,    rcode) ! momentum tracer
!       idvKzTr    = ncvdef(idDatFile,'KzTr',    NCFLOAT,4,dims,    rcode)
!       idvfb      = ncvdef(idDatFile,'fricb',   NCFLOAT,4,dims,    rcode)
!       idvpvt     = ncvdef(idDatFile,'pvt',     NCFLOAT,4,dims,    rcode) ! potential vorticity
!       idvpv3     = ncvdef(idDatFile,'pv3',     NCFLOAT,4,dims,    rcode) ! pv in the vertical
!       idvp       = ncvdef(idDatFile,'pre',     NCFLOAT,4,dims,    rcode)
!
!       idvrhbar   = ncvdef(idDatFile,'rhobar',  NCFLOAT,4,dims,    rcode)
!       idvn2bar   = ncvdef(idDatFile,'n2bar',   NCFLOAT,4,dims,    rcode)
!       idvpsiv    = ncvdef(idDatFile,'psivbar', NCFLOAT,4,dims,    rcode)
!       idvpsiw    = ncvdef(idDatFile,'psiwbar', NCFLOAT,4,dims,    rcode)
!       idvbysection = ncvdef(idDatFile,'by',    NCFLOAT,4,dims,    rcode) ! horizontal buoyancy
!       idvby      = ncvdef(idDatFile,'bybar',   NCFLOAT,4,dims,    rcode)
!       idvbz      = ncvdef(idDatFile,'bzbar',   NCFLOAT,4,dims,    rcode)
!       idvvb      = ncvdef(idDatFile,'vbbar',   NCFLOAT,4,dims,    rcode)
!       idvwb      = ncvdef(idDatFile,'wbbar',   NCFLOAT,4,dims,    rcode)
!       idvpe      = ncvdef(idDatFile,'psieddy', NCFLOAT,4,dims,    rcode)
!       idvpv      = ncvdef(idDatFile,'pvbar',   NCFLOAT,4,dims,    rcode)
!       idvwpv     = ncvdef(idDatFile,'wpvbar',  NCFLOAT,4,dims,    rcode)
!
!       idvdivfreyn=ncvdef(idDatFile,'divfreyn', NCFLOAT,4,dims,    rcode)
!       idvdivfeddy=ncvdef(idDatFile,'divfeddy', NCFLOAT,4,dims,    rcode)
!
!       idvcy      = ncvdef(idDatFile,'cybar',   NCFLOAT,4,dims,    rcode)
!       idvcz      = ncvdef(idDatFile,'czbar',   NCFLOAT,4,dims,    rcode)
!       idvvc      = ncvdef(idDatFile,'vcbar',   NCFLOAT,4,dims,    rcode)
!       idvwc      = ncvdef(idDatFile,'wcbar',   NCFLOAT,4,dims,    rcode)
!       idvpec     = ncvdef(idDatFile,'psieddytracer',NCFLOAT,4,dims,rcode)                                

      CALL ncendf(idDatFile,rcode) 

   else 
      !---------------------------------------
      !                  ID
      !---------------------------------------
      idDatFile = ncopn(TRIM(dirout)//outname, NCWRITE,rcode) 
!       call ncredf(idDatFile,rcode)
      idvtim     = NCVID(idDatFile,'day',     RCODE) 
!       idvstressx = NCVID(idDatFile,'stressx', RCODE)
!       idvstressy = NCVID(idDatFile,'stressy', RCODE)
!       idvswr     = NCVID(idDatFile,'swr',     RCODE)
!       idvqloss   = NCVID(idDatFile,'qloss',   RCODE)
!       idvcon     = NCVID(idDatFile,'consump', RCODE)
!       idvtr      = NCVID(idDatFile,'tr',      RCODE)
!       idvtbar    = NCVID(idDatFile,'trbar',   RCODE)

      idvs       = NCVID(idDatFile,'s',       RCODE)         
      idvt       = NCVID(idDatFile,'temp',    RCODE)
!       idvsini    = NCVID(idDatFile,'sini',    RCODE) !!! ADDED
!       idvtini    = NCVID(idDatFile,'tini',    RCODE) !!! ADDED
!       idvsrestore= NCVID(idDatFile,'srestore',RCODE) !!! ADDED
!       idvtrestore= NCVID(idDatFile,'trestore',RCODE) !!! ADDED

      idvrho = NCVID(idDatFile, 'rho', RCODE) 
      idvu   = NCVID(idDatFile, 'u' , RCODE) 
!       idvug  = NCVID(idDatFile, 'ug', RCODE) 
      idvv   = NCVID(idDatFile, 'v' , RCODE) 
!       idvvg  = NCVID(idDatFile, 'vg', RCODE) 
      idvw   = NCVID(idDatFile, 'w' , RCODE) 
!       idvvor = NCVID(idDatFile, 'vor', RCODE)
!       idvshear = NCVID(idDatFile, 'shear', RCODE)
!       idvn2  = NCVID(idDatFile, 'n2', RCODE)
!       idvKzMom  = NCVID(idDatFile, 'KzMom', RCODE)
!       idvKzTr  = NCVID(idDatFile, 'KzTr', RCODE)
!       idvfb = NCVID(idDatFile, 'fricb', RCODE)
!       idvpvt = NCVID(idDatFile, 'pvt', RCODE)
!       idvpv3 = NCVID(idDatFile, 'pv3', RCODE)
!       idvp   = NCVID(idDatFile, 'pre', RCODE)
!
!       idvrhbar = NCVID(idDatFile, 'rhobar', RCODE)
!       idvn2bar = NCVID(idDatFile, 'n2bar', RCODE)
!       idvpsiv = NCVID(idDatFile, 'psivbar', RCODE)
!       idvpsiw = NCVID(idDatFile, 'psiwbar', RCODE)
!       idvbysection = NCVID(idDatFile, 'by', RCODE)
!       idvby = NCVID(idDatFile, 'bybar', RCODE)
!       idvbz = NCVID(idDatFile, 'bzbar', RCODE)
!       idvvb = NCVID(idDatFile, 'vbbar', RCODE)
!       idvwb = NCVID(idDatFile, 'wbbar', RCODE)
!       idvpe = NCVID(idDatFile, 'psieddy', RCODE)
!       idvpv = NCVID(idDatFile, 'pvbar', RCODE)
!       idvwpv = NCVID(idDatFile, 'wpvbar', RCODE)
!
!       idvdivfreyn = ncvid(idDatFile,'divfreyn',rcode)
!       idvdivfeddy = ncvid(idDatFile,'divfeddy',rcode)
!
!       idvcy = NCVID(idDatFile,'cybar',rcode)
!       idvcz = NCVID(idDatFile,'czbar',rcode)
!       idvvc = NCVID(idDatFile,'vcbar',rcode)
!       idvwc = NCVID(idDatFile,'wcbar',rcode)
!       idvpec = NCVID(idDatFile,'psieddytracer',rcode)

   endif 

! =================================================================================================
   !---------------------------------------
   !           Count & start
   !---------------------------------------
   counttim =1 
   count1d= NK 

   count(1)= 1
   count(2)= NJ 
   count(3)= NK 
   count(4)= 1

   countwind(1)= 1
   countwind(2)= NJ
   countwind(3)= 1

   count4(1)= ntr 
   count4(2)= NJ 
   count4(3)= NK 
   count4(4)= 1 

   count4c(1)= nconsume
   count4c(2)= NJ 
   count4c(3)= NK 
   count4c(4)= 1 

   !---------------------------------------
   start2d(1)= 1 
   start2d(2)= 1 

   start(1)= 1 
   start(2)= 1 
   start(3)= 1
   start(4)= counter_2d 

   startwind(1)= 1
   startwind(2)= 1
   startwind(3)= counter_2d

   start4(1)= 1 
   start4(2)= 1 
   start4(3)= 1 
   start4(4)= counter_2d

   !---------------------------------------
   !                Write
   !---------------------------------------
   if (counter_2d.eq.1) then 
      CALL ncvpt(idDatFile,idvy,     start(2),  count(2), yslice,     rcode) 
!       CALL ncvpt(idDatFile,idvzave,  start(1),  count1d,  zcave,      rcode)
      CALL ncvpt(idDatFile,idvz,     start,     count,    zslice,     rcode) 
   endif 

   CALL ncvpt(idDatFile,idvtim,      start(4),  counttim, time_days,  rcode) 
!    CALL ncvpt(idDatFile,idvstressx,startwind, countwind,stressxsp,    rcode)
!    CALL ncvpt(idDatFile,idvstressy,startwind, countwind,stressysp,    rcode)
!    CALL ncvpt(idDatFile,idvswr,    startwind, countwind,swrsl,        rcode)
!    CALL ncvpt(idDatFile,idvqloss,  startwind, countwind,qlosssl,      rcode)
!    CALL ncvpt(idDatFile,idvcon,    start4,    count4c,  conslice,     rcode)
!    CALL ncvpt(idDatFile,idvtr,     start4,    count4,   Trslice,      rcode)
!    CALL ncvpt(idDatFile,idvtbar,   start4,    count4,   Trbarsl,      rcode)

   CALL ncvpt(idDatFile,idvs,        start, count, sslice,            rcode) 
   CALL ncvpt(idDatFile,idvt,        start, count, Tslice,            rcode) 
!    CALL ncvpt(idDatFile,idvsini,     start, count, s_ini_no_ghost,    rcode)  !!! ADDED
!    CALL ncvpt(idDatFile,idvtini,     start, count, T_ini_no_ghost,    rcode)  !!! ADDED
!    CALL ncvpt(idDatFile,idvsrestore, start, count, s_restore_no_ghost,rcode)  !!! ADDED
!    CALL ncvpt(idDatFile,idvtrestore, start, count, T_restore_no_ghost,rcode)  !!! ADDED

   CALL ncvpt(idDatFile,idvrho,      start, count, rslice,            rcode) 
   CALL ncvpt(idDatFile,idvu ,       start, count, uslice ,           rcode) 
!    CALL ncvpt(idDatFile,idvug,       start, count, ugslice,           rcode)
   CALL ncvpt(idDatFile,idvv ,       start, count, vslice ,           rcode) 
!    CALL ncvpt(idDatFile,idvvg,       start, count, vgslice,           rcode) 
   CALL ncvpt(idDatFile,idvw ,       start, count, wslice ,           rcode) 
!    CALL ncvpt(idDatFile,idvvor,      start, count, vorslice,          rcode)
!    CALL ncvpt(idDatFile,idvshear,    start, count, shearslice,        rcode)
!    CALL ncvpt(idDatFile,idvn2,       start, count, n2slice,           rcode)
!    CALL ncvpt(idDatFile,idvKzMom,    start, count, KzMomsl,           rcode)
!    CALL ncvpt(idDatFile,idvKzTr,     start, count, KzTrsl,            rcode)
!    CALL ncvpt(idDatFile,idvfb,       start, count, fricbsl,           rcode)
!    CALL ncvpt(idDatFile,idvpvt,      start, count, pvtslice,          rcode)
!    CALL ncvpt(idDatFile,idvpv3,      start, count, pv3slice,          rcode)
!    CALL ncvpt(idDatFile,idvp,        start, count, pslice,            rcode)
!
!    CALL ncvpt(idDatFile,idvrhbar,    start, count, rhobarslice,       rcode)
!    CALL ncvpt(idDatFile,idvn2bar,    start, count, n2barslice,        rcode)
!    CALL ncvpt(idDatFile,idvpsiv,     start, count, psivslice,         rcode)
!    CALL ncvpt(idDatFile,idvpsiw,     start, count, psiwslice,         rcode)
!    CALL ncvpt(idDatFile,idvbysection,start, count, bysection,         rcode)
!    CALL ncvpt(idDatFile,idvby,       start, count, byslice,           rcode)
!    CALL ncvpt(idDatFile,idvbz,       start, count, bzslice,           rcode)
!    CALL ncvpt(idDatFile,idvvb,       start, count, vbbarsl,           rcode)
!    CALL ncvpt(idDatFile,idvwb,       start, count, wbbarsl,           rcode)
!    CALL ncvpt(idDatFile,idvpe,       start, count, psieddysl,         rcode)
!    CALL ncvpt(idDatFile,idvpv,       start, count, pvbarsl,           rcode)
!    CALL ncvpt(idDatFile,idvwpv,      start, count, wpvbarsl,          rcode)
!
!    CALL ncvpt(idDatFile,idvdivfreyn, start, count, Freyndivsl,        rcode)
!    CALL ncvpt(idDatFile,idvdivfeddy, start, count, Feddydivsl,        rcode)
!
!    CALL ncvpt(idDatFile,idvcy,       start4,count4,cybarsl,           rcode)
!    CALL ncvpt(idDatFile,idvcz,       start4,count4,czbarsl,           rcode)
!    CALL ncvpt(idDatFile,idvvc,       start4,count4,vcbarsl,           rcode)
!    CALL ncvpt(idDatFile,idvwc,       start4,count4,wcbarsl,           rcode)
!    CALL ncvpt(idDatFile,idvpec,      start4,count4,csfeddysl,         rcode)

   CALL ncclos(idDatFile,rcode) 

   return 
END                                           
