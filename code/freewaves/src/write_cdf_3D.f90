subroutine write_cdf_3D(stepl,n)

   USE header,only: NI,NJ,NK,xc,yc,zc,p,h,consump,T,s,rho,Tr,u,v,w,vor,pv3,uf,vf,wf,Kz,conv,con100,nconsume,dirout,rc_kind,ntr,&
                    zf,Jac,UL,WL,LEN,DL,Kz,Kx,Ky,iv_compute_kh,iv_compute_kz  ! drpx,drpy
! #ifdef extra_output
   USE header,only: term_advcT, term_advcU, term_advcV, term_advcW, term_restT, term_restU,& !!! ADDED
       term_dissTh,term_dissUh,term_dissVh,term_dissWh,term_dissTv,term_dissUv,term_dissVv,term_dissWv!!! ADDED
! #endif

#include "netcdf.inc"                                                                        

   integer i,j,k,n,stepl,nstp,it 
   character(LEN=150) outname_data,outname_face
   REAL(kind=rc_kind) :: rcode
!    REAL(kind=rc_kind) :: Trwrite(0:NI+1,0:NJ+1,0:NK+1,ntr)
   REAL(kind=rc_kind) :: z(0:NI+1,0:NJ+1,0:NK+1) 
   REAL(kind=rc_kind) :: smax,Tmax,umax,vmax,wmax,smin,Tmin,umin,vmin,wmin 

   integer start(3),  start2d(2),start4(4),start2dtr(3)
   integer count(3),  count2d(2),count4(4),count2dtr(3),count4consump(4)
   integer countuf(3),countvf(3),countwf(3)                                                   !! ADDED 
   integer dims(3),   dims2d(2), dims4(4), dims2dtr(3), dimsconsump(4)
   integer dimuf(3),  dimvf(3),dimwf(3)

   integer :: iddatfile, idigit, idudx, idudy, idudz, idvbysection, idvby, idvbz, idvcon, idvcy, idvbx
   integer :: idvcz, idvc, idvdivfeddy, idvdivfreyn, idvdx, idvdz, idvd, idvfb, idvh, idvn2bar
   integer :: idvn2, idvnsq100m, idvnsq30m, idvpe,idvpsiv,idvpsiw,idvpv, idvp,idvrhbar,idvrho,idvrnk
   integer :: idvstrain,idvstress,idvstr,idvs,idvtbar,idvtemp,idvtim,idvtr,idvt,idvu,idvvb,idvvc,idvvor,idvv
   integer :: idvwb,idvwc,idvwpv,idvw,idvy,idvzsave,idvz,idvz3,idwdx,idwdy,idwdz,iimday,ipos
   integer :: idvx,idvcon100,idFaceFile,idvuf,idvvf,idvwf,idvKf
   
   integer :: idvzf,idvJac,idvptfx,idvptfy,idvptfz!,idvptc                                 !! ADDED
   integer :: idvKz,idvKx,idvKy                                                               !! ADDED
!    integer :: idvdrpy, idvdrpx,idvugeo,idvvgeo,idvdum1,idvdum2,idvdum3,idvfdiv              !! ADDED
! #ifdef extra_output
   integer :: idvadvcT, idvadvcU, idvadvcV, idvadvcW, idvrestT, idvrestU                      !! ADDED
   integer :: idvdissTh,idvdissUh,idvdissVh,idvdissWh,idvdissTv,idvdissUv,idvdissVv,idvdissWv !! ADDED
! #endif

!    REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1) :: ugeo, vgeo                     !! ADDED
   REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1) :: ptc                            !! ADDED
   REAL(kind=rc_kind), dimension(    0:NI,    NJ  ,   NK  ) :: ptfx                           !! ADDED
   REAL(kind=rc_kind), dimension(      NI,  0:NJ,     NK  ) :: ptfy                           !! ADDED
   REAL(kind=rc_kind), dimension(      NI,    NJ,   0:NK  ) :: ptfz                           !! ADDED

! =================================================================================================
!                                       Initialization
! =================================================================================================
   DATA start /1, 1, 1/ 
!    DATA start4 /1, 1, 1, 1/
   DATA start2d /1, 1/
!    DATA start2dtr /1, 1, 1/
! 
!    do k=0,NK+1
!    do j=0,NJ+1
!    do i=0,NI+1
!       z(i,j,k)= zc(i,j,k)
!       do it=1,ntr
!          Trwrite(i,j,k,it)= Tr(it,i,j,k,n)
!       end do
!    end do
!    end do
!    end do
! 
!---------------------------------------
!     For the sake of better plots 
!---------------------------------------
!   do j=0,NJ+1
!     do i=0,NI+1
!       z(i,j,NK+1)= 0.d0 ;   z(i,j,0)= 0.5*(z(i,j,0)+z(i,j,1));
!       s(i,j,NK+1,n)= s(i,j,NK,n); s(i,j,0,n)= s(i,j,1,n);
!       T(i,j,NK+1,n)= T(i,j,NK,n); T(i,j,0,n)= T(i,j,1,n);
!       rho(i,j,NK+1)= rho(i,j,NK); rho(i,j,0)= rho(i,j,1);
!       u(i,j,NK+1,n)= u(i,j,NK,n); u(i,j,0,n)= u(i,j,1,n);
!       v(i,j,NK+1,n)= v(i,j,NK,n); v(i,j,0,n)= v(i,j,1,n);
!       w(i,j,NK+1,n)= w(i,j,NK,n); w(i,j,0,n)= w(i,j,1,n);
!       vor(i,j,NK+1)= vor(i,j,NK); vor(i,j,0)= vor(i,j,1);
!       pv(i,j,NK+1) = pv(i,j,NK) ; pv(i,j,0) = pv(i,j,1) ;
!       do it=1,ntr
!         Trwrite(i,j,NK+1,it)= Trwrite(i,j,NK,it); Trwrite(i,j,0,it)= Trwrite(i,j,1,it);
!       end do
!     end do
!   end do

   call diag_pv(n); ! pvt=pv1+pv2+pv3;
   call diag_ptc(ptc)                 !! ADDED
   call intpol_pt(ptc,ptfx,ptfy,ptfz) !! ADDED
!    call gradup3df                     !! ADDED
!    call uv_geostrophic(ugeo,vgeo)

!---------------------------------------
!          Output file names
!---------------------------------------                                             
   WRITE(outname_data,'("full_",I5.5,".cdf")') stepl ! Cell centers 
   WRITE(outname_face,'("face_",I5.5,".cdf")') stepl ! Cell faces
!    Overwrite the face file, so only the last one is saved
!    WRITE(outname_face,'("face_values.cdf")')       ! Cell faces

! =================================================================================================
!                                        Center values
! =================================================================================================
   idDatFile =  nccre(TRIM(dirout)//outname_data,NCCLOB,rcode)
  
   !---------------------------------------
   !              Dimension
   !---------------------------------------
   ! 3D dimensions for most variables
   count(1)= NI+2; count(2)= NJ+2; count(3)= NK+2;
   dims(1) = ncddef(idDatFile,'x',    NI+2,rcode);  
   dims(2) = ncddef(idDatFile,'y',    NJ+2,rcode); 
   dims(3) = ncddef(idDatFile,'sigma',NK+2,rcode);

   ! 2D dimensions (for h) 
   count2d(1)= NJ+2 ; count2d(2)= NK+2; 
   dims2d(1)= dims(2); dims2d(2)= dims(3);

!   ! 2D dimensions for tracer (currently not used)
!   count2dtr(1)= NI+2; count2dtr(2)= NJ+2; count2dtr(3)= ntr;
!   dims2dtr(1)= dims(1); dims2dtr(2)= dims(2); dims2dtr(3)= dims4(4);
! 
!   ! 3D dimensions (for tracer)                                          
!   count4(1)= NI+2; count4(2)= NJ+2; count4(3)= NK+2; count4(4)= ntr;
!   dims4(1) = dims(1); dims4(2) = dims(2); dims4(3) = dims(3); dims4(4) = ncddef(idDatFile,'ntr',ntr,rcode);
! 
!   ! 3D dimensions (for consump)
!   count4consump(1)= NI+2; count4consump(2)= NJ+2; count4consump(3)= NK+2; count4consump(4)= nconsume;
!   dimsconsump(1) = dims(1); dimsconsump(2) = dims(2);
!   dimsconsump(3) = dims(3); dimsconsump(4) = ncddef(idDatFile,'ntrcon',nconsume,rcode);
! 
!   ! 3D dimensions for divergence of energy flux
!   dimfdiv(1) = ncddef(idDatFile,'fdiv',1, rcode);   !! ADDED

   !---------------------------------------
   !                  ID
   !---------------------------------------
!    if (stepl.lt.0.5) then
!    idvx     =ncvdef(idDatFile,'xc',     NCDOUBLE,1,dims(1),rcode)
!    idvy     =ncvdef(idDatFile,'yc',     NCDOUBLE,1,dims(2),rcode)
!    idvz     =ncvdef(idDatFile,'zc',     NCDOUBLE,3,dims,   rcode)
!    idvzf    =ncvdef(idDatFile,'zf',     NCDOUBLE,3,dims,   rcode)
!    endif
!    idvh     =ncvdef(idDatFile,'h',      NCDOUBLE,2,dims2d, rcode)
!    idvc     =ncvdef(idDatFile,'consump',NCDOUBLE,4,dimsconsump,rcode)
!    idvtr    =ncvdef(idDatFile,'tr',     NCDOUBLE,4,dims4,  rcode)
!    idvs     =ncvdef(idDatFile,'s',      NCDOUBLE,3,dims,   rcode)
!    idvt     =ncvdef(idDatFile,'temp',   NCDOUBLE,3,dims,   rcode)
   idvrho   =ncvdef(idDatFile,'rho',    NCDOUBLE,3,dims,   rcode)
!    idvp     =ncvdef(idDatFile,'p',      NCDOUBLE,3,dims,   rcode)
   idvu     =ncvdef(idDatFile,'u',      NCDOUBLE,3,dims,   rcode)
   idvv     =ncvdef(idDatFile,'v',      NCDOUBLE,3,dims,   rcode)
!    idvugeo  =ncvdef(idDatFile,'ugeo',   NCDOUBLE,3,dims,   rcode)
!    idvvgeo  =ncvdef(idDatFile,'vgeo',   NCDOUBLE,3,dims,   rcode)
   idvw     =ncvdef(idDatFile,'w',      NCDOUBLE,3,dims,   rcode)
!    idvz3    =ncvdef(idDatFile,'vor',    NCDOUBLE,3,dims,   rcode)
   idvpv    =ncvdef(idDatFile,'pv',     NCDOUBLE,3,dims,   rcode)
!    idvcon   =ncvdef(idDatFile,'conv',   NCshort, 3,dims,   rcode)
!    idvcon100=ncvdef(idDatFile,'con100', NCshort, 3,dims,   rcode)
!    idvJac   =ncvdef(idDatFile,'Jac',    NCDOUBLE,3,dims,   rcode) !! ADDED
!    idvptc   =ncvdef(idDatFile,'ptc',    NCDOUBLE,3,dims,   rcode) !! ADDED
!    idvdum1  =ncvdef(idDatFile,'dum1',   NCDOUBLE,3,dims,   rcode) !! ADDED
!    idvdum2  =ncvdef(idDatFile,'dum2',   NCDOUBLE,3,dims3d, rcode) !! ADDED
!    idvdum3  =ncvdef(idDatFile,'dum3',   NCDOUBLE,3,dims3d, rcode) !! ADDED
!    idvfdiv  =ncvdef(idDatFile,'fdiv',   NCDOUBLE,1,dimfdiv,rcode) !! ADDED
!    idvdrpy  =ncvdef(idDatFile,'drpy',   NCDOUBLE,3,dims,   rcode) !! ADDED
!    idvdrpx  =ncvdef(idDatFile,'drpx',   NCDOUBLE,3,dims,   rcode) !! ADDED
if (iv_compute_kh==1) then
   idvKx    =ncvdef(idDatFile,'Kx',     NCDOUBLE,3,dims,   rcode) !! ADDED
   idvKy    =ncvdef(idDatFile,'Ky',     NCDOUBLE,3,dims,   rcode) !! ADDED
endif
if (iv_compute_kz==1) then
   idvKz    =ncvdef(idDatFile,'Kz',     NCDOUBLE,3,dims,   rcode) !! ADDED
endif
! #ifdef extra_output
   idvadvcT =ncvdef(idDatFile,'advcT',  NCDOUBLE,3,dims,   rcode) !! ADDED
   idvadvcU =ncvdef(idDatFile,'advcU',  NCDOUBLE,3,dims,   rcode) !! ADDED
   idvadvcV =ncvdef(idDatFile,'advcV',  NCDOUBLE,3,dims,   rcode) !! ADDED
   idvadvcW =ncvdef(idDatFile,'advcW',  NCDOUBLE,3,dims,   rcode) !! ADDED
   idvrestT =ncvdef(idDatFile,'restT',  NCDOUBLE,2,dims2d, rcode) !! ADDED
   idvrestU =ncvdef(idDatFile,'restU',  NCDOUBLE,2,dims2d, rcode) !! ADDED
   idvdissTh=ncvdef(idDatFile,'dissTh', NCDOUBLE,3,dims,   rcode) !! ADDED
   idvdissUh=ncvdef(idDatFile,'dissUh', NCDOUBLE,3,dims,   rcode) !! ADDED
   idvdissVh=ncvdef(idDatFile,'dissVh', NCDOUBLE,3,dims,   rcode) !! ADDED
   idvdissWh=ncvdef(idDatFile,'dissWh', NCDOUBLE,3,dims,   rcode) !! ADDED
   idvdissTv=ncvdef(idDatFile,'dissTv', NCDOUBLE,3,dims,   rcode) !! ADDED
   idvdissUv=ncvdef(idDatFile,'dissUv', NCDOUBLE,3,dims,   rcode) !! ADDED
   idvdissVv=ncvdef(idDatFile,'dissVv', NCDOUBLE,3,dims,   rcode) !! ADDED
   idvdissWv=ncvdef(idDatFile,'dissWv', NCDOUBLE,3,dims,   rcode) !! ADDED
! #endif
   !---------------------------------------
   !                Write
   !---------------------------------------
   ! 
!    if (stepl.lt.0.5) then
   CALL ncendf(idDatFile,rcode)
!    CALL ncvpt(idDatFile,idvx,     start(1),count(1),xc*DL,      rcode)
!    CALL ncvpt(idDatFile,idvy,     start(2),count(2),yc*DL,      rcode)
!    CALL ncvpt(idDatFile,idvz,     start,   count,   zc*DL,      rcode)
!    CALL ncvpt(idDatFile,idvzf,    start,   count,   zf(:,:,-1:NK)*DL,rcode)
!    endif
!    CALL ncvpt(idDatFile,idvh,     start2d, count2d, h,          rcode)
!    CALL ncvpt(idDatFile,idvc,     start4, count4consump,consump,rcode)
!    CALL ncvpt(idDatFile,idvtr,    start4,  count4,  Trwrite,    rcode)
!    CALL ncvpt(idDatFile,idvs,     start,   count,   s,          rcode)
!    CALL ncvpt(idDatFile,idvt,     start,   count,   T,          rcode)
   CALL ncvpt(idDatFile,idvrho,   start,   count,   rho,        rcode)
!    CALL ncvpt(idDatFile,idvp,     start,   count,   p,          rcode)
   CALL ncvpt(idDatFile,idvu,     start,   count,   u*UL,       rcode)
   CALL ncvpt(idDatFile,idvv,     start,   count,   v*UL,       rcode)
!    CALL ncvpt(idDatFile,idvugeo,  start,   count,   ugeo*UL,    rcode)
!    CALL ncvpt(idDatFile,idvvgeo,  start,   count,   vgeo*UL,    rcode)
   CALL ncvpt(idDatFile,idvw,     start,   count,   w*WL,       rcode)
!    CALL ncvpt(idDatFile,idvz3,    start,   count,   vor,        rcode)
   CALL ncvpt(idDatFile,idvpv,    start,   count,   pv3,        rcode)
!    CALL ncvpt(idDatFile,idvcon,   start,   count,   conv,       rcode)
!    CALL ncvpt(idDatFile,idvcon100,start,   count,   con100,     rcode)
!    CALL ncvpt(idDatFile,idvJac,   start,   count,   Jac,        rcode) !! ADDED
!    CALL ncvpt(idDatFile,idvptc,   start,   count,   ptc,        rcode) !! ADDED
!    CALL ncvpt(idDatFile,idvdum1,  start,   count,   dum1,       rcode) !! ADDED
!    CALL ncvpt(idDatFile,idvdum2,  start,   countdiv,dum2,       rcode) !! ADDED
!    CALL ncvpt(idDatFile,idvdum3,  start,   countdiv,dum3,       rcode) !! ADDED
!    CALL ncvpt(idDatFile,idvfdiv,  start,   1,       fdiv,       rcode) !! ADDED
!    CALL ncvpt(idDatFile,idvdrpy,  start,   count,   drpx,       rcode) !! ADDED
!    CALL ncvpt(idDatFile,idvdrpx,  start,   count,   drpy,       rcode) !! ADDED
if (iv_compute_kh==1) then
   CALL ncvpt(idDatFile,idvKx,    start,   count,   Kx,         rcode) !! ADDED
   CALL ncvpt(idDatFile,idvKy,    start,   count,   Ky,         rcode) !! ADDED
endif
if (iv_compute_kz==1) then
   CALL ncvpt(idDatFile,idvKz,    start,   count,   Kz,         rcode) !! ADDED
endif

! #ifdef extra_output
   CALL ncvpt(idDatFile,idvadvcT, start,   count,   term_advcT, rcode) !! ADDED
   CALL ncvpt(idDatFile,idvadvcU, start,   count,   term_advcU, rcode) !! ADDED
   CALL ncvpt(idDatFile,idvadvcV, start,   count,   term_advcV, rcode) !! ADDED
   CALL ncvpt(idDatFile,idvadvcW, start,   count,   term_advcW, rcode) !! ADDED
   CALL ncvpt(idDatFile,idvrestT, start2d, count2d, term_restT, rcode) !! ADDED
   CALL ncvpt(idDatFile,idvrestU, start2d, count2d, term_restU, rcode) !! ADDED
   CALL ncvpt(idDatFile,idvdissTh,start,   count,   term_dissTh,rcode) !! ADDED
   CALL ncvpt(idDatFile,idvdissUh,start,   count,   term_dissUh,rcode) !! ADDED
   CALL ncvpt(idDatFile,idvdissVh,start,   count,   term_dissVh,rcode) !! ADDED
   CALL ncvpt(idDatFile,idvdissWh,start,   count,   term_dissWh,rcode) !! ADDED
   CALL ncvpt(idDatFile,idvdissTv,start,   count,   term_dissTv,rcode) !! ADDED
   CALL ncvpt(idDatFile,idvdissUv,start,   count,   term_dissUv,rcode) !! ADDED
   CALL ncvpt(idDatFile,idvdissVv,start,   count,   term_dissVv,rcode) !! ADDED
   CALL ncvpt(idDatFile,idvdissWv,start,   count,   term_dissWv,rcode) !! ADDED
! #endif

CALL ncclos(idDatFile,rcode)

! =================================================================================================
!                                         Face values
! =================================================================================================
   idFaceFile =  nccre(TRIM(dirout)//outname_face,NCCLOB,rcode)

   !---------------------------------------
   !               Dimension
   !---------------------------------------
   ! 3D dimensions for face values
   countuf(1)= NI+1; countuf(2)= NJ  ; countuf(3)= NK  ;
   countvf(1)= NI  ; countvf(2)= NJ+1; countvf(3)= NK  ;
   countwf(1)= NI  ; countwf(2)= NJ  ; countwf(3)= NK+1;

!    dimuf(1) = ncddef(idFaceFile,'xi-ew',   NI+1,rcode)
!    dimuf(2) = ncddef(idFaceFile,'eta-ew',  NJ,  rcode)
!    dimuf(3) = ncddef(idFaceFile,'sigma-ew',NK,  rcode)

!    dimvf(1) = ncddef(idFaceFile,'xi-ns',   NI,  rcode)
!    dimvf(2) = ncddef(idFaceFile,'eta-ns',  NJ+1,rcode)
!    dimvf(3) = ncddef(idFaceFile,'sigma-ns',NK,  rcode)

   dimwf(1) = ncddef(idFaceFile,'xi-tb',   NI,  rcode)
   dimwf(2) = ncddef(idFaceFile,'eta-tb',  NJ,  rcode)
   dimwf(3) = ncddef(idFaceFile,'sigma-tb',NK+1,rcode)

   !---------------------------------------
   !                  ID
   !---------------------------------------
!    idvKf   = ncvdef(idFaceFile,'Kz',  NCDOUBLE,3,dimwf,rcode)
!    idvuf   = ncvdef(idFaceFile,'uf',  NCDOUBLE,3,dimuf,rcode)
!    idvvf   = ncvdef(idFaceFile,'vf',  NCDOUBLE,3,dimvf,rcode)
   idvwf   = ncvdef(idFaceFile,'wf',  NCDOUBLE,3,dimwf,rcode)
!    idvptfx = ncvdef(idFaceFile,'ptfx',NCDOUBLE,3,dimuf,rcode) !! ADDED
!    idvptfy = ncvdef(idFaceFile,'ptfy',NCDOUBLE,3,dimvf,rcode) !! ADDED
   idvptfz = ncvdef(idFaceFile,'ptfz',NCDOUBLE,3,dimwf,rcode)   !! ADDED

   !---------------------------------------
   !                Write
   !---------------------------------------
   CALL ncendf(idFaceFile,rcode)
!  CALL ncvpt(idFaceFile,idvKf,  start,countwf,Kz,           rcode)
!  CALL ncvpt(idFaceFile,idvuf,  start,countuf,uf*UL*LEN*DL, rcode)
!  CALL ncvpt(idFaceFile,idvvf,  start,countvf,vf*UL*LEN*DL, rcode)
   CALL ncvpt(idFaceFile,idvwf,  start,countwf,wf*UL*LEN*DL, rcode)
!  CALL ncvpt(idFaceFile,idvptfx,start,countuf,ptfx,         rcode) !! ADDED
!  CALL ncvpt(idFaceFile,idvptfy,start,countvf,ptfy,         rcode) !! ADDED
   CALL ncvpt(idFaceFile,idvptfz,start,countwf,ptfz,         rcode) !! ADDED
   CALL ncclos(idFaceFile,rcode)

! =================================================================================================

   return
end
