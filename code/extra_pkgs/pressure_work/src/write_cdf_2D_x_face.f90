subroutine write_cdf_2D_x_face(islice,counter_2d)

USE header, ONLY: NI, NJ, NK, zf, uf, vf, wf, dirout, rc_kind, DL, Jac, ux, vy, wz, EPS

#include "netcdf.inc"

character(LEN=150) outname

INTEGER :: i,j,k,islice,counter_2d
INTEGER :: idFaceFile, idvuf, idvvf, idvwf, idvptfx, idvptfy, idvptfz, idvzf
INTEGER :: start(3), countfx(3), countfy(3), countfz(3), dimfx(3), dimfy(3), dimfz(3), dim2z(2)

REAL(kind=rc_kind) :: rcode
REAL(kind=4      ) ::   zfslice(NJ,0:NK)
REAL(kind=4      ) ::   ufslice(NJ,NK),   vfslice(0:NJ,NK),   wfslice(NJ,0:NK)
REAL(kind=4      ) :: ptfxslice(NJ,NK), ptfyslice(0:NJ,NK), ptfzslice(NJ,0:NK)

REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1)  :: ptc
REAL(kind=rc_kind), dimension(    0:NI,    NJ  ,   NK  )  :: ptfx
REAL(kind=rc_kind), dimension(      NI,  0:NJ,     NK  )  :: ptfy
REAL(kind=rc_kind), dimension(      NI,    NJ,   0:NK  )  :: ptfz

!-------------------------------------------
! Write face values

  call diag_ptc (ptc)
  call intpol_pt(ptc,ptfx,ptfy,ptfz)

      i=islice

      do k=1,NK
      do j=1,NJ
        ufslice  (j,k)=uf  (i,j,k)
	ptfxslice(j,k)=ptfx(i,j,k)
      enddo
      enddo

      do k=1,NK
      do j=0,NJ
        vfslice  (j,k)=vf  (i,j,k)
	ptfyslice(j,k)=ptfy(i,j,k)
      enddo
      enddo

      do k=0,NK
      do j=1,NJ
        wfslice  (j,k)=wf  (i,j,k)
	ptfzslice(j,k)=ptfz(i,j,k)
!        zfslice  (j,k)=zf  (i,j,k)*DL
      enddo
      enddo
      
      WRITE(outname,'("xslice_face_",I3.3,".cdf")') islice

      if (counter_2d.eq.1) then 
      
         idFaceFile =  nccre(TRIM(dirout)//outname,NCCLOB,rcode) 

         dimfx(1)  = ncddef(idFaceFile,'xi-ew' ,NJ,rcode) 
         dimfx(2)  = ncddef(idFaceFile,'eta-ew',NK,rcode) 
         dimfx(3)  = ncddef(idFaceFile,'time'  ,NCUNLIM,rcode) 

         dimfy(1)  = ncddef(idFaceFile,'xi-ns' ,NJ+1,rcode) 
         dimfy(2)  = ncddef(idFaceFile,'eta-ns',NK  ,rcode) 
         dimfy(3)  = dimfx(3)

         dimfz(1)  = ncddef(idFaceFile,'xi-tb' ,NJ  ,rcode) 
         dimfz(2)  = ncddef(idFaceFile,'eta-tb',NK+1,rcode) 
	 dimfz(3)  = dimfx(3)

         dim2z(1) = dimfz(1)
         dim2z(2) = dimfz(2)
         
 !        idvzf   = ncvdef(idFaceFile,'zf',NCFLOAT,2,dim2z,rcode)

         idvuf   = ncvdef(idFaceFile,'uf',NCFLOAT,3,dimfx,rcode)
         idvvf   = ncvdef(idFaceFile,'vf',NCFLOAT,3,dimfy,rcode)
         idvwf   = ncvdef(idFaceFile,'wf',NCFLOAT,3,dimfz,rcode)

         idvptfx = ncvdef(idFaceFile,'ptfx',NCFLOAT,3,dimfx,rcode)
         idvptfy = ncvdef(idFaceFile,'ptfy',NCFLOAT,3,dimfy,rcode)
         idvptfz = ncvdef(idFaceFile,'ptfz',NCFLOAT,3,dimfz,rcode)

         CALL ncendf(idFaceFile,rcode) 
                                                                        
      else 

         idFaceFile = ncopn(TRIM(dirout)//outname, NCWRITE,rcode) 
         
	 idvuf      = NCVID(idFaceFile,'uf',RCODE) 
	 idvvf      = NCVID(idFaceFile,'vf',RCODE) 
	 idvwf      = NCVID(idFaceFile,'wf',RCODE) 

	 idvptfx    = NCVID(idFaceFile,'ptfx',RCODE) 
	 idvptfy    = NCVID(idFaceFile,'ptfy',RCODE) 
	 idvptfz    = NCVID(idFaceFile,'ptfz',RCODE) 

      endif

      start(1)=1 
      start(2)=1 
      start(3)=counter_2d 

      countfx(1)=NJ
      countfx(2)=NK
      countfx(3)=1

      countfy(1)=NJ+1
      countfy(2)=NK
      countfy(3)=1

      countfz(1)=NJ
      countfz(2)=NK+1
      countfz(3)=1

  !    if (counter_2d .eq. 1) then
  !      CALL ncvpt(idFaceFile,idvzf, start, countfz, zfslice, rcode) 
  !    endif

      CALL ncvpt(idFaceFile,idvuf, start, countfx, ufslice, rcode) 
      CALL ncvpt(idFaceFile,idvvf, start, countfy, vfslice, rcode) 
      CALL ncvpt(idFaceFile,idvwf, start, countfz, wfslice, rcode) 

      CALL ncvpt(idFaceFile,idvptfx, start, countfx, ptfxslice, rcode) 
      CALL ncvpt(idFaceFile,idvptfy, start, countfy, ptfyslice, rcode) 
      CALL ncvpt(idFaceFile,idvptfz, start, countfz, ptfzslice, rcode) 

      CALL ncclos(idFaceFile,rcode)

RETURN
END
