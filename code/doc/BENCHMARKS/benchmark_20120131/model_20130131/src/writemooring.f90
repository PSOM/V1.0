      subroutine writemooring(imooring,jmooring,counter_1d) 
!     ------------------------------------------------------------      
                                                                        
!      implicit none                                                    
  USE header, ONLY : NI, NJ, NK, ntr, u, v, w, vor, rho, HL, h, dx, dy, zc, DL, time_seconds, dirout, rc_kind

#include "netcdf.inc" 

      integer i,j,k,it,imooring,jmooring,counter_1d 
      REAL(kind=rc_kind) drhodz,dz,rcode
                                                                        
      REAL(kind=rc_kind) ugeo(0:NI+1,0:NJ+1,0:NK+1),vgeo(0:NI+1,0:NJ+1,0:NK+1)
      REAL(kind=4) zprof(NK),n2prof(NK),uprof(NK),vprof(NK),wprof(NK),     &
           ugeoprof(NK),vgeoprof(NK),vorprof(NK),rhoprof(NK)
                        

 integer :: iddatfile, idigit, idudx, idudy, idudz, idvbysection, idvby, idvbz, idvcon, idvcy, idvbx
  integer :: idvcz, idvc, idvdivfeddy, idvdivfreyn, idvdx, idvdz, idvd, idvfb, idvh, idvn2bar
  integer :: idvn2, idvnsq100m, idvnsq30m, idvpe,idvpsiv,idvpsiw,idvpv, idvp,idvrhbar,idvrho,idvrnk
  integer :: idvstrain,idvstress,idvstr,idvs,idvtbar,idvtemp,idvtim,idvtr,idvt,idvu,idvvb,idvvc,idvvor,idvv
  integer :: idvwb,idvwc,idvwpv,idvw,idvy,idvzsave,idvz,idwdx,idwdy,idwdz,ipos

  REAL(kind=rc_kind) ::  vol
  integer :: idvx,idvug,idvvg


                                                
      character(LEN=150) outname 
                                                                        
      integer start(2),count(2),dims(2) 

!      call uv_geostrophic(ugeo,vgeo)

      i=imooring                                                           
      j=jmooring                                                       
      do k=1,NK                                                        
         rhoprof(k)= rho(i,j,k)                                        
         zprof(k) = zc(i,j,k)*DL                                       
         uprof(k) = u(i,j,k,0)                                         
         vprof(k) = v(i,j,k,0)                                         
         wprof(k) = w(i,j,k,0)                                         
         ugeoprof(k) = ugeo(i,j,k)                                     
         vgeoprof(k) = vgeo(i,j,k)                                     
         vorprof(k) = vor(i,j,k) 
      end do                                                           
!      hxy= 0.5*HL*(0.5*(h(i+1,j+1) - h(i-1,j+1))/dx                    
!     &     -0.5*(h(i+1,j-1) - h(i-1,j-1))/dx )/dy                      
                                                                        
      do k=1,NK 
         if (k.eq.1) then 
            dz = (zprof(k+1)-zprof(k)) 
            drhodz= (rhoprof(k+1)-rhoprof(k))/dz 
         else if (k.eq.NK) then 
            dz = (zprof(k)-zprof(k-1)) 
            drhodz= (rhoprof(k)-rhoprof(k-1))/dz 
         else 
            dz = (zprof(k+1)-zprof(k-1)) 
            drhodz= (rhoprof(k+1)-rhoprof(k-1))/dz 
         end if 
         n2prof(k) = -drhodz*9.81/(rhoprof(k)+1000.) 
      end do 
                                                                        
!     write to a netcdf file                                            
                
     WRITE(outname,'("mooring_",I3.3,"_",I3.3,".cdf")') imooring,jmooring                                                        
                                                                        
      if (counter_1d.eq.1) then 
         idDatFile =  nccre(TRIM(dirout)//outname,NCCLOB,rcode) 
                                                                        
         dims(1) = ncddef(idDatFile,'z',NK,rcode) 
         dims(2) = ncddef(idDatFile,'time',NCUNLIM,rcode) 
                                                                        
         idvz = ncvdef(idDatFile,'zc',NCFLOAT,1,dims(1),rcode) 
!         idvhxy = ncvdef(idDatFile,'hxy',NCFLOAT,1,dims(2),rcode)     
         idvrho = ncvdef(idDatFile,'rho',NCFLOAT,2,dims,rcode) 
         idvn2 = ncvdef(idDatFile,'n2',NCFLOAT,2,dims,rcode) 
         idvu = ncvdef(idDatFile,'u',NCFLOAT,2,dims,rcode) 
         idvv= ncvdef(idDatFile,'v',NCFLOAT,2,dims,rcode) 
         idvug = ncvdef(idDatFile,'ugeo',NCFLOAT,2,dims,rcode) 
         idvvg= ncvdef(idDatFile,'vgeo',NCFLOAT,2,dims,rcode) 
         idvw= ncvdef(idDatFile,'w',NCFLOAT,2,dims,rcode) 
         idvvor= ncvdef(idDatFile,'vor',NCFLOAT,2,dims,rcode) 
                                                                        
         CALL ncendf(idDatFile,rcode) 
                                                                        
      else 
                                                                        
         idDatFile = ncopn(TRIM(dirout)//outname, NCWRITE,rcode) 
                                                                        
!         idvhxy = NCVID(idDatFile, 'hxy', RCODE)                       
         idvrho = NCVID(idDatFile, 'rho', RCODE) 
         idvn2 = NCVID(idDatFile, 'n2', RCODE) 
         idvu = NCVID(idDatFile,'u', RCODE) 
         idvv= NCVID(idDatFile,'v', RCODE) 
         idvug = NCVID(idDatFile,'ugeo', RCODE) 
         idvvg= NCVID(idDatFile,'vgeo', RCODE) 
         idvw= NCVID(idDatFile,'w', RCODE) 
         idvvor= NCVID(idDatFile,'vor', RCODE) 
                                                                        
      endif 
      count(1)= NK 
      count(2)= 1 
                                                                        
      start(1)= 1 
      start(2)= counter_1d
                                                                        
      if (counter_1d.eq.1) then 
         CALL ncvpt(idDatFile,idvz, start(1), count(1), zprof, rcode) 
      endif 
!      CALL ncvpt(idDatFile,idvhxy, start, count(2), hxy, rcode)        
      CALL ncvpt(idDatFile,idvrho, start, count, rhoprof, rcode) 
      CALL ncvpt(idDatFile,idvn2, start, count, n2prof, rcode) 
      CALL ncvpt(idDatFile,idvu, start, count, uprof, rcode) 
      CALL ncvpt(idDatFile,idvv, start, count, vprof, rcode) 
      CALL ncvpt(idDatFile,idvug, start, count, ugeoprof, rcode) 
      CALL ncvpt(idDatFile,idvvg, start, count, vgeoprof, rcode) 
      CALL ncvpt(idDatFile,idvw, start, count, wprof, rcode) 
      CALL ncvpt(idDatFile,idvvor, start, count, vorprof, rcode) 
                                                                        
      CALL ncclos(idDatFile,rcode) 
                                                                        
      return 
      END                                           


! subroutine uv_geostrophic(ugeo,vgeo)
! !     ---------------------------------------------                     
!   !   Finds the geostrophic velocities ugeo,vgeo
!   USE header
!   USE rpgrads
! !     modified for periodicew bc                                        
! !     ---------------------------                                       
! !     Sets up the initial velocities so that they are in geostrophic bal
! !     with the initial density field.                                   
! !      implicit logical (a-z)                                           
!       integer i,j,k,n,imax,jmax,kmax,m 
!       double precision uxi,vyj,hxi,heta,hx,hy,px,py,ujfc 
!       double precision res,resmax
!       double precision ainv,be2,fac2,wzsk,wxsk,wysk,Jack,pxi,peta,      &
!      &     psig,pgrad,con,pz                                            
!       real*8 ugeo(0:NI+1,0:NJ+1,0:NK+1),vgeo(0:NI+1,0:NJ+1,0:NK+1)
!                                                                         
!       call rpevalgrad(n) 
!                                                                         
!       kaph1= 1.d0 -kappah 
!       con = 1.0 -qpr 
!                                                                         
!       do j=1,NJ 
!          do i=1,NI 
!             hxi= 0.5d0*( h(i+1,j)-h(i-1,j) ) 
!             heta= 0.5d0*( h(i,j+1)-h(i,j-1) ) 
!             hx= ux(i,j)*hxi +vx(i,j)*heta 
!             hy= uy(i,j)*hxi +vy(i,j)*heta 
!             do k=1,NK 
!                pxi= 0.5d0*(p(i+1,j,k)-p(i-1,j,k)) 
!                peta= 0.5d0*(p(i,j+1,k)-p(i,j-1,k)) 
!                psig= 0.5d0*(p(i,j,k+1)-p(i,j,k-1)) 
!                px= (ux(i,j)*pxi +vx(i,j)*peta +wx(i,j,k)*psig) 
!                py= (uy(i,j)*pxi +vy(i,j)*peta +wy(i,j,k)*psig) 
!                                                                         
!                ugeo(i,j,k) = - (qpr*py +drpy(i,j,k)+gpr*hy)/             &
!                     (ffc(i,j))                                          
!                vgeo(i,j,k) = (qpr*px +drpx(i,j,k) +gpr*hx)/              &
!                     (ffc(i,j))                                          
!             end do
!          end do 
!          do k=1,NK 
!             ugeo(0,j,k)= ugeo(NI,j,k) 
!             vgeo(0,j,k)= vgeo(NI,j,k) 
!             ugeo(NI+1,j,k) = ugeo(1,j,k) 
!             vgeo(NI+1,j,k) = vgeo(1,j,k) 
!          end do 
!       end do 

!       return
!       end


                                                                        
