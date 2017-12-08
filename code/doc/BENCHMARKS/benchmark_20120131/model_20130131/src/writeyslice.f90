subroutine writeyslice(jslice,counter_2d,n)
  !     ------------------------------------------------------------      
USE header,ONLY : NI,NJ,NK,ntr,s,T,rho,u,v,w,vor,xc,zc,nconsume,DL,time_seconds,dirout,rc_kind
  !     reads in a netcdf file and interpolates the data from the sigma gr
  !     onto a z level                                                    
                                                                        
#include "netcdf.inc" 
      integer NI2,NJ2,NK2,i,j,k,n
      parameter ( NI2=NI+2,NJ2=NJ+2,NK2=NK+2) 
      REAL(kind=rc_kind) :: svert(NK),Tvert(NK),rvert(NK),uvert(NK),         &
     &     vvert(NK),wvert(NK),zvert(NK),seval                          
                                                                        
      REAL(kind=4) ::  xslice(NI),zslice(NI,0:NK+1),sslice(ntr,NI,0:NK+1),        &
     &     Tslice(ntr,NI,0:NK+1),rslice(NI,0:NK+1),uslice(NI,0:NK+1),   &
     &     vslice(NI,0:NK+1),wslice(NI,0:NK+1),vorslice(NI,0:NK+1)      
!                                                                       

  integer :: itr,jslice
   REAL(kind=rc_kind) ::  rcode

                                                                        
      character(LEN=150) outname 
                      
      REAL :: time_days
      INTEGER :: counter_2d
                                                  
      integer start(3),count(3),dims(3),dims4(4),start2d(2),count4(4),  &
     &     start4(4)                                                    
                           

 integer :: iddatfile, idigit, idudx, idudy, idudz, idvbysection, idvby, idvbz, idvcon, idvcy, idvbx
  integer :: idvcz, idvc, idvdivfeddy, idvdivfreyn, idvdx, idvdz, idvd, idvfb, idvh, idvn2bar
  integer :: idvn2, idvnsq100m, idvnsq30m, idvpe,idvpsiv,idvpsiw,idvpv, idvp,idvrhbar,idvrho,idvrnk
  integer :: idvstrain,idvstress,idvstr,idvs,idvtbar,idvtemp,idvtim,idvtr,idvt,idvu,idvvb,idvvc,idvvor,idvv
  integer :: idvwb,idvwc,idvwpv,idvw,idvy,idvzsave,idvz,idwdx,idwdy,idwdz,iimday,ipos
  integer :: idvzave,idvshear


     time_days=time_seconds/86400.d0                                             
                                    
      j=jslice                                    
      do k=1,NK 
         do i=1,NI 
            zslice(i,k)= zc(i,j,k)*DL 
            xslice(i)= xc(i) 
                                                                        
!=            sslice(i,k)= s(i,j,k) +S0                                 
!=            Tslice(i,k)= T(i,j,k) +T0                                 
            do itr=1,ntr 
               sslice(itr,i,k)= s(i,j,k,n) 
               Tslice(itr,i,k)= T(i,j,k,n) 
            end do 
            rslice(i,k)= rho(i,j,k) 
            uslice(i,k)= u(i,j,k,n) 
            vslice(i,k)= v(i,j,k,n) 
            wslice(i,k)= w(i,j,k,n) 
            vorslice(i,k)= vor(i,j,k) 
                                                                        
         end do 
      end do 
!     FOR THE SAKE OF BETTER PLOTS                                      
      do i=1,NI 
         zslice(i,NK+1)= 0.d0 
         zslice(i,0)= 0.5*DL*(zc(i,j,0)+zc(i,j,1)) 
         rslice(i,NK+1)= rslice(i,NK) 
         rslice(i,0)= rslice(i,1) 
         uslice(i,NK+1)= uslice(i,NK) 
         uslice(i,0)= 2.0*uslice(i,1) -uslice(i,2) 
         vslice(i,NK+1)= 2.0*vslice(i,NK) -vslice(i,NK-1) 
         vslice(i,0)= 2.0*vslice(i,1) -vslice(i,2) 
         wslice(i,NK+1)= 2.0*wslice(i,NK) -wslice(i,NK-1) 
         wslice(i,0)= 0.0 
         vorslice(i,NK+1)= 2.0*vorslice(i,NK) -vorslice(i,NK-1) 
         vorslice(i,0)= 2.0*vorslice(i,1) -vorslice(i,2) 
         do itr=1,ntr 
            sslice(itr,i,NK+1)= sslice(itr,i,NK) 
            sslice(itr,i,0)= sslice(itr,i,1) 
            Tslice(itr,i,NK+1)= Tslice(itr,i,NK) 
            Tslice(itr,i,0)= Tslice(itr,i,1) 
         end do 
      end do 
                                                                        
!     write to a netcdf file      

      WRITE(outname,'("yslice_",I3.3,".cdf")') jslice                                      
                                                                        
                                                                        
      if (counter_2d.eq.1) then 
         idDatFile =  nccre(TRIM(dirout)//outname,NCCLOB,rcode) 
                                                                        
         dims(1) = ncddef(idDatFile,'x',NI,rcode) 
         dims(2) = ncddef(idDatFile,'z',NK+2,rcode) 
         dims(3) = ncddef(idDatFile,'time',NCUNLIM,rcode) 
                                                                        
         dims4(1) = dims(1) 
         dims4(2) = ncddef(idDatFile,'ntr',ntr,rcode) 
         dims4(3) = dims(2) 
         dims4(4) = dims(3) 
                                                                        
         idvy = ncvdef(idDatFile,'xc',NCFLOAT,1,dims(1),rcode) 
         idvz = ncvdef(idDatFile,'zc',NCFLOAT,2,dims,rcode) 
         idvs = ncvdef(idDatFile,'consump',NCFLOAT,4,dims4,rcode) 
!=         idvt = ncvdef(idDatFile,'temp',NCFLOAT,3,dims,rcode)         
         idvt = ncvdef(idDatFile,'tr',NCFLOAT,4,dims4,rcode) 
         idvrho = ncvdef(idDatFile,'rho',NCFLOAT,3,dims,rcode) 
         idvu = ncvdef(idDatFile,'u',NCFLOAT,3,dims,rcode) 
         idvv = ncvdef(idDatFile,'v',NCFLOAT,3,dims,rcode) 
         idvw = ncvdef(idDatFile,'w',NCFLOAT,3,dims,rcode) 
         idvvor = ncvdef(idDatFile,'vor',NCFLOAT,3,dims,rcode) 
                                                                        
         CALL ncendf(idDatFile,rcode) 
                                                                        
      else 
         idDatFile = ncopn(TRIM(dirout)//outname, NCWRITE,rcode) 
!         call ncredf(idDatFile,rcode)                                  
         idvs = NCVID(idDatFile, 'consump', RCODE) 
!=         idvt = NCVID(idDatFile, 'temp', RCODE)                       
         idvt = NCVID(idDatFile, 'tr', RCODE) 
         idvrho = NCVID(idDatFile, 'rho', RCODE) 
         idvu = NCVID(idDatFile, 'u', RCODE) 
         idvv = NCVID(idDatFile, 'v', RCODE) 
         idvw = NCVID(idDatFile, 'w', RCODE) 
         idvvor = NCVID(idDatFile, 'vor', RCODE) 
      endif 
                                                                        
      count(1)= NI 
      count(2)= NK+2 
      count(3)= 1 
                                                                        
      count4(1)= NI 
      count4(2)= ntr 
      count4(3)= NK+2 
      count4(4)= 1 
                                                                        
      start2d(1)= 1 
      start2d(2)= 1 
                                                                        
      start(1)= 1 
      start(2)= 1 
      start(3)= counter_2d
                                                                        
      start4(1)= 1 
      start4(2)= 1 
      start4(3)= 1 
      start4(4)= counter_2d 
                                                                        
      if (counter_2d.eq.1) then 
         CALL ncvpt(idDatFile,idvy, start(1), count(1), xslice, rcode) 
         CALL ncvpt(idDatFile,idvz, start, count, zslice, rcode) 
      endif 
      CALL ncvpt(idDatFile,idvs, start4, count4, sslice, rcode) 
!=      CALL ncvpt(idDatFile,idvT, start, count, Tslice, rcode)         
      CALL ncvpt(idDatFile,idvT, start4, count4, Tslice, rcode) 
      CALL ncvpt(idDatFile,idvrho, start, count, rslice, rcode) 
      CALL ncvpt(idDatFile,idvu, start, count, uslice, rcode) 
      CALL ncvpt(idDatFile,idvv, start, count, vslice, rcode) 
      CALL ncvpt(idDatFile,idvw, start, count, wslice, rcode) 
      CALL ncvpt(idDatFile,idvvor, start, count, vorslice, rcode) 
                                                                        
      CALL ncclos(idDatFile,rcode) 
                                                                        
      END                                           
