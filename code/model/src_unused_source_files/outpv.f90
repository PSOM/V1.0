subroutine outpv(vor1,vor2,vor3,pv1,pv2,pv3) 
  !----------------------------------------------------------             
  USE header, only : NI,NJ,NK,dirout,rc_kind
  !     outputs the solution as a netCDF file                             
      integer i,j,k,step,nstp,it 
      REAL(kind=rc_kind) :: vor1(NI,NJ,NK),vor2(NI,NJ,NK),vor3(NI,NJ,NK) 
      REAL(kind=rc_kind) :: pv1(NI,NJ,NK),pv2(NI,NJ,NK),pv3(NI,NJ,NK) 
#include "netcdf.inc"                                                                        
!      INCLUDE '/usr/local/include/netcdf.inc' 
                                                                        


 integer :: idigit, idudx, idudy, idudz, idvbysection, idvby, idvbz, idvcon, idvcy, idvbx
  integer :: idvcz, idvc, idvdivfeddy, idvdivfreyn, idvdx, idvdz, idvd, idvfb, idvh, idvn2bar
  integer :: idvn2, idvnsq100m, idvnsq30m, idvpe,idvpsiv,idvpsiw,idvpv, idvp,idvrhbar,idvrho,idvrnk
  integer :: idvstrain,idvstress,idvstr,idvs,idvtbar,idvtemp,idvtim,idvtr,idvt,idvu,idvvb,idvvc,idvvor,idvv
  integer :: idvwb,idvwc,idvwpv,idvw,idvy,idvzsave,idvz,idvz3,idwdx,idwdy,idwdz,iimday,ipos

  REAL(kind=rc_kind) ::  rcode,vol
  integer :: idvx,idvcon100,idFaceFile,idvuf, idvvf, idvwf, idvKf
    
  integer :: idDatFile,idv1, idv2, idv3, idpv1, idpv2, idpv3

                                                                        
      integer start(3), count(3) 
      integer dims(3) 
                                                                        
                                                                        
      DATA start /1, 1, 1/ 
                                                                        
      count(1)= NI 
      count(2)= NJ 
      count(3)= NK 
!     name the output file                                              
!                                                                       
      idDatFile =  nccre(TRIM(dirout)//'pv.cdf',NCCLOB,rcode)                                                                       
 
      dims(1) = ncddef(idDatFile,'x',NI,rcode) 
      dims(2) = ncddef(idDatFile,'y',NJ,rcode) 
      dims(3) = ncddef(idDatFile,'sigma',NK,rcode) 
                                                                        
      idv1 = ncvdef(idDatFile,'vor1',NCDOUBLE,3,dims,rcode) 
      idv2 = ncvdef(idDatFile,'vor2',NCDOUBLE,3,dims,rcode) 
      idv3 = ncvdef(idDatFile,'vor3',NCDOUBLE,3,dims,rcode) 
      idpv1 = ncvdef(idDatFile,'pv1',NCDOUBLE,3,dims,rcode) 
      idpv2 = ncvdef(idDatFile,'pv2',NCDOUBLE,3,dims,rcode) 
      idpv3 = ncvdef(idDatFile,'pv3',NCDOUBLE,3,dims,rcode) 
                                                                        
      CALL ncendf(idDatFile,rcode) 
                                                                        
      CALL ncvpt(idDatFile,idv1, start, count, vor1, rcode) 
      CALL ncvpt(idDatFile,idv2, start, count, vor2, rcode) 
      CALL ncvpt(idDatFile,idv3, start, count, vor3, rcode) 
      CALL ncvpt(idDatFile,idpv1, start, count, pv1, rcode) 
      CALL ncvpt(idDatFile,idpv2, start, count, pv2, rcode) 
      CALL ncvpt(idDatFile,idpv3, start, count, pv3, rcode) 
                                                                        
      CALL ncclos(idDatFile,rcode) 
                                                                        
      return 
      END                                           
