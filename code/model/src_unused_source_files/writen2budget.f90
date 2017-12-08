subroutine writen2budget(frame_int,step,dtf,TL,stressmax,mldn2,   &
     zbtop,zbbot,advecPV,friction,diabatic,frictop,diatop,        &
     fricbot,diabot)                                              
  !     -------------------------------------------                       
  USE dimensions
  !     write out n2budget in a netCDF file                               
#include 'netcdf.inc'                                                                        
!      INCLUDE '/sw/lib/netcdf-gfortran/include/netcdf.inc'             
                                                                        
      character*80 outname 
      integer step,frame_int,start(2),count(2),dims(2),starttim,        &
     &     counttim,dimstim                                             
      REAL(kind=rc_kind) :: dtf,TL,timday,stressmax,                         &
     &     mldn2(3),zbtop(3),zbbot(3),advecpv(3),friction(3),           &
     &     diabatic(3),frictop(3),diatop(3),fricbot(3),diabot(3)        
                                                                        
      outname= 'n2budgetmovie.cdf' 
      timday= dtf*TL*dble(step)/86400.d0 
                                                                        
      if (step.eq.0) then 
         idDatFile =  nccre(outname,NCCLOB,rcode) 
         dims(1) = ncddef(idDatFile,'levels',3,rcode) 
         dims(2) = ncddef(idDatFile,'time',NCUNLIM,rcode) 
                                                                        
         dimstim= dims(2) 
                                                                        
         idvtim = ncvdef(idDatFile,'day',NCDOUBLE,1,dimstim,rcode) 
         idvwind=ncvdef(idDatFile,'windstress',NCDOUBLE,1,dimstim,rcode) 
         idvn = ncvdef(idDatFile,'mldn2',NCDOUBLE,2,dims,rcode) 
         idvfrictop = ncvdef(idDatFile,'frictop',NCDOUBLE,2,dims,rcode) 
         idvdiatop = ncvdef(idDatFile,'diatop',NCDOUBLE,2,dims,rcode) 
         idvfricbot = ncvdef(idDatFile,'fricbot',NCDOUBLE,2,dims,rcode) 
         idvdiabot = ncvdef(idDatFile,'diabot',NCDOUBLE,2,dims,rcode) 
         idvzbtop = ncvdef(idDatFile,'zbtop',NCDOUBLE,2,dims,rcode) 
         idvzbbot = ncvdef(idDatFile,'zbbot',NCDOUBLE,2,dims,rcode) 
         idvadv = ncvdef(idDatFile,'advecpv',NCDOUBLE,2,dims,rcode) 
         idvfric = ncvdef(idDatFile,'friction',NCDOUBLE,2,dims,rcode) 
         idvdia = ncvdef(idDatFile,'diabatic',NCDOUBLE,2,dims,rcode) 
         CALL ncendf(idDatFile,rcode) 
                                                                        
      else 
         idDatFile = ncopn(outname, NCWRITE,rcode) 
         idvtim = NCVID(idDatFile,'day',RCODE) 
                                                                        
         idvwind=NCVID(idDatFile,'windstress',rcode) 
         idvn = NCVID(idDatFile,'mldn2',rcode) 
         idvfrictop = NCVID(idDatFile,'frictop',rcode) 
         idvdiatop = NCVID(idDatFile,'diatop',rcode) 
         idvfricbot = NCVID(idDatFile,'fricbot',rcode) 
         idvdiabot = NCVID(idDatFile,'diabot',rcode) 
         idvzbtop = NCVID(idDatFile,'zbtop',rcode) 
         idvzbbot = NCVID(idDatFile,'zbbot',rcode) 
         idvadv = NCVID(idDatFile,'advecpv',rcode) 
         idvfric = NCVID(idDatFile,'friction',rcode) 
         idvdia = NCVID(idDatFile,'diabatic',rcode) 
                                                                        
                                                                        
      endif 
                                                                        
      counttim =1 
      starttim= step/frame_int +1 
                                                                        
      count(1)= 3 
      count(2)= 1 
                                                                        
      start(1)= 1 
      start(2)= step/frame_int +1 
                                                                        
      CALL ncvpt(idDatFile,idvtim,starttim,counttim,timday,rcode) 
      CALL ncvpt(idDatFile,idvwind,starttim,counttim,stressmax,rcode) 
      CALL ncvpt(idDatFile,idvn,start,count,mldn2,rcode) 
      CALL ncvpt(idDatFile,idvfrictop,start,count,frictop,rcode) 
      CALL ncvpt(idDatFile,idvdiatop,start,count,diatop,rcode) 
      CALL ncvpt(idDatFile,idvfricbot,start,count,fricbot,rcode) 
      CALL ncvpt(idDatFile,idvdiabot,start,count,diabot,rcode) 
      CALL ncvpt(idDatFile,idvzbtop,start,count,zbtop,rcode) 
      CALL ncvpt(idDatFile,idvzbbot,start,count,zbbot,rcode) 
      CALL ncvpt(idDatFile,idvadv,start,count,advecpv,rcode) 
      CALL ncvpt(idDatFile,idvfric,start,count,friction,rcode) 
      CALL ncvpt(idDatFile,idvdia,start,count,diabatic,rcode) 
                                                                        
                                                                        
      CALL ncclos(idDatFile,rcode) 
                                                                        
      return 
      END                                           
