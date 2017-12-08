subroutine windstress 
  !     =======================                                           
  USE header
     implicit none
      integer i,j,k,m,step,nday, nbegin
      REAL(kind=rc_kind) :: rhoair,vspd,uspd,Cd , yr,day
      REAL(kind=rc_kind) :: taux(365),tauy(365) 
      common/realwind/ taux,tauy 
      parameter(rhoair=1.293) 
                                                                        
      TL= LEN/UL 
!     READ IN WINDS                                                     
!     ===================                                               
      if (step.eq.17501) then 
                                                                        
  101 open(unit=22, file='across_yr99winds.in') 
      open(unit=32, file='alongshlf_yr99winds.in') 
                                                                        
      if ( (dtf*TL*dble(nsteps-nbegin)/86400.d0).gt.365.d0) then 
         write(6,*) 'length of wind record is insufficient' 
         stop 
      end if 
                                                                        
      do nday=1,365 
         read(22,*) yr,day,vspd 
         if (day.ne.nday) then 
            write(6,*) 'data problems in along-shelf wind field' 
            stop 
         end if 
                                                                        
         if (abs(vspd).le.6.d0) then 
            Cd= 1.1d-3 
         else 
            Cd= 1.d-3*(0.61 + 0.063*abs(vspd)) 
         end if 
         tauy(nday) = rhoair*Cd*vspd*abs(vspd) 
      end do 
                                                                        
      do nday=1,365 
         read(32,*) yr,day,uspd 
         if (day.ne.nday) then 
            write(6,*) 'data problems in across-shelf wind field' 
            stop 
         end if 
         if (abs(uspd).le.6.d0) then 
            Cd= 1.1d-3 
         else 
            Cd= 1.d-3*(0.61 + 0.063*abs(uspd)) 
         end if 
         taux(nday) = rhoair*Cd*uspd*abs(uspd) 
!         write(6,*) 'uspd',day,uspd,taux(nday),tauy(nday)              
      end do 
      close(22); close (32) 
      write(6,*) 'wind read' 
      endif 
!     ================================                                  
                                                                        
      return 
      END                                           
