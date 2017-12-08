subroutine diffusion(sdif,Tdif,Trdif,udif,vdif,m,step) 
  !     ---------------------------------------------                     
  USE header
  !     Wind Stress can be specified here                                 
  !     Kz*du/dz = -Tx/rho,   Kz*dv/dz = -Ty/rho                          
  !     use level m                                                       
  !     computes d2s/dz2 at the cell centers.                             
  !     fricu and fricv are u,v diffusion terms (in header)               
  !     that need to be saved (if m.eq.0) for n2budget)                   

  implicit none 
      REAL(kind=rc_kind) :: dSin,posneg,facreverse,ywind,edgeice,edgeshelf
      integer i,j,k,m,step,nday 
!     stressmax is common in header and written out in main.f           
      REAL(kind=rc_kind) :: zdep,yw,yset,ymid,timday,                        &
     &     stress,ustar,Ekdepth,ff                                      
                                                                        
!     The viscosity is Kz*Kzmax                                         
!     Use the same value for mixing tracers and momentum                
      REAL(kind=rc_kind) :: dsdzfc(0:NK),dTdzfc(0:NK),drdzfc(0:NK),          &
     &     dudzfc(0:NK),dvdzfc(0:NK),                                   &
     &     sdif(NI,NJ,NK),Tdif(NI,NJ,NK),Trdif(ntr,NI,NJ,NK),           &
     &     udif(NI,NJ,NK),                                              &
     &     vdif(NI,NJ,NK),Kdudzb,Kdvdzb,Kdudzt,Kdvdzt,wzkth,fact,Cp     
      REAL(kind=rc_kind) :: fac,facb,dfactor,Kzmax,facy,rhoinv,diff,      &
     &     stressy,wgt,day,ts,ycenter, ymax,KzmaxTr,                    &
     &     ztransit,zextent,thy                                         
!     taux and tauy are the wind stress in N/m2                         
!     tau= 0.025 N/m2 =0.25 dynes/cm2 corresponds to approx wind sp of 4
!     K du/dz = tau / rho  and tau = rho Cd *U*U, Cd =(approx) 1.d-3    
!     negative taux implies upshore wind, neg tauy implies onshore      
!c      parameter (tauy= 0.025, taux=0.0)                               
!      parameter (taux= 0.0, tauy=0.0)                                  
!      parameter (taux= -0.025, tauy=0.0)                               
                                                                        
!     mldepth is the MLD                                                
                                                                        
!     Kzmax= 1.d-3 (Chapman and Lentz value)  Gawarkiewicz-Chap use 2.d-
!     RR is the bottom friction parameter in m/s                        
!=== const kz                                                           
!==usual      parameter (RR= 0.0005, Kzmax= 1.d-5)                      
!==  TRY HIGH VALUE OF KZ                                               
!      parameter (RR= 0.0005, Kzmax= 5.d-3)                             
!     for KPP use higher Kzmax                                          
!=      parameter (RR= 0.0005, Kzmax= 5.d-4)  - tried this              
!==      parameter (RR= 0.005, Kzmax= 1.d-5)  !quad drag                
!     for KPP, using Kzmax= 1.d-3, turn off convective adjustment       
!==      parameter (RR= 0.0005, Kzmax= 1.d-3) (Try lower for no wind cas
!=      parameter (RR= 0.0005, Kzmax= 1.d-3) ! high value with KPP      
!     no bottom friction.                                               
!==      parameter (RR= 0.d0, Kzmax= 1.d-5)                             
!forEkres run                                                           
!workswell      parameter (RR= 0.d0, Kzmax= 1.d-2, KzmaxTr=1.d-5)       
!=      parameter (RR= 0.0005d0, Kzmax= 1.d-3, KzmaxTr=1.d-5) !Tanh prof
!     mix buoyancy too                                                  
                                                            !Tanh profil
!!      parameter (RR= 0.0005d0, Kzmax= 1.d-3, KzmaxTr=1.d-3) 
      parameter (Kzmax= 1.d-3, KzmaxTr=1.d-3) 
!      parameter (RR= 0.0005d0, Kzmax= 1.d-2, KzmaxTr=1.d-2) !Tanh profi
!==usual      parameter (RR= 0.0005d0, Kzmax= 1.d-5, KzmaxTr=1.d-5)     
!=ramped- calc      parameter (stressmax= 0.4d0 )                       
                                                                        
!     Kv*du/dz = r*u  at the bottom                                     
!                                                                       
!     rho must be evaluated to compute the bouyancy freq N, if          
!     we are using the KPP scheme (subroutine viscosity).               
                                                                        
!      TL= LEN/UL 

!   RR=0.0005d0                                                      

      timday= dble(step)*dtf*TL/86400.d0 
!==   No wind                                                           
!      do j=1,NJ                                                        
!         stressx(j) =0.0                                               
!      end do                                                           
!      goto 201                                                         
!==                                                                     
!=      goto 101                                                        
                                                                        
!      stressmax=  0.025d0                                              
!      stressmax=  -0.1d0     
      stressmax=  0.1d0                                                
      ! stressmax=  0.0d0 
                                                                        
!     Sinusoidal wind                                                   
!     =======================                                           
!c     Ramped in Time                                                   
!     --------------                                                    
!c     upfront after 25 days                                            
!      if (timday.ge.25.d0) then                                        
!         stressmax= -0.1d0                                             
!c     dnfront before 20 days                                           
!      else if (timday.le.20.d0) then                                   
!        stressmax= 0.2d0                                               
!      else                                                             
!         stressmax = -0.1d0*(timday-20.d0)/5.d0                        
!     &        + 0.2*(25.d0-timday)/5.d0                                
!      end if                                                           
!c     =====================                                            
!     Sin wind - 10 day positive, 10 day negative                       
!     Sin wind - 20 day positive, 20 day negative                       
!==      stressmax= 0.2d0*dSin(PI*timday/10.d0)                         
!     Square wind - 10 day on/off                                       
!      posneg= 0.2d0*dSin(PI*timday/10.d0)                              
!      if (posneg.gt.0.d0) then                                         
!         stressmax= 0.2d0                                              
!      else                                                             
!         stressmax= 0.d0                                               
!      end if                                                           
                                                                        
                                                                        
      ycenter = 0.5*(yc(NJ)+yc(1)) 
      ymax = yc(NJ-10)-ycenter 
!      do j=1,NJ/2
!         stressx(j)= stressmax
!         diff = (yc(j)-ycenter) 
!         if (dabs(diff).gt.ymax) then 
!            stressx(j)= 0.d0 
!         else 
!            stressx(j)= stressmax*dcos((diff/ymax)*0.5*PI) 
!         end if 
!      end do 
      edgeshelf=0.03
      ywind= 90.0 ! 25 km from coast
      !do j=1,NJ
        !  stressx(j)= stressmax*0.5*(tanh(edgeshelf*(yc(j)-ywind)*PI)+1.)
! stressx(j)= 0+0.*stressmax*(1.-(0.5*(tanh(edgeshelf*(yc(j)-ywind-40)*PI)+1.)+(1.-0.5*(tanh(edgeshelf*(yc(j)-ywind+40)*PI)+1.))))
      !end do
!      ywind= yc(NJ)-25.0 !no ice  - wind ends 25 km from bndry
!      ywind= 200.0     ! ice  - wind ends close to shelf
! reduce tightnness by 3/4
!      edgeice=0.03
!      do j=NJ/3,NJ
!         stressx(j)= -stressmax*0.5*(tanh(edgeice*(yc(j)-ywind)*PI)-1.d0)
!      end do
!      do j=1,NJ
!         write(60,*) j,stressx(j)                                       
!      end do
      goto 201 
!                                                                       
!     FLAT WIND - tapered at boundaries (Tanh profile)                  
!     --------------------------------                                  
  101  continue 
      yw= 10.d0 
      ymid= 0.5* dble(NJ) 
      yset =20.d0 
      do j=1,NJ/2 
         yc(j)= dble(j)-0.5 
         stressx(j) = stressmax* 0.5*                                   &
     &        (1.d0 +tanh((( yc(j)-yset)/yw)*PI ))                      
      end do 
                                                                        
      yset= NJ-20 
      do j=NJ/2+1,NJ 
         yc(j)= dble(j)-0.5 
         stressx(j) = stressmax*0.5*                                    &
     &        (1.d0 +tanh((-( yc(j)-yset)/yw)*PI ))                     
      end do 
!==                                                                     
                                                                        
      do j=1,NJ 
         if (abs(stressx(j)-stressmax).lt.0.001) stressx(j)=stressmax 
         if (abs(stressx(j)).lt.0.001) stressx(j)=0.d0 
      end do 
                                                                        
                                                                        
!     TAPER WINDS at boundaries                                         
!=      do j=1,10                                                       
!=         stressx(j)= 0.d0                                             
!=         stressx(NJ-j+1) = 0.d0                                       
!=      end do                                                          
!     =      do j= 11,30                                                
!      do j=1,20                                                        
!=         stressx(j) =  (dble(j-10)/20.d0)*stressmax                   
!=         stressx(NJ-j+1) = (dble(j-10)/20.d0)*stressmax               
!      end do                                                           
!      do j=1,NJ                                                        
!         write(6,*) j,stressx(j)                                       
!      end do                                                           
!      stop                                                             
                                                                        
  201 stressy= 0.d0 
                                                                        
      ts=  dble(step)/4000.d0 
!      stressx = 0.025*dsin(ts*2.d0*PI)                                 
!=      stressx = 0.075*dsin(ts*2.d0*PI)                                
!===Not needed for constant viscosity case, needed for KPP              
!=      call evalrho(rho,m)                                             
                                                                        
!     for linear drag, facb = RR*DL, for quadratic drag, facb = RR*DL*UL
!     Linear Drag                                                       
      facb = RR*DL 
!     Quadratic drag                                                    
!==      facb = RR*DL*UL                                                
!     --------                                                          
      fac= 1.d0/(UL*DL*delta) 
      fact = DL/UL 
                                                                        
      do j=1,NJ 
          if (yc(j).gt.yfront) then 
            facy = 0.1 
         else 
            facy =  max((yfront - yc(j))/yfront, 0.1d0) 
         end if 
         do i=1,NI 
            do k=1,NK-1 
               wzkth= wz(i,j,k) 
               dsdzfc(k)= wzkth*(s(i,j,k+1,m)-s(i,j,k,m)) 
               dTdzfc(k)= wzkth*(T(i,j,k+1,m)-T(i,j,k,m)) 
               dudzfc(k)= wzkth*(u(i,j,k+1,m)-u(i,j,k,m)) 
               dvdzfc(k)= wzkth*(v(i,j,k+1,m)-v(i,j,k,m)) 
               drdzfc(k)= wzkth*(rho(i,j,k+1)-rho(i,j,k)) 
            end do 
            dsdzfc(0)= 0.d0 
            dTdzfc(0)= 0.d0 
            dudzfc(0)= 0.d0 
            dvdzfc(0)= 0.d0 
            drdzfc(0)= 0.0 
! changed on Aug 7,2002, for wind stress. If no wind, dudzfc(NK)=0      
            dsdzfc(NK)= 0.0 
            dTdzfc(NK)= 0.0 
            dudzfc(NK)= 0.0 
            dvdzfc(NK)= 0.0 
            drdzfc(NK)= 0.0 
!            dsdzfc(NK)= dsdzfc(NK-1)                                   
!            dudzfc(NK)= dudzfc(NK-1)                                   
!            dvdzfc(NK)= dvdzfc(NK-1)                                   
!            drdzfc(NK)= drdzfc(NK-1)                                   
!                                                                       
!     I don't need dudzfc(k=0),dvdzfc(k=0), but do need dsdzfc(k=0).    
!     Kz(k=0) is used only in the s,T equations.                        
!                                                                       
                                                                        
!     Const viscosity, then, evalrho & viscosity need not be called.    
!     call viscosity (routine below) for KPP mixing, but set Kzmax to 1.
!=            call viscosity(dudzfc,dvdzfc,drdzfc,i,j)                  
!= else... use const viscosity Kz=1                                     
!= const visc. kz=1                                                     
!     Use Tanh profile for Kz, transition between 120 and 220m, i.e. cen
!     at 170m depth,  for mld=100m                                      
!     Center at ztransit, zextent is the vertical extent of the variatio
!                                                                       
!=            zextent= 100.d0                                           
!=            ztransit= -mldepth -20. -0.5*zextent                      
!            ztransit= -170.d0 for mldepth=100                          
!            ztransit= -270.d0 for mldepth=100                          
!     Transition starts zextent above ztransit and ends zextent below   
!     at ztransit, Value is half the max value                          
!     Set ztransit as the Ekman depth                                   
            stress= sqrt(stressx(j)*stressx(j)+ stressy*stressy) 
            ustar = sqrt( stress /R0 ) 
            ff= ffc(i,j)*FPAR 
            Ekdepth= 0.4d0*ustar/ff 
!     Set zextent                                                       
            zextent= 0.5d0*Ekdepth 
            ztransit = -Ekdepth 
!            zextent= 0.5*mldepth                                       
!            ztransit= -(mldepth -0.25*mldepth)                         
            do k=NK,0,-1 
               Kz(i,j,k)= 1.d0 
               zdep = zf(i,j,k)*DL 
!==     Tanh profile in Kz,  Use KzmaxTr=1.d-3                          
               thy = (1.d0 +tanh(((zdep-ztransit)/zextent)*PI))         &
     &              *0.5d0                                              
!     thy varies from +1 to 0                                           
               Kz(i,j,k)= max(0.01d0,thy) 
!               if ((j.eq.(NJ/2)).and. (i.eq.(NI/2)))                   
!     &              write(6,*) k,Kz(i,j,k)                             
!==                                                                     
!===  end of Tanh profile in Kz                                         
!====CONSTANT Kz introduced                                             
!               Kz(i,j,k)= 1.d-2                                        
!====CONST Kz, Pr=1                                                     
                                                                        
!exp               Kz(i,j,k) = dmax1(0.001d0,(dexp(-zdep/25.d0)))       
!     to taper Kz across channel in accordance with wind                
!==               Kz(i,j,k) = dmax1(0.001d0,(dexp(-zdep/25.d0)*         
!==     &              (stressx(j)/stressmax)**2) )                     
!               write(6,*) k,zdep,Kz(i,j,k)                             
!               Kz(i,j,k)=  10.0*facy                                   
!               if ((j.eq.10).and.(i.eq.10))                            
!=               if ((j.eq.NJ/2).and. (i.eq.NI/2))                      
!     &              write(6,*) k,zdep,thy,Kz(i,j,k)                    


            end do 
                                                                        
                                                                        
!     Quadratic drag                                                    
!=            Kdudzb= facb*u(i,j,1,m)*abs(u(i,j,1,m))                  
!=            Kdvdzb= facb*v(i,j,1,m)*abs(v(i,j,1,m))                  
!     Linear Drag                                                       
            Kdudzb= facb*u(i,j,1,m) 
            Kdvdzb= facb*v(i,j,1,m) 
                                                                        
!      fac= 1.d0/(UL*DL*delta)                                          
!      fact = DL/UL                                                     
            rhoinv = 1.d0/rho(i,j,NK) 
!            Kdudzt= taux(nday)*rhoinv*fact                             
!            Kdvdzt= tauy(nday)*rhoinv*fact                             
            Kdudzt= stressx(j)*rhoinv*fact 
            Kdvdzt= stressy*rhoinv*fact 
!     fac= 1.d0/(UL*DL*delta)                                           
            k=1 
!          -----                                                        
            dfactor=  fac*Jac(i,j,k)*wz(i,j,k) 

            sdif(i,j,k)= dfactor*KzmaxTr                                &
     &           *(dsdzfc(k)*Kz(i,j,k)- dsdzfc(k-1)*Kz(i,j,k-1) )       
            Tdif(i,j,k)= dfactor*KzmaxTr                                &
     &           *(dTdzfc(k)*Kz(i,j,k)- dTdzfc(k-1)*Kz(i,j,k-1) )       
            udif(i,j,k)= dfactor                                        &
     &           *(dudzfc(k)*Kz(i,j,k)*Kzmax -Kdudzb )                  
            vdif(i,j,k)= dfactor                                        &
     &           *(dvdzfc(k)*Kz(i,j,k)*Kzmax -Kdvdzb )                  
!=                                                                      
!            fricu(i,j,k)=                                              
!     &           (dudzfc(k)*Kz(i,j,k)*Kzmax -Kdudzb )                  
!     &           *wz(i,j,k)*UL/(DL*DL)                                 
!                                                                       
            do k=2,NK-1 
               dfactor=  fac*KzmaxTr*Jac(i,j,k)*wz(i,j,k) 
               sdif(i,j,k)= dfactor                                     &
     &              *(dsdzfc(k)*Kz(i,j,k)- dsdzfc(k-1)*Kz(i,j,k-1) )    
               Tdif(i,j,k)= dfactor                                     &
     &              *(dTdzfc(k)*Kz(i,j,k)- dTdzfc(k-1)*Kz(i,j,k-1) )    
               dfactor=  fac*Kzmax*Jac(i,j,k)*wz(i,j,k) 
               udif(i,j,k)= dfactor                                     &
     &              *(dudzfc(k)*Kz(i,j,k)- dudzfc(k-1)*Kz(i,j,k-1) )    
               vdif(i,j,k)= dfactor                                     &
     &              *(dvdzfc(k)*Kz(i,j,k)- dvdzfc(k-1)*Kz(i,j,k-1) )    
                                                                        
!=                                                                      
!               fricu(i,j,k)= Kzmax*                                    
!     &              (dudzfc(k)*Kz(i,j,k)- dudzfc(k-1)*Kz(i,j,k-1) )    
!     &              *wz(i,j,k)*(UL/DL/DL)                              
!==                                                                     
!            if ((m.eq.0).and.(i.eq.NI/2).and.(j.eq.NJ/2))              
!     &                 write(6,*) k, 'dudz',dudzfc(k)                  
            end do 
                                                                        
            k=NK 
!           ----                                                        
!     fac= 1.d0/(UL*DL*delta)                                           
            dfactor=  fac*Jac(i,j,k)*wz(i,j,k) 
            sdif(i,j,k)= dfactor*KzmaxTr                                &
     &           *(dsdzfc(k)*Kz(i,j,k)- dsdzfc(k-1)*Kz(i,j,k-1) )       
            Tdif(i,j,k)= dfactor*KzmaxTr                                &
     &           *(dTdzfc(k)*Kz(i,j,k)- dTdzfc(k-1)*Kz(i,j,k-1) )       
            udif(i,j,k)= dfactor                                        &
     &           *(Kdudzt - dudzfc(k-1)*Kz(i,j,k-1)*Kzmax )             
            vdif(i,j,k)= dfactor                                        &
     &           *(Kdvdzt - dvdzfc(k-1)*Kz(i,j,k-1)*Kzmax )             
                                                                        
!=                                                                      
!            fricu(i,j,k)= (Kdudzt - dudzfc(k-1)*Kz(i,j,k-1)*Kzmax )    
!     &           *wz(i,j,k)*(UL/DL/DL)                                 
         end do 
      end do 
                                                                        
!     save udif,vdif if m.eq.0                                          
      if (m.eq.0) then 
         do k=1,NK 
            do j=1,NJ 
               do i=1,NI 
!     facreverse undoes the part done for solving non-dim equations     
!     fricu, fricv in m/s^2                                             
                  facreverse=UL*UL/(Jac(i,j,k)*LEN) 
                  fricu(i,j,k)= udif(i,j,k)*facreverse 
                  fricv(i,j,k)= vdif(i,j,k)*facreverse 
                  fricw(i,j,k)= 0.d0 
                  fricb(i,j,k)= -(gpr*10.d0*sdif(i,j,k)*                &
     &                 UL)/(R0*Jac(i,j,k)*LEN)                          
!                  if ((i.eq.NI/2).and.(j.eq.NJ/2))                     
!     &                 write(6,*) k, fricu(i,j,k)                      
               end do 
            end do 
         end do 
      end if 
!      do k=1,NK                                                        
!         write(6,*) 'vert',k,zc(i,NJ/2,k),Kz(i,NJ/2,k)                 
!      end do                                                           
!      do j=1,NJ                                                        
!         write(6,*) 'hor',j,Kz(NI/2,j,NK)                              
!      end do                                                           
!c      stop                                                            
                                                                        
      return 
      END                                           
                                                                        
                                                                        
                                                                        
                                                                        
      subroutine viscosity(dudz,dvdz,drdz,i,j) 
!     ---------------------------------------------                     
        USE header
!     Compute a Richardson number based vertical viscosity              
!     coefficient Kz using the final velocities u,v and density         
!     with the n=0 subscript                                            
!     Kz is situated at k cell faces                                    
!     The algorithm is from Rev Geophys., p 373, 1994. Large et al.     
!                                                                       
!     Assumes that rho is evaluated just prior to this call             
                                                                        
      integer i,j,k 
                                                                        
      integer n1 
      parameter (n1=3) 
      REAL(kind=rc_kind) :: fac,bvfreq,grho,RiCr,Ri,vshear 
      REAL(kind=rc_kind) :: dudz(0:NK),dvdz(0:NK),drdz(0:NK) 
      parameter (RiCr= 0.7d0) 
                                                                        
      grho= 9.81/R0 
!     fac= DL/(UL*UL)                                                   
      DLinv = 1.0/DL 
      fac = UL*UL/(DL*DL) 
!                                                                       
!     Set kz(k=0)= 1. It is required only for the s,T equations, since  
!     for the momentum equations,  K*du/dz= ru.                         
!     Kz is not needed at k=0 and NK since stress bcs are used.         
      Kz(i,j,0)= 0.d0 
      Kz(i,j,NK)= 0.d0 
!      if  (i.eq.16) write(100,*) 'j = ', j                             
      do k=1,NK-1 
         bvfreq=  -grho*drdz(k)*DLinv 
!     BVfreq is N**2 and is in s^-2 if re-dim by DLinv                  
         vshear= (dudz(k)*dudz(k) + dvdz(k)*dvdz(k))*fac 
         if (vshear.eq.0.d0) then 
!         if (vshear.le.1.d-12) then                                    
            Ri= 100.d0 
         else 
            Ri= bvfreq/vshear 
         end if 
!     unstable density profile => Ri < 0                                
         if (Ri.le.0.d0) then 
            Kz(i,j,k)= 1.d0 
         else if (Ri.lt.RiCr) then 
            Kz(i,j,k)= (1.d0 - (Ri*Ri)/(RiCr*RiCr))**n1 
         else 
            Kz(i,j,k)= 0.d0 
         endif 
!         if ((i.eq.NI/2).and.(j.eq.NJ/2)) write(6,*) Ri ,bvfreq,vshear,
!     &        Kz(i,j,k)                                                
      end do 
                                                                        
!     The value of Kz to be used is  Kz(k) *Kzmax                       
                                                                        
      return 
      END                                           
