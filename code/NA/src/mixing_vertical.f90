      subroutine mixing_vertical(sdif,Tdif,Trdif,udif,vdif,m,step)

!     ---------------------------------------------
!     Wind Stress can be specified here
!     Kz*du/dz = -Tx/rho,   Kz*dv/dz = -Ty/rho
!     use level m 
!     computes d2s/dz2 at the cell centers.
!     fricu and fricv are u,v diffusion terms (in header) 
!     that need to be saved (if m.eq.0) for n2budget)
!    
      USE header 
      implicit none
      double precision dSin,posneg,facreverse
      integer i,j,k,m,step,nday,nw1,nw2
!     stressmax is common in header and written out in main.f
      double precision zdep,yw,yset,ymid,timday,yearday,dj,dmid, &
     &     stress,ustar,Ekdepth,ff,drhoml,rhocrit,Kdrhodz
      integer nwork,iseed,num
      parameter (nwork=10000)
      double precision da,dwork(nwork),am,sd

!     The viscosity is Kz*Kzmax 
!     Use the same value for mixing tracers and momentum
      double precision dsdzfc(0:NK),drdzfc(0:NK),&
     &     dudzfc(0:NK),dvdzfc(0:NK),&
     &     sdif(NI,NJ,NK),Tdif(ntr,NI,NJ,NK),Trdif(ntr,NI,NJ,NK),udif(NI,NJ,NK),&
     &     vdif(NI,NJ,NK),Kdudzb,Kdvdzb,Kdudzt,Kdvdzt,wzkth,fact,Cp
      double precision fac,facb,dfactor,Kzmax,facy,rhoinv,diff,&
     &     stressy,wgt,day,ts,ycenter, ymax,KzmaxTr,&
     &     ztransit,zextent,thy
      double precision rhoflux(NI,NJ)
      parameter (drhoml=0.01)
!     taux and tauy are the wind stress in N/m2
!     tau= 0.025 N/m2 =0.25 dynes/cm2 corresponds to approx wind sp of 4.5 m/s
!     K du/dz = tau / rho  and tau = rho Cd *U*U, Cd =(approx) 1.d-3
!     negative taux implies upshore wind, neg tauy implies onshore
!c      parameter (tauy= 0.025, taux=0.0)
!      parameter (taux= 0.0, tauy=0.0)
!      parameter (taux= -0.025, tauy=0.0) 

!     mldepth is the MLD

!     Kzmax= 1.d-3 (Chapman and Lentz value)  Gawarkiewicz-Chap use 2.d-3
!     RR is the bottom friction parameter in m/s
!=== const kz      
!==usual      parameter (RR= 0.0005, Kzmax= 1.d-5)
!==  TRY HIGH VALUE OF KZ
!      parameter (RR= 0.0005, Kzmax= 5.d-3)
!     for KPP use higher Kzmax
!=      parameter (RR= 0.0005, Kzmax= 5.d-4)  - tried this
!==      parameter (RR= 0.005, Kzmax= 1.d-5)  !quad drag
!     for KPP, using Kzmax= 1.d-3, turn off convective adjustment 
!==      parameter (RR= 0.0005, Kzmax= 1.d-3) (Try lower for no wind case)
!=      parameter (RR= 0.0005, Kzmax= 1.d-3) ! high value with KPP
!     no bottom friction.
!==      parameter (RR= 0.d0, Kzmax= 1.d-5) 
!forEkres run
!workswell      parameter (RR= 0.d0, Kzmax= 1.d-2, KzmaxTr=1.d-5) 
!=      parameter (RR= 0.0005d0, Kzmax= 1.d-3, KzmaxTr=1.d-5) !Tanh profile kz
!     mix buoyancy too
!!      parameter (RR= 0.0005d0, Kzmax= 1.d-2, KzmaxTr=1.d-2) !exp profile kz
!! Kzmax and KzmaxTr are made dependent on stressmax
!      parameter (RR= 0.0005d0) !exp profile in Kz
!==usual      parameter (RR= 0.0005d0, Kzmax= 1.d-5, KzmaxTr=1.d-5) 
!=ramped- calc      parameter (stressmax= 0.4d0 )

!     Kv*du/dz = r*u  at the bottom
!
!     rho must be evaluated to compute the bouyancy freq N, if 
!     we are using the KPP scheme (subroutine viscosity).


      timday= dble(step)*dtf*TL/86400.d0
      yearday= timday + yrday_start   ! simulation starts yearday 30

      call cooling(rhoflux,yearday)

! FIND MLD - needed if mixing depth= mld/3
      do j=1,NJ
         do i=1,NI
!     Find MLD mld
            rhocrit =rho(i,j,NK)+ drhoml
            do k=NK-1,1,-1
               if ( rho(i,j,k).ge.rhocrit ) then 
                  mld(i,j)= -DL*( (zc(i,j,k)-zc(i,j,k+1))* &
     &            (rhocrit -rho(i,j,k+1))/(rho(i,j,k)-rho(i,j,k+1)) &
     &            + zc(i,j,k+1)  )
                  goto 51
               end if
            end do
!     if the do loop is computed, and the ML base is not found, set it to D
            mld(i,j)= -zf(i,j,0)*DL
 51         continue
         end do
      end do
!===================found MLD---------------

!      TL= LEN/UL
      timday= dble(step)*dtf*TL/86400.d0
      yearday= timday + yrday_start   ! simulation starts yearday 30
!==   No wind
!      do j=1,NJ
!         stressx(j) =0.0 
!      end do
!      goto 201
!==
!=      goto 101

!      stressmax=  0.1d0
!=      stressmax=  0.025d0
!      stressmax=  0.0d0

!     Sinusoidal wind 
!     =======================
!c     Ramped in Time
!     --------------
!     upfront after 25 days
!      if (timday.ge.25.d0) then
!         stressmax= -0.1d0
!     dnfront before 20 days
!      else if (timday.le.20.d0) then
!         stressmax= 0.15d0
!      else 
!         stressmax = -0.1d0*(timday-20.d0)/5.d0
!     &        + 0.15d0*(25.d0-timday)/5.d0
!      end if

!     FOR MODEL START AT YRDAY 50 - RUN UPTO YRDAY 120
!     upfront -0.1 after yr day 85, start yrday=50, so after 35 days
      if (yearday.ge.85.d0) then
         stressmax= -0.1d0
!     dnfront 0.2 prior to yrday 70, i.e. before day 20
      else if (yearday.le.70.d0) then
         stressmax= 0.2d0
      else 
         stressmax = -0.1d0*(yearday-70.d0)/15.d0 &
     &        + 0.2d0*(85.d0-yearday)/15.d0
      end if

!!!USE THIS     for constant
!      stressmax=0.2d0           !constant downfront

!=========added temporarily  3.12.12
!     CONSTANT wind + perturb
!      iseed= step
!      num=1
!      am=  0.d0
!      sd= 0.1d0
!      da= 0.d0
!      call rannw (am, sd, iseed, da, num, dwork, nwork)
!      stressmax= 0.1d0*(1.d0 +da)
!=========================================

!      Kzmax= 1.d-2
!stop Oct 3,2011      Kzmax= dmax1(1.d-5, 1.d-2* (dabs(stressmax) / 0.2d0))
      Kzmax= dmax1(1.d-3, 1.d-2* (abs(stressmax) / 0.2d0))
      KzmaxTr= Kzmax

!     DONN-WIND
!     ===============
!=      stressmax=0.25d0    ! by=0.45 x1e-7 , H=300m,  r=1
!      stressmax= 0.167      ! by=0.3 x1e-7 , H=300m,  r=1
!      stressmax=0.2d0           ! by=0.45 x1e-7 , H=250m,  r=1, tau=0.173
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

!        
!     =====================
!     Windstress - Tanh profile in y
!-      do j=1,njf1
!-         stressx(j)=  stressmax* 0.5*(1.+ 
!-     &        dtanh( pi*dble(j-njf1/2)/dble(njf1/2)) )
!-      end do
!-      do j=njf1+1,njf2
!-         stressx(j)= stressmax * 1.d0
!-      end do
!-      dj= 0.5d0*dble(NJ-njf2)
!-      dmid= dble(njf2)+dj
!-      do j=njf2+1,NJ
!-         stressx(j)= stressmax* 0.5*(1.-
!-     &        dtanh( pi*(dble(j)-dmid)/dj ) )
!-      end do
!-      goto 201
!     Sin wind
!     Extent of wind (where wind is flat)
!      nw1=njf1/2                 ! 5 fronts , NJ=640
!      nw2=njf4+(nj-njf4)/2       ! 5 fronts , NJ=640
!!      nw1=njf1/4                 ! 3 fronts, NJ=480
!!      nw2=njf4+3*(nj-njf4)/4     ! 3 fronts, NJ=480
!     Fix the extent of region where there is a wind-stress curl
      nw1=0.1*nj
      nw2=0.9*nj
      nw1=48
      nw2= 432

      ycenter=yc(nw1)
      ymax = ycenter-yc(1+10)
      do j=1,nw1
         diff = (yc(j)-ycenter)
         if (abs(diff).gt.ymax) then
            stressx(j)= 0.d0
         else
            stressx(j)= stressmax*dcos((diff/ymax)*0.5*PI)
         end if
      end do
      do j=nw1+1,nw2
         stressx(j)= stressmax
      end do
      ycenter=yc(nw2)
      ymax = yc(nj-10)-ycenter
      do j=nw2 +1,NJ
         diff = (yc(j)-ycenter)
         if (abs(diff).gt.ymax) then
            stressx(j)= 0.d0
         else
            stressx(j)= stressmax*dcos((diff/ymax)*0.5*PI)
         end if
      end do
!      write(6,*) 'njf',njf1,njf2,njf3,njf4
!      do j=1,NJ
!         write(60,*) j,stressx(j)
!      end do
!      write(6,*) 'stop in diffusion'
!      stop
      goto 201
!
!     FLAT WIND - tapered at boundaries (Tanh profile)
!     --------------------------------
 101   continue
      yw= 10.d0
      ymid= 0.5* dble(NJ)
      yset =20.d0
      do j=1,NJ/2
         yc(j)= dble(j)-0.5
         stressx(j) = stressmax* 0.5* &
     &        (1.d0 +tanh((( yc(j)-yset)/yw)*PI ))
      end do

      yset= NJ-20
      do j=NJ/2+1,NJ
         yc(j)= dble(j)-0.5
         stressx(j) = stressmax*0.5* &
     &        (1.d0 +tanh((-( yc(j)-yset)/yw)*PI ))
      end do

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

 201  stressy= 0.d0
      
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
!          if (yc(j).gt.yfront) then
!            facy = 0.1
!         else
!            facy =  dmax1((yfront - yc(j))/yfront, 0.1d0)
!         end if
         do i=1,NI
            do k=1,NK-1
               wzkth= wz(i,j,k)
               dsdzfc(k)= wzkth*(s(i,j,k+1,m)-s(i,j,k,m))
               dudzfc(k)= wzkth*(u(i,j,k+1,m)-u(i,j,k,m))
               dvdzfc(k)= wzkth*(v(i,j,k+1,m)-v(i,j,k,m))
               drdzfc(k)= wzkth*(rho(i,j,k+1)-rho(i,j,k))
            end do
            dsdzfc(0)= 0.d0
            dudzfc(0)= 0.d0
            dvdzfc(0)= 0.d0
            drdzfc(0)= 0.0
! changed on Aug 7,2002, for wind stress. If no wind, dudzfc(NK)=0
            dsdzfc(NK)= 0.0
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
!           
!     Const viscosity, then, evalrho & viscosity need not be called.
!     call viscosity (routine below) for KPP mixing, but set Kzmax to 1.d-3
!=            call viscosity(dudzfc,dvdzfc,drdzfc,i,j)
!= else... use const viscosity Kz=1
!= const visc. kz=1
!     Use Tanh profile for Kz, transition between 120 and 220m, i.e. center
!     at 170m depth,  for mld=100m
!     Center at ztransit, zextent is the vertical extent of the variation
!     
!=            zextent= 100.d0
!=            ztransit= -mldepth -20. -0.5*zextent
!            ztransit= -170.d0 for mldepth=100
!            ztransit= -270.d0 for mldepth=100
!     Transition starts zextent above ztransit and ends zextent below 
!     at ztransit, Value is half the max value
!     Set ztransit as the Ekman depth
            stress= dsqrt(stressx(j)*stressx(j)+ stressy*stressy)
            ustar = dsqrt( stress /R0 )
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
               thy = (1.d0 +tanh(((zdep-ztransit)/zextent)*PI) )&
     &              *0.5d0 
!=               if ((j.eq.(NJ/2)).and. (i.eq.(NI/2))) write(6,*) 
!     thy varies from +1 to 0
!     EXP
!== Ekman depth               Kz(i,j,k)= dmax1(1.d-3,dexp(zdep/Ekdepth))
! use mld/3 or mld/4.6 as the e-folding depth
              Kz(i,j,k)= dmax1(1.d-3,dexp(4.6d0*zdep/mld(i,j)))
! use Ekdep/2 as the e-folding depth
!               Kz(i,j,k)= dmax1(1.d-3,dexp(2.d0*zdep/Ekdepth))
!=TANH used in stratpaper                Kz(i,j,k)= dmax1(0.01d0,thy)
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
!=            Kdudzb= facb*u(i,j,1,m)*dabs(u(i,j,1,m))
!=            Kdvdzb= facb*v(i,j,1,m)*dabs(v(i,j,1,m))
!     Linear Drag
            Kdudzb= facb*u(i,j,1,m)
            Kdvdzb= facb*v(i,j,1,m)

!      fac= 1.d0/(UL*DL*delta)
!      fact = DL/UL
!            rhoinv = 1.d0/rho(i,j,NK)
            rhoinv = 1.d0/R0
!            Kdudzt= taux(nday)*rhoinv*fact
!            Kdvdzt= tauy(nday)*rhoinv*fact
            Kdudzt= stressx(j)*rhoinv*fact
            Kdvdzt= stressy*rhoinv*fact
!     fac= 1.d0/(UL*DL*delta)
            k=1 
!          -----
            dfactor=  fac*Jac(i,j,k)*wz(i,j,k)
            sdif(i,j,k)= dfactor*KzmaxTr&
     &           *(dsdzfc(k)*Kz(i,j,k)- dsdzfc(k-1)*Kz(i,j,k-1) )
            udif(i,j,k)= dfactor&
     &           *(dudzfc(k)*Kz(i,j,k)*Kzmax -Kdudzb )
            vdif(i,j,k)= dfactor&
     &           *(dvdzfc(k)*Kz(i,j,k)*Kzmax -Kdvdzb )
!=     
!            fricu(i,j,k)= 
!     &           (dudzfc(k)*Kz(i,j,k)*Kzmax -Kdudzb )
!     &           *wz(i,j,k)*UL/(DL*DL)
!
            do k=2,NK-1
               dfactor=  fac*KzmaxTr*Jac(i,j,k)*wz(i,j,k)
               sdif(i,j,k)= dfactor&
     &              *(dsdzfc(k)*Kz(i,j,k)- dsdzfc(k-1)*Kz(i,j,k-1) )
               dfactor=  fac*Kzmax*Jac(i,j,k)*wz(i,j,k)
               udif(i,j,k)= dfactor&
     &              *(dudzfc(k)*Kz(i,j,k)- dudzfc(k-1)*Kz(i,j,k-1) )
               vdif(i,j,k)= dfactor&
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
!     rhoflux is in kg/m^2/s  Mult. by DL to make in  (kg/m)/s
            Kdrhodz= rhoflux(i,j)*DL
!     fac= 1.d0/(UL*DL*delta)
            dfactor=  fac*Jac(i,j,k)*wz(i,j,k)
            sdif(i,j,k)= dfactor&
     &           *(Kdrhodz - dsdzfc(k-1)*Kz(i,j,k-1)*KzmaxTr )
            udif(i,j,k)= dfactor&
     &           *(Kdudzt  - dudzfc(k-1)*Kz(i,j,k-1)*Kzmax )
            vdif(i,j,k)= dfactor&
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
                  fricb(i,j,k)= -(gpr*10.d0*sdif(i,j,k)*&
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
!      write(6,*) 'stop in diffusion'
!      stop


      return
      end




      subroutine viscosity(dudz,dvdz,drdz,i,j)
!     ---------------------------------------------
!     Compute a Richardson number based vertical viscosity 
!     coefficient Kz using the final velocities u,v and density
!     with the n=0 subscript
!     Kz is situated at k cell faces
!     The algorithm is from Rev Geophys., p 373, 1994. Large et al.
!     
!     Assumes that rho is evaluated just prior to this call

       USE header
      implicit none
      integer i,j,k

      integer n1
      parameter (n1=3)
      double precision fac,bvfreq,grho,RiCr,Ri,shearl
      double precision dudz(0:NK),dvdz(0:NK),drdz(0:NK)
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
         shearl= (dudz(k)*dudz(k) + dvdz(k)*dvdz(k))*fac
         if (shearl==0.d0) then
!         if (shear.le.1.d-12) then
            Ri= 100.d0
         else
            Ri= bvfreq/shearl
         end if
!     unstable density profile => Ri < 0
         if (Ri.le.0.d0) then
            Kz(i,j,k)= 1.d0
         else if (Ri.lt.RiCr) then
            Kz(i,j,k)= (1.d0 - (Ri*Ri)/(RiCr*RiCr))**n1
         else
            Kz(i,j,k)= 0.d0
         endif
!         if ((i.eq.NI/2).and.(j.eq.NJ/2)) write(6,*) Ri ,bvfreq,shear,
!     &        Kz(i,j,k)
      end do

!     The value of Kz to be used is  Kz(k) *Kzmax

      return
      end

