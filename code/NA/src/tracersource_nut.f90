       subroutine tracersource(n,step,dtimel)
                                  
       USE header
       integer i,j,k,n,m,step,int_yd
       double precision dtimel,lightmld,light0,l0min,l0max,weight,  &
     &     KzmaxTr,chl2c,chlinteg,chl,dzsum,dz,             &
     &     alpha,death,wsink,dpdz,deltime,depthmix,drhoml,timday, &
     &     yearday,sw,swpar,swconvert,ydl1,ydl2,wphy,mort,KChl,Kwater,  &
     &     rhocrit,mumax,temp,netprod(NK),oldT(NK),fluxdiv(NK),   &
     &     Tave,trseed,zdep,fac
      parameter (trseed=1.2)
      double precision avemldlight(1:NI,1:NJ)
      double precision light(1:NI,1:NJ,-1:NK)
      double precision parlight(0:150)
      data  parlight /& 
     &          0.0,             0.0,             0.0,             0.0,&
     &          0.0,             0.0,             0.0,             0.0,&
     &          0.0,             0.0,             0.0,             0.0,&
     &          0.0,             0.0,             0.0,             0.0,&
     &          0.0,   7.2833333,   7.8500000,   8.4166667,            &
     &   8.9833333,   9.5500000,   1.0116667d+01,   1.0683333d+01,     &
     &   1.1250000d+01,  1.1816667d+01,  1.2383333d+01,  1.2950000d+01,&
     &   1.3516667d+01,  1.4083333d+01,  1.4650000d+01,  1.5216667d+01,&
     &   1.5783333d+01,  1.6350000d+01,  1.6916667d+01,  1.7483333d+01,&
     &   1.8050000d+01,  1.8616667d+01,  1.9183333d+01,  1.9750000d+01,&
     &   2.0316667d+01,  2.0883333d+01,  2.1450000d+01,  2.2016667d+01,&
     &   2.2583333d+01,  2.3150000d+01,  2.3716667d+01,  2.4466667d+01,&
     &   2.5400000d+01,  2.6333333d+01,  2.7266667d+01,  2.8200000d+01,&
     &   2.9133333d+01,  3.0066667d+01,  3.1000000d+01,  3.1933333d+01,&
     &   3.2866667d+01,  3.3800000d+01,  3.4733333d+01,  3.5666667d+01,&
     &   3.6600000d+01,  3.7533333d+01,  3.8466667d+01,  3.9400000d+01,&
     &   4.0333333d+01,  4.1266667d+01,  4.2200000d+01,  4.3133333d+01,&
     &   4.4066667d+01,  4.5000000d+01,  4.5933333d+01,  4.6866667d+01,&
     &   4.7800000d+01,  4.8733333d+01,  4.9666667d+01,  5.0600000d+01,&
     &   5.1533333d+01,  5.2770492d+01,  5.4311475d+01,  5.5852459d+01,&
     &   5.7393443d+01,  5.8934426d+01,  6.0475410d+01,  6.2016393d+01,&
     &   6.3557377d+01,  6.5098361d+01,  6.6639344d+01,  6.8180328d+01,&
     &   6.9721311d+01,  7.1262295d+01,  7.2803279d+01,  7.4344262d+01,&
     &   7.5885246d+01,  7.7426230d+01,  7.8967213d+01,  8.0508197d+01,&
     &   8.2049180d+01,  8.3590164d+01,  8.5131148d+01,  8.6672131d+01,&
     &   8.8213115d+01,  8.9754098d+01,  9.1295082d+01,  9.2836066d+01,&
     &   9.4377049d+01,  9.5918033d+01,  9.7459016d+01,  9.9000000d+01,&
     &   1.0040984d+02,  1.0181967d+02,  1.0322951d+02,  1.0463934d+02,&
     &   1.0604918d+02,  1.0745902d+02,  1.0886885d+02,  1.1027869d+02,&
     &   1.1168852d+02,  1.1309836d+02,  1.1450820d+02,  1.1591803d+02,&
     &   1.1732787d+02,  1.1873770d+02,  1.2014754d+02,  1.2155738d+02,&
     &   1.2296721d+02,  1.2437705d+02,  1.2578689d+02,  1.2719672d+02,&
     &   1.2860656d+02,  1.3001639d+02,  1.3142623d+02,  1.3283607d+02,&
     &   1.3424590d+02,  1.3565574d+02,  1.3706557d+02,  1.3847541d+02,&
     &   1.3988525d+02,  1.4129508d+02,  1.4240984d+02,  1.4322951d+02,&
     &   1.4404918d+02,  1.4486885d+02,  1.4568852d+02,  1.4650820d+02,&
     &   1.4732787d+02,  1.4814754d+02,  1.4896721d+02,  1.4978689d+02,&
     &   1.5060656d+02,  1.5142623d+02,  1.5224590d+02 / 

! light0= PAR of 30 Ein/m2/day (from April Climatology -NAtl) is 30x1e6/86400 uEin/m2/s
! Muliply by 0.21 to get PAR=72 watts/m2
!     For newlight, inc light0 to 90.(too much) previously light0=32(too lit)
! mumax (per d) based on mu0 at 10deg
! alpha in m2/watt/day

      parameter (alpha=0.0538, mumax=0.536d0,chl2c=0.0667d0,wphy=1.2, &
     &     mort=0.0748,KzmaxTr=1.d-2, drhoml=0.01,swconvert=0.43, &
     &     KChl=0.041,Kwater=0.059)  

!old      parameter (alpha=0.046, mumax=0.72d0,chl2c=0.02d0,
!old    &     KzmaxTr=1.d-2, drhoml=0.01,swconvert=0.43)  


      double precision z(NK),growth(NK),flux(0:NK)
!     drhoml is the change in density from surf for defining MLD

      deltime=dtimel*TL
      timday= dble(step)*dtimel*TL/86400.d0
      DLinv= 1.d0/DL
      wsink= wphy/86400.d0       ! sinking vel 1.7 m/day

      death= mort/(86400.0)     ! loss rate from Katja -no grazing accounted
!      lightdepth= 1.d0/0.059    ! 16.9 m attenuation depth
!      chldepth= 1.d0/0.041      ! 25 m attenuation depth for Chl
      yearday= timday + yrday_start   ! simulation starts yearday 50
!      sw= 106.5 + 36.5*dtanh(40.d0*((yearday-103.d0)/103.d0)*pi/2.d0)
!     Sw is 70w/m^2 before yearday 103, and 143 w/m^2 after-transition in 5d
!      ydl1=103  !used until dec 8,2010
!      ydl2=113  !used until dec 8,2010
!      ydl1=100  
!      ydl2=120
!      if (yearday.le.ydl1) then
!         sw=70.d0
!      else if (yearday.ge.ydl2) then 
!         sw=143.d0
!      else 
!         sw = 143.d0*(yearday-ydl1)/(ydl2-ydl1)
!     &        + 70.d0*(ydl2-yearday)/(ydl2-ydl1)
!      end if 
!== Fixlight
!      sw= 70.d0
!      if (step.eq.1) then
!         open(unit=30,file='PAR_NOC.dat')
!         read(30,*) parlight
!         close(30)
!      end if
      int_yd= int(yearday)
      weight= yearday-dble(int_yd)
      sw= weight*parlight(int_yd+1) + (1.d0-weight)* parlight(int_yd)

! fix SW at 100      for null case
!      sw= 100.d0

      swpar=sw
      light0=swconvert*sw     

      do j=1,NJ
         do i=1,NI
!     Find MLD mld
            rhocrit =rho(i,j,NK)+ drhoml
            do k=NK-1,1,-1
               if ( rho(i,j,k).ge.rhocrit ) then
                  mld(i,j)= -DL*( (zc(i,j,k)-zc(i,j,k+1))*    &
     &            (rhocrit -rho(i,j,k+1))/(rho(i,j,k)-rho(i,j,k+1)) &
     &            + zc(i,j,k+1)  )
                  goto 51
               end if
            end do
!     if the do loop is computed, and the ML base is not found, set it to D
            mld(i,j)= -zf(i,j,0)*DL
 51         chlinteg=0.d0
!     Find light
            do k=NK,1,-1
               chl= Tr(1,i,j,k,0)*chl2c
               chlinteg= chlinteg + chl*(zf(i,j,k)-zf(i,j,k-1))*DL
               z(k)= zf(i,j,k)*DL  ! z is negative
               light(i,j,k)= light0*dexp( Kwater*z(k) &
     &              -chlinteg*KChl )
            end do
!     Find ML light
            lightmld = 0.d0
            dzsum= 0.d0
            do k=NK,1,-1
               if (-1.*z(k).le.mld(i,j)) then
                  dz= (zf(i,j,k)-zf(i,j,k-1))*DL
                  lightmld= lightmld+ light(i,j,k)*dz
                  dzsum= dzsum+ dz
               else
                  lightmld= lightmld/dzsum
                  avemldlight(i,j)= lightmld
                  goto 101
               end if
            end do

               
!     Average light in the mixed layer obtained by taking the integral/ mld
! 101        lightmld= light0 *(1.d0- dexp(-mld(i,j)/lightdepth))*
!     &           lightdepth / mld(i,j)

 101        do k=NK,1,-1

               z(k)= zf(i,j,k)*DL  ! z is negative
               if (-z(k).le.mld(i,j)) then    !ave light over MLD
                  light(i,j,k)=lightmld
               end if
               temp=alpha*alpha*light(i,j,k)*light(i,j,k)
               growth(k)= mumax*alpha*light(i,j,k) &
     &              /dsqrt(mumax*mumax + temp) /86400.d0
            end do

            flux(0)= 0.0
            flux(NK)= 0.0
            do k=NK-1,1,-1
!     Kz is at cell faces
               dpdz= (Tr(1,i,j,k+1,0)-Tr(1,i,j,k,0))*wz(i,j,k)*DLinv
               flux(k) = KzmaxTr*Kz(i,j,k)*dpdz + wsink*Tr(1,i,j,k+1,0)
            end do

            do k=1,NK
!     Since everything is in dim form, deltime is in seconds
!               if ((j.eq.NJ/2).and.(i.eq.NI/2)) then
!                  oldT(k)=T(1,i,j,k,0)  
!                  netprod(k)=(growth(k)-death)*deltime*T(1,i,j,k,0)  
!             fluxdiv(k)= deltime*(flux(k)-flux(k-1))*wz(i,j,k)*DLinv 
!               end if
               Tr(1,i,j,k,0)= Tr(1,i,j,k,0) + deltime*(       &
     &              (flux(k)-flux(k-1))*wz(i,j,k)*DLinv +   &
     &              (growth(k)-death)*Tr(1,i,j,k,0)  )
               consump(i,j,k,1)=(growth(k)-death)*Tr(1,i,j,k,0)  
! consump is in mg C per m^3 per second
            end do
         end do
      end do

      fac=deltime/86400.d0    ! restore on time scale of 1 day
      do j=1,NJ
         do i=1,NI
            Tave= 0.d0
            dzsum= 0.d0
            do k=NK,1,-1
               zdep= zf(i,j,k)*DL ! z is negative
               if (-zdep.le.mld(i,j) ) then 
                  dz= (zf(i,j,k+1)-zf(i,j,k))*DL
                  Tave= Tr(1,i,j,k,0)*dz+ Tave
                  dzsum= dzsum+ dz
               end if
            end do
            Tave= Tave/dzsum
            if (Tave.lt.trseed) then
               do k=NK,1,-1
                  zdep= zf(i,j,k)*DL ! z is negative
                  if (-zdep.le.mld(i,j)) then 
                     Tr(1,i,j,k,0)= Tr(1,i,j,k,0)+ (trseed-Tave)*fac
                  end if
               end do
            end if
         end do
      end do

!
!=      call phymix
!
      do m=0,1
         do k=0,NK+1
            do j=0,NJ+1
               do it=1,ntr
                  Tr(it,0,j,k,m)= Tr(it,NI,j,k,m)
                  Tr(it,NI+1,j,k,m)= Tr(it,1,j,k,m)
               end do
            end do
         end do
      end do
!
      return

      end


      subroutine phymix
!     Homogenizes Tr(0,...) in the ML
      USE header
      implicit none
      integer i,j,k
      double precision physum,dzsum,dz,z(NK)

!
!     Mix PHY
      do j=1,NJ
         do i=1,NI
            physum= 0.d0
            dzsum= 0.d0
            
            do k=NK,1,-1

               z(k)= -zf(i,j,k)*DL ! z is negative
               dz= zf(i,j,k)-zf(i,j,k-1)
               if (z(k).le.mld(i,j)) then    !mixover MLD
                  physum=physum+ Tr(1,i,j,k,0)*dz
                  dzsum= dz+ dzsum
               else
                  goto 301
               end if
            end do
 301        physum= physum/dzsum
            do k=NK,1,-1
               if (z(k).lt.mld(i,j)) then
                 Tr(1,i,j,k,0)= physum
               else
                  goto 302
               end if
            end do
 302        continue
         end do
      end do

      return
      end







