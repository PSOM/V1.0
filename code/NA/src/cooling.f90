      subroutine cooling(rhoflux,yrday)
!     ---------------------------------------------
!     heat flux to alter rho at k=NK
!     
      USE header
      implicit none
      integer i,j,k,iseed,num,nwork
      integer n1,n2
!     qflux declared in header as common
!      double precision alpha,Cp,rhoflux(NI,NJ),qfluxday(0:120) increase by 30 days
      double precision alpha,Cp,rhoflux(NI,NJ),qfluxday(0:150)
!      parameter (alpha=2.5d-4, Cp= 4127.d0)   ! old qflux=50.d0)
!      parameter (alpha=1.d-4, Cp= 3988.d0)   ! 
      parameter (alpha=1.6d-4, Cp= 3988.d0)   ! 
!     qflux and alpha are both negative, and so sign is neglected
!     random number
      parameter (nwork=100000)
      double precision dtemp,wflux,Tflux,rflux,fun_potdens,delflux_nor, &
     &     delflux_sou
      double precision da(NI,NJ),dwork(nwork),am,sd,yrday
      double precision Tflux_south,rflux_south,wts,                 &
     &     Tflux_north,rflux_north
!     NOC Heat fluxes, from yr day 30 to 150., constant before
      data qfluxday / &          ! yrday 0-20 - wrong vals. not using this portion
     &        0.0,             0.0,             0.0,             0.0,&
     &        0.0,             0.0,             0.0,             0.0,&
     &        0.0,             0.0,             0.0,             0.0,&
     &        0.0,             0.0,             0.0,             0.0,&
     &        0.0, -2.6301667d+02, -2.6105000d+02, -2.5908333d+02,   &
     &  -2.5711667d+02, -2.5515000d+02, -2.5318333d+02, -2.5121667d+02,&
     &  -2.4925000d+02, -2.4728333d+02, -2.4531667d+02, -2.4335000d+02,&
     &  -2.4138333d+02, -2.3941667d+02, -2.3745000d+02, -2.3548333d+02,&
     &  -2.3351667d+02, -2.3155000d+02, -2.2958333d+02, -2.2761667d+02,&
     &  -2.2565000d+02, -2.2368333d+02, -2.2171667d+02, -2.1975000d+02,&
     &  -2.1778333d+02, -2.1581667d+02, -2.1385000d+02, -2.1188333d+02,&
     &  -2.0991667d+02, -2.0795000d+02, -2.0598333d+02, -2.0480000d+02,&
     &  -2.0440000d+02, -2.0400000d+02, -2.0360000d+02, -2.0320000d+02,&
     &  -2.0280000d+02, -2.0240000d+02, -2.0200000d+02, -2.0160000d+02,&
     &  -2.0120000d+02, -2.0080000d+02, -2.0040000d+02, -2.0000000d+02,&
     &  -1.9960000d+02, -1.9920000d+02, -1.9880000d+02, -1.9840000d+02,&
     &  -1.9800000d+02, -1.9760000d+02, -1.9720000d+02, -1.9680000d+02,&
     &  -1.9640000d+02, -1.9600000d+02, -1.9560000d+02, -1.9520000d+02,&
     &  -1.9480000d+02, -1.9440000d+02, -1.9400000d+02, -1.9360000d+02,&
     &  -1.9320000d+02, -1.9096721d+02, -1.8690164d+02, -1.8283607d+02,&
     &  -1.7877049d+02, -1.7470492d+02, -1.7063934d+02, -1.6657377d+02,&
     &  -1.6250820d+02, -1.5844262d+02, -1.5437705d+02, -1.5031148d+02,&
     &  -1.4624590d+02, -1.4218033d+02, -1.3811475d+02, -1.3404918d+02,&
     &  -1.2998361d+02, -1.2591803d+02, -1.2185246d+02, -1.1778689d+02,&
     &  -1.1372131d+02, -1.0965574d+02, -1.0559016d+02, -1.0152459d+02,&
     &  -9.7459016d+01, -9.3393443d+01, -8.9327869d+01, -8.5262295d+01,&
     &  -8.1196721d+01, -7.7131148d+01, -7.3065574d+01, -6.9000000d+01,&
     &  -6.5688525d+01, -6.2377049d+01, -5.9065574d+01, -5.5754098d+01,&
     &  -5.2442623d+01, -4.9131148d+01, -4.5819672d+01, -4.2508197d+01,&
     &  -3.9196721d+01, -3.5885246d+01, -3.2573770d+01, -2.9262295d+01,&
     &  -2.5950820d+01, -2.2639344d+01, -1.9327869d+01, -1.6016393d+01,&
     &  -1.2704918d+01, -9.3934426d+00, -6.0819672d+00, -2.7704918d+00,&
     &   5.4098361d-01,  3.8524590d+00,  7.1639344d+00,  1.0475410d+01,&
     &   1.3786885d+01,  1.7098361d+01,  2.0409836d+01,  2.3721311d+01,&
     &   2.7032787d+01,  3.0344262d+01,  3.2852459d+01,  3.4557377d+01,&
     &   3.6262295d+01,  3.7967213d+01,  3.9672131d+01,  4.1377049d+01,&  
     &   4.3081967d+01,  4.4786885d+01,  4.6491803d+01,  4.8196721d+01,&
     &   4.9901639d+01,  5.1606557d+01,  5.3311475d+01/

!     qflux is the cooling flux (negative heat flux)
!     NOC heat flux
      n1=yrday 
      n2=n1+1
      qflux=(yrday-dble(n1))*qfluxday(n2)+ (dble(n2)-yrday)* &
     &     qfluxday(n1) 

      qflux =qflux + 40.d0

!     constant cooling
!      qflux= -100.d0
!     ==========================

      delflux_nor= 25.d0 ! this amount is subtracted at N.bndry
      delflux_sou= 25.d0 ! this amount is  added  at S.bndry

!     Zero heat flux
!     ===============
!      qflux= 0.d0
!      qflux= -150.d0            !constant cooling
!      qflux=dmin1(qflux,-150.d0) !max hf= -150w/m2
!     ===== for no fronts, and no flux
!      delflux_nor= 0.d0
!      delflux_sou= 0.d0
!     ========

      iseed= 44294
      num=NI*NJ
      am=  0.d0
      sd= 0.01d0
      call rannw (am, sd, iseed, da, num, dwork, nwork)
      call rannw (am, sd, iseed, da, num, dwork, nwork)
!     This generates a ran number array with mean 0 and sd=.01
!      do j=1,NJ
!         do i=1,NI
!            write(60,*) da(i,j)
!         end do
!      end do
!      write(6,*) 'stop after random number'
!      stop
!     da (random - white noise) adds a 1% variation to qflux

!     rhoflux is in (kg/m^2) per second, i.e. (kg/m3)*(m/s)

      Tflux=  qflux/(R0*Cp)
!      wflux= (dztop*DL)/(dtf*TL)
!      dtemp= Tflux/wflux
!      rflux=(fun_potdens(s0,T0+dtemp)-fun_potdens(S0,T0))*wflux
      rflux= -Tflux*alpha*R0

      Tflux_south=  (qflux+ delflux_sou)/(R0*Cp)
      rflux_south= -Tflux_south*alpha*R0      

      Tflux_north=  (qflux- delflux_nor)/(R0*Cp)
      rflux_north= -Tflux_north*alpha*R0      

      do j=1,NJ
         wts=dble(NJ-j)/dble(NJ) 
         do i=1,NI
!            rhoflux(i,j)=(1.d0+da(i,j)*1.d-2)*rflux ! add a random perturb
!           Add latitudinal variation
            rhoflux(i,j)=((1.d0-wts)*rflux_north + wts*rflux_south)* &
     &           (1.d0+da(i,j)*1.d-2)
         end do
      end do

      return
      end


