      subroutine ini_st
!     --------------------
!     initializes s as pot.density with a single vertical profile
!     Across-front use the TANH profile with tightness as a measure of frontal
!     spread
!     Larger the factor, tighter the front and larger b_xx.
!     tightness=4, gives the needed b_xx to make b_xx H/ (2f^2) < 1.
      USE header
      implicit none
      integer  i,j,k,iseed,npm1,iter,iter2,npl
      parameter (npl=1001)
      double precision rhovert(npl),Tvert(npl),dep(npl),depoff(npl)
      double precision bs1(npl),cs1(npl),ds1(npl),bT1(npl),cT1(npl),dT1(npl), &
     &     z,seval,zmm,sbkgrnd,z1,z2,zprev,yfront1,yfront2, &
     &     yfront3
      double precision slfac,dscl,rdscl,yarg,ex2y,thy,ran3, &
     &     perturb,slfacnew,dz,bfsqbkgrnd,wiggles,amplitude


! ADDITION Dec. 19
     !double precision  drhofront,njf1, njf2, njf3, njf4
! CORRECTION 27Mar2014
     double precision  drhofront
     integer njf1, njf2, njf3, njf4



!     tightness = 10 represents a very tight front, =1 loose front
!     mldepth is the desired ML depth
!= usual      parameter (mldepth= 0.d0, tightness=0.03) ! dy=1km
!      parameter (mldepth= 100.d0 ) !tightness below ,used in inith_fixdh 
!      parameter (mldepth= 200.d0 ) !tightness below ,used in inith_fixdh 
!     Average rho profile for 40N from DATA/ts/Levitus_Atlas.nc
!      data dep / -5500.d0, -5000.d0, -4500.d0, -4000.d0, -3500.d0,
!     &      -3000.d0, -2500.d0, -2000.d0, -1750.d0, -1500.d0, -1400.d0,
!     &      -1300.d0, -1200.d0, -1100.d0, -1000.d0, -900.d0, -800.d0,
!     &      -700.d0, -600.d0, -500.d0, -400.d0, -300.d0, -250.d0,
!     &      -200.d0, -150.d0, -125.d0, -100.d0, -75.d0, -50.d0,
!     &      -30.d0, -20.d0, -10.d0, -0.d0 / 
!      data rhovert /  1027.7694, 1027.7852, 1027.7894, 1027.7886, 
!     &     1027.7897, 1027.8050, 1027.8165, 1027.7848, 1027.7630,
!     &     1027.7297, 1027.7115, 1027.6846, 1027.6727, 1027.6372,
!     &     1027.5988, 1027.5540, 1027.5026,
!     &     1027.4337, 1027.3380, 1027.2524, 1027.1506, 1027.0222,
!     &     1026.9653, 1026.8773, 1026.7458, 1026.6521, 
!     &     1026.5293, 1026.3756, 1026.0901, 1025.7001, 1025.5029,
!     &     1025.3558, 1025.3215 / 
!
!      data dep /-1100.,-800.,-600.,-450.,-350.,-305.,-275.,
!     &    -250.,-240.,-230.,-220.,-210.,-200.,-190.,-180.,
!     &    -170.,-160.,-150.,-140.,-130.,-120.,-110.,-100.,
!     &    -90.,-80.,-70.,-60.,-50.,-40.,-30.,-20.,-10., 0./
!         
!     double n2  (temporarily)
!==   ---------------------
!      do k=np-1,1,-1
!         drho(k)= rhovert(k+1)-rhovert(k)
!      end do
!      do k=np-1,1,-1
!         rhovert(k)= rhovert(k+1)-2.d0*drho(k)
!      end do
!     ----------------
!     Read in density profile
      open(unit=31,file='NA/NABprofile-200m.in')
      read(31,*) dep,rhovert
      close(31)
      do k=1,npl
         rhovert(k)=rhovert(k)+1000.d0   
      end do

!     -----------------
!     Increase density by 2.1 ! not needed with NAB data
!      do i=1,np
!         rhovert(i)= rhovert(i)+2.1
!      end do
!     Specify MLD
!      mldepth= 200.d0
!      mldepth= 300.d0
      mldepth= 250.d0
!used with nj640      mldepth= 250.d0   - 
      z1= mldepth
!      z2= mldepth+100.d0
!!!      z2= 500.d0  ORIG VALUE USED IN NAB08 and Science paper
      z2= 800.d0  ! increased to get larger eddies for C export
!      z2= mldepth+50.d0

!      tightness= 0.3d0   !used in inith_fixdh
!      used throughout      tightness= 0.03d0   !for larger domain
      tightness= 0.03d0   ! by=.9e-7 drho=0.2 
!      tightness = 0.02d0    ! by=0.3x10^-7 for  dho=0.1
!=used      tightness = 0.01d0    ! wide front
!      bfsqbkgrnd = 1.d-6  !reduce further
      bfsqbkgrnd = 0.3d-6   ! not used any more
!      write(6,*) 'ML background stratific N2 = ', bfsqbkgrnd
!     no longer offset - just fit to winter profile 
!      do k=1,np
!         depoff(k)= dep(k)- mldepth
!      end do
      
!      call spline (np,depoff,rhovert,bs1,cs1,ds1)
      call spline (npl,dep,rhovert,bs1,cs1,ds1)

!     sclwidth is the scale width of the channel (i.e. orig width = 48km)
!     yfront is the distance of the front from the coast
      sclwidth = 48.0
!      yfront = 0.5*(yc(NJ+1) +yc(0))
      yfront2 = 0.5*(yc(NJ-1) +yc(0))
      yfront1 = 0.5*yfront2
      yfront3 = 1.5*yfront2
!-offset      yfront = (yc(NJ+1) +yc(0))/3.d0  !used for many MLI runs

!     ADD a PERTURBATION to the width of the front
      iseed= 44294
      dum = ran3(iseed)

!     z1 is the depth of the ml (diff rho on both sides above this depth)
!     z2 is the vertical extent of the density anamoly (it is gradually 
!        linearly anihillated with depth).
!     Orig vals. z1= 50. z2=250.
!nolonger      z1= mldepth - 10.d0
!      z2= mldepth + 10.d0
!      slfac= -0.15d0
!      slfac= 0.1                ! for 100 m deep ML (using convect)
      slfac= 0.04   !3fronts  ! for 200 m deep ML (NAB), 300m MLD, y=480 km 
!      slfac= 0.024  !5fronts  ! for 200 m deep ML (NAB), 300m MLD, y=480 km 
!===
!      slfac= 0.d0    ! NO FRONTS
!--

!      slfac= 0.033               ! for 300 m deep ML (using convect)
      drhofront=slfac

!      slfac= 0.025               ! for 100 m deep ML (using convect)

!     0.12 in pot dens, 0.15 in salinity
      do j=0,NJ+1
         do i=0,NI+1
            do k=0,NK+1
               z= DL*zc(i,j,k)
!                  if (z.ge.-mldepth) then
!                     sbkgrnd = 
!     &               seval(np,-1.*mldepth,depoff,rhovert,bs1,cs1,ds1)
!     &                       - (z+mldepth)*bfsqbkgrnd*R0/(gpr*10.)
!                  else
!                     sbkgrnd = 
!     &                    seval(np,z,depoff,rhovert,bs1,cs1,ds1) 
!                  end if
               z=-z
                  sbkgrnd =  seval(npl,z,dep,rhovert,bs1,cs1,ds1) 

!==
                  s(i,j,k,0)=  sbkgrnd 
!                  if (k.eq.NK) then 
!                     slfacnew= slfac*2.d0
!                     if ((-1.d0*slfacnew*thy).gt.0.d0) 
!     &                    s(i,j,k,0)= -slfacnew*thy + s(i,j,k,0)
!                  end if
!
!-               else
!-                  s(i,j,k,0)= 
!-     &                 seval(np,z,dep,svert,bs1,cs1,ds1) -S0
!-               end if
            end do
         end do
      end do


!      do iter2=1,40
!         call conadjust(0)
!      end do


!      njf1=NJ/3
!      njf2=2*NJ/3
!      njf1=2*NJ/5
!      njf2=3*NJ/5
      njf1=0.375*NJ
      njf2=0.625*NJ
      njf3=njf1
      njf4=njf2
      write(6,*) 'njf1 , njf2 = ', njf1, njf2,njf3,njf4

!      wiggles=10.d0
      wiggles=2.d0
      amplitude= 0.1d0

!     FRONT 3
!     ---------
      do j=njf2+1,NJ+1
!     WIGGLE in FRONT
!         yfront= 0.5d0*(yc(NJ+1) +yc(njf2)) 
!         yfront= 0.5d0*(yc(NJ+1) +yc(njf2)) 
         yfront= yfront3 !+ amplitude*dsin(2.d0*PI*xc(i)/(0.5*(xc(NI+1)+xc(NI))/wiggles)  )
         thy = tanh(tightness*(yc(j)-yfront)*PI)
         do i=0,NI+1
            do k=1,NK
               z= DL*zc(i,j,k)
!               if (z.gt.-mldepth) then
               if (z.gt.-z2) then
                  
                  if (z.ge.-z1) then 
                     slfacnew = slfac
                  else if (z.ge.-z2) then
                     slfacnew = slfac*(z+z2)/(z2-z1)
                  else
                     slfacnew = 0.d0
                  end if

                  s(i,j,k,0)= slfacnew*(thy-1.d0) + s(i,NJ,k,0)
               end if
            end do
         end do
      end do
!----------
!     FRONT 2
!     ---------
      do j=njf1+1,njf2
!     WIGGLE in FRONT
!         yfront= 0.5d0*(yc(njf2) +yc(njf1)) 
         yfront= yfront2 !+ amplitude*dsin(2.d0*PI*xc(i)/(0.5*(xc(NI+1)+xc(NI))/wiggles)  )
         thy = tanh(tightness*(yc(j)-yfront)*PI)
         do i=0,NI+1
            do k=1,NK
               z= DL*zc(i,j,k)
!               if (z.gt.-mldepth) then
               if (z.gt.-z2) then
                  
                  if (z.ge.-z1) then 
                     slfacnew = slfac
                  else if (z.ge.-z2) then
                     slfacnew = slfac*(z+z2)/(z2-z1)
                  else
                     slfacnew = 0.d0
                  end if

                  s(i,j,k,0)= slfacnew*(thy-1.d0) + s(i,njf2+1,k,0)
               end if
            end do
         end do
      end do
!----------
!     FRONT 1
!     ---------
      do j=0,njf1
!     WIGGLE in FRONT
!         yfront= 0.5d0*(yc(1) +yc(njf1)) 
         yfront= yfront1 !+ amplitude*dsin(2.d0*PI*xc(i)/(0.5*(xc(NI+1)+xc(NI))/wiggles)  )
         thy = tanh(tightness*(yc(j)-yfront)*PI)
         do i=0,NI+1
            do k=1,NK
               z= DL*zc(i,j,k)
!               if (z.gt.-mldepth) then
               if (z.gt.-z2) then
                  
                  if (z.ge.-z1) then 
                     slfacnew = slfac
                  else if (z.ge.-z2) then
                     slfacnew = slfac*(z+z2)/(z2-z1)
                  else
                     slfacnew = 0.d0
                  end if

                  s(i,j,k,0)= slfacnew*(thy-1.d0) + s(i,njf1+1,k,0)
               end if
            end do
         end do
      end do

!      write(6,*) 'rho i=24,k=24', (s(24,j,24,0),j=1,NJ)
!      stop

!         do iter2=1,10
!            call conadjust(0)
!         end do

      do k=0,NK+1
         do i=1,NI
            s(i,0,k,0)= s(i,1,k,0)
            s(i,NJ+1,k,0)= s(i,NJ,k,0)
         end do
!     periodicew
         do j=0,NJ+1
            s(0,j,k,0)= s(NI,j,k,0)
            s(NI+1,j,k,0)= s(1,j,k,0)
         end do
      end do


      return
      end

