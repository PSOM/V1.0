subroutine stprofile 
  !     --------------------                                              
  USE header
!     initializes s as pot.density with a single vertical profile       
!     Across-front use the TANH profile with tightness as a measure of f
!     spread                                                            
!     Larger the factor, tighter the front and larger b_xx.             
!     tightness=4, gives the needed b_xx to make b_xx H/ (2f^2) < 1.    
      IMPLICIT NONE 
      integer  i,j,k,iseed,nplm1,iter,iter2 
      integer, parameter :: npl=33
      REAL(kind=rc_kind) :: rhovert(npl),Tvert(npl),svert(npl),dep(npl),         &
           depoff(npl) 
      REAL(kind=rc_kind) :: bs1(npl),cs1(npl),ds1(npl),bT1(npl),cT1(npl),dT1(npl), &
     &     z,seval,zmm,sbkgrnd,Tbkgrnd,z1,z2,zprev,alpha,dtdz,dsdz
      REAL(kind=rc_kind) :: slfac,dscl,rdscl,yarg,ex2y,thy,ran3,         &
     &     perturb,slfacnew,dz,bfsqbkgrnd,wiggles,amplitude             
!     tightness = 10 represents a very tight front, =1 loose front      
!     mldepth is the desired ML depth                                   
!= usual      parameter (mldepth= 0.d0, tightness=0.03) ! dy=1km        
!      parameter (mldepth= 100.d0 ) !tightness below ,used in inith_fixd
!      parameter (mldepth= 200.d0 ) !tightness below ,used in inith_fixd
      parameter (alpha=1.6d-4)   ! thermal expansion coeff
!     Average rho profile for 40N from DATA/ts/Levitus_Atlas.nc         
      data dep / -5500.d0, -5000.d0, -4500.d0, -4000.d0, -3500.d0,      &
     &      -3000.d0, -2500.d0, -2000.d0, -1750.d0, -1500.d0, -1400.d0, &
     &      -1300.d0, -1200.d0, -1100.d0, -1000.d0, -900.d0, -800.d0,   &
     &      -700.d0, -600.d0, -500.d0, -400.d0, -300.d0, -250.d0,       &
     &      -200.d0, -150.d0, -125.d0, -100.d0, -75.d0, -50.d0,         &
     &      -30.d0, -20.d0, -10.d0, -0.d0 /                             
!      data rhovert /  1027.7694, 1027.7852, 1027.7894, 1027.7886,       &
!     &     1027.7897, 1027.8050, 1027.8165, 1027.7848, 1027.7630,       &
!     &     1027.7297, 1027.7115, 1027.6846, 1027.6727, 1027.6372,       &
!     &     1027.5988, 1027.5540, 1027.5026,                             &
!     &     1027.4337, 1027.3380, 1027.2524, 1027.1506, 1027.0222,       &
!     &     1026.9653, 1026.8773, 1026.7458, 1026.6521,                  &
!     &     1026.5293, 1026.3756, 1026.0901, 1025.7001, 1025.5029,       &
!     &     1025.3558, 1025.3215 /                                       
      data svert / 34.76,  34.87, 34.91, 34.94,  35.00,  35.05, 35.08,  &
           35.11,  35.13,  35.13,  35.11,  35.05,  34.99,  34.94,       &
           34.89,  34.91,  34.89,  34.89,  34.88,  34.89,  34.90,       &
           34.88,  34.88,  34.88,  34.88,  34.88,  34.87,  34.81,       &
           34.76,  34.75,  34.75,  34.74,  34.71 /                            
      data tvert / 1.61, 1.72, 1.77, 1.78, 1.87, 2.18, 2.61, 3.05,      &
           3.28, 3.62, 3.80, 4.06, 4.32, 4.57, 4.84, 5.29, 5.71,        &
           6.37, 6.96, 7.83, 8.75, 9.82, 10.42, 11.00, 11.71, 12.12,    &
           12.63, 13.28, 14.46, 16.01, 16.76, 17.04, 17.25/

!     -----------------                                                 
!     Specify MLD                                                       
      mldepth= 50.d0 
!      mldepth= 200.d0                                                  
                                                                        
!      tightness= 0.3d0   !used in inith_fixdh                          
                          !for larger domain                            
      tightness= 0.03d0 
      bfsqbkgrnd = 1.d-6 
      dtdz= bfsqbkgrnd/(gpr*10.d0)/alpha
      dsdz = 0.d0

      do k=1,npl 
         depoff(k)= dep(k)- mldepth 
      end do 
                                                                        
      call spline (npl,depoff,svert,bs1,cs1,ds1) 
      call spline (npl,depoff,tvert,bT1,cT1,dT1) 
                                                                        
!     sclwidth is the scale width of the channel (i.e. orig width = 48km
!     yfront is the distance of the front from the coast                
      sclwidth = 48.0 
      yfront = 0.5*(yc(NJ+1) +yc(0)) 
!-offset      yfront = (yc(NJ+1) +yc(0))/3.d0  !used for many MLI runs  
                                                                        
!     ADD a PERTURBATION to the width of the front                      
      iseed= 44294 
      dum = ran3(iseed) 
                                                                        
!     z1 is the depth of the ml (diff rho on both sides above this depth
!     z2 is the vertical extent of the density anamoly (it is gradually 
!        linearly anihillated with depth).                              
!     Orig vals. z1= 50. z2=250.                                        
      z1= mldepth - 10.d0 
!      z2= mldepth + 10.d0 
      z2= mldepth + 200.d0 
                      !drho=0.2     ! for 100-200 m deep ML (using conve
!      slfac= 0.
      slfac= -(0.3/1025.d0)/alpha 
!      slfac= 0.15d0                ! for 100 m deep ML (using convect) 
                                                                        
!     0.12 in pot dens, 0.15 in salinity                                
      do j=0,NJ+1 
         do i=0,NI+1 
            do k=0,NK+1 
               z= DL*zc(i,j,k) 
                  if (z.ge.-mldepth) then 
                     Tbkgrnd =                                          &
     &               seval(npl,-1.*mldepth,depoff,tvert,bT1,cT1,dT1)     &
     &                       - (z+mldepth)*dtdz
                     sbkgrnd =                                          &
     &               seval(npl,-1.*mldepth,depoff,svert,bs1,cs1,ds1)     &
     &                       - (z+mldepth)*dsdz
!     &               seval(npl,-1.*mldepth,depoff,svert,bs1,cs1,ds1)     &
!     &                       - (z+mldepth)*bfsqbkgrnd*R0/(gpr*10.)      
                  else 
                     Tbkgrnd =                                          &
     &                    seval(npl,z,depoff,tvert,bT1,cT1,dT1)        
                     sbkgrnd =                                          &
     &                    seval(npl,z,depoff,svert,bs1,cs1,ds1)        
!     &                    seval(npl,z,depoff,svert,bs1,cs1,ds1)        
                  end if 
!==                                                                     
                  T(i,j,k,0)=  Tbkgrnd 
                  s(i,j,k,0)=  sbkgrnd 

!                  if (k.eq.NK) then                                    
!                     slfacnew= slfac*2.d0                              
!                     if ((-1.d0*slfacnew*thy).gt.0.d0)                 
!     &                    s(i,j,k,0)= -slfacnew*thy + s(i,j,k,0)       
!                  end if                                               
!                                                                       
!-               else                                                   
!-                  s(i,j,k,0)=                                         
!-     &                 seval(npl,z,dep,svert,bs1,cs1,ds1) -S0          
!-               end if                                                 
            end do 
         end do 
      end do 
 
                                                                        
!      do iter2=1,40                                                    
!         call conadjust(0)                                             
!      end do                                                           
                                                                        
!     WIGGLE in FRONT                                                   
      wiggles=1.d0 
      amplitude= 2.0d0 
      do j=0,NJ+1 
                                                                        
         do i=0,NI+1 
            yfront= 0.5d0*(yc(NJ+1) +yc(0)) + amplitude*dsin(2.d0*PI*xc(i)/(0.5*(xc(NI+1)+xc(NI))/wiggles)  )                     
            thy = tanh(tightness*(yc(j)-yfront)*PI) 
  
            do k=1,NK 
               z= DL*zc(i,j,k) 
!c=               if (z.gt.-mldepth) then 
                                                                        
                  if (z.ge.-z1) then 
                     slfacnew = slfac 
                  else if (z.ge.-z2) then 
                     slfacnew = slfac*(z+z2)/(z2-z1) 
                  else 
                     slfacnew = 0.d0 
                  end if 
!                  if ((i.eq.1).and.(j.eq.1)) write(6,*) 'slfc',k,slfacnew
                                                                        
                  T(i,j,k,0)= slfacnew*(thy-1.d0) + T(i,NJ,k,0) 
                                                                        
!     if (z.ge.-z1) then                                                
!                  sbkgrnd =                                            
!     &              seval(npl,-1.*z1,depoff,rhovert,bs1,cs1,ds1)        
!     &              - (z+z1)*bfsqbkgrnd*R0/(gpr*10.)                   
!               s(i,j,k,0)= slfac*(1.d0+thy) + sbkgrnd                  
!     &               seval(npl,-1.*mldepth,depoff,rhovert,bs1,cs1,ds1)  
!     &                    -(z+mldepth)*bfsqbkgrnd*R0/(gpr*10.)         
                                                                        
!c=               end if 
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
            T(i,0,k,0)= T(i,1,k,0) 
            T(i,NJ+1,k,0)= T(i,NJ,k,0) 
         end do 
!     periodicew                                                        
         do j=0,NJ+1 
            s(0,j,k,0)= s(NI,j,k,0) 
            s(NI+1,j,k,0)= s(1,j,k,0) 
            T(0,j,k,0)= T(NI,j,k,0) 
            T(NI+1,j,k,0)= T(1,j,k,0) 
         end do 
      end do 
                                                                        
                                                                        
      return 
      END                                           
