subroutine staticsigma 

#include "cppdefs.h"
USE header
! ------------------------------------------------------------------
!     FOR THE REGION BELOW THE TOPMOST LAYER                            
!     modified a for periodicew bc                                      
!     This subroutine updates all those quantities that are a function  
!     of time and  (time and depth). It is called  every time step.     
!     The surface bc  wfbct is updated here.                            
! ------------------------------------------------------------------
   implicit none
   integer i,j,k
   REAL(kind=rc_kind) :: zpd,hpd,hpdinv,hu,hv,hx,hy,hxpdx,hypdy,z,temp, &
      be2,wxk,wyk,d2,sig,ep,epm1,dnkm1,dnf,Df,              &
      g13(0:NI+1,0:NJ+1,0:NK+1),g23(0:NI+1,0:NJ+1,0:NK+1)
                          
   be2= beta*EPS*EPS 
   dnkm1= dble(NK-1)
   dnf= Nf
   Df= (-depmean_dim+distance)/DL
   
#ifdef sigma_stretch
   epm1= exp(pfac)-1.d0 
#endif

   ! ------------------------------------------------------------------
   !				        Derivatives of sigma (Upper Layers)
   ! ------------------------------------------------------------------
 
      do j=0,NJ+1
         do i=0,NI+1

            ! All these variables are functions of time
            hpd= dztop +D(i,j)
            hpdinv= 1.d0/hpd
         
            ! -------------------------------

             do 20 k=Nf+1,NK-1 ! upper layer
                sig= dble(k)-0.5d0
                z= zc(i,j,k)
                zpd= z +dztop
#ifdef sigma_stretch
                wz(i,j,k)= -epm1*(dnkm1-dnf)/pfac/( epm1*zpd +(Df+dztop) ) !!!
#else
                wz(i,j,k)= (dnkm1-dnf)/(-dztop-Df) !!!
#endif
                Jac(i,j,k)= J2d(i,j)/wz(i,j,k)

                ! wx(i,j,k)= wz(i,j,k)*(z+dztop)/(Df+dztop)*Ddx(i,j) ! evenly between -dztop and Df
                ! wy(i,j,k)= wz(i,j,k)*(z+dztop)/(Df+dztop)*Ddy(i,j) ! evenly between -dztop and Df
                wx(i,j,k)= 0.d0
                wy(i,j,k)= 0.d0
                g13(i,j,k)= ux(i,j)*wx(i,j,k) +uy(i,j)*wy(i,j,k)
                g23(i,j,k)= vx(i,j)*wx(i,j,k) +vy(i,j)*wy(i,j,k)
       20    continue

            ! ----------------------- Z FACE  -----------------------

            do k=Nf+1,NK-1   ! upper layer
               sig= dble(k)
               z= zc(i,j,k)  !!! use z at center
               zpd= z +dztop
               wt(i,j,k)= 0.d0
#ifdef sigma_stretch
               wzk(i,j,k)= -epm1*(dnkm1-dnf)/pfac/( epm1*zpd+(Df+dztop) ) !!!
#else
               wzk(i,j,k)= (dnkm1-dnf)/(-dztop-Df) !!!
#endif
            end do ! end k

            wzk(i,j,NK-1)= 0.5*( wz(i,j,NK) + wz(i,j,NK-1) )   ! Why take the average ???

            ! -----------------------  X, Y FACE  -----------------------
            do 25 k=Nf+2,NK  ! upper layer
               sig= dble(k)
               z= zf(i,j,k)
               zpd= z +dztop

   			   ! temp= epm1*(dnkm1-dnf)/pfac/( epm1*zpd +(Df+dztop) ) *(z+dztop)/(Df+dztop)
        	      ! wxk= Ddx(i,j)*temp
               ! wyk= Ddy(i,j)*temp

               wxk= 0.d0
               wyk= 0.d0
               gqk(i,j,k,1)= qpr*Jac(i,j,k)*(ux(i,j)*wxk +uy(i,j)*wyk)
               gqk(i,j,k,2)= qpr*Jac(i,j,k)*(vx(i,j)*wxk +vy(i,j)*wyk)
               gqk(i,j,k,3)= Jac(i,j,k)*(qpr*(wxk*wxk +wyk*wyk) +       &
      &           be2*wz(i,j,k)*wz(i,j,k))
      25    continue

      enddo ! end i
      enddo ! end j

! ------------------------------------------------------------------
!				        Derivatives of sigma (Lower Layers)
! ------------------------------------------------------------------
#ifdef fixed_bottom_thickness
      dnf= dble(Nf-1)    !!! ADDED
#else
      dzbot = 0.d0
#endif

! ------------------------------------------------------------------ 
   do j=0,NJ+1 
      do i=0,NI+1 

         ! All these variables are functions of time                         
         hpd= dztop +D(i,j) 
         hpdinv= 1.d0/hpd 

         ! ---------------------  CENTER  -----------------------

         do 10 k=0,Nf    ! lower layer
            sig= dble(k)-0.5d0
            z= zc(i,j,k) 
            zpd= z +dztop 
#ifdef sigma_stretch
            wz(i,j,k)= -epm1*dnf/pfac/( epm1*(z-Df)+(D(i,j)+dzbot-Df) ) !!!
#else
            wz(i,j,k)= dnf/(Df-D(i,j)-dzbot) !!!
#endif
            Jac(i,j,k)= J2d(i,j)/wz(i,j,k) 

            wx(i,j,k)= wz(i,j,k)*(z-Df)/(D(i,j)+dzbot-Df)*Ddx(i,j) ! evenly between Df and D
            wy(i,j,k)= wz(i,j,k)*(z-Df)/(D(i,j)+dzbot-Df)*Ddy(i,j) ! evenly between Df and D
            g13(i,j,k)= ux(i,j)*wx(i,j,k) +uy(i,j)*wy(i,j,k) 
            g23(i,j,k)= vx(i,j)*wx(i,j,k) +vy(i,j)*wy(i,j,k) 
   10    continue

         ! ----------------------- Z FACE  -----------------------

         do k=0,Nf    ! lower layer
            sig= dble(k)
            z= zc(i,j,k)  !!! Why use z at center ??? 
            zpd= z +dztop
            wt(i,j,k)= 0.d0
#ifdef sigma_stretch
            wzk(i,j,k)= -epm1*dnf/pfac/( epm1*(z-Df)+(D(i,j)+dzbot-Df) ) !!!
#else
            wzk(i,j,k)= dnf/(Df-D(i,j)-dzbot)
#endif
         end do ! end k

         ! -----------------------  X, Y FACE  -----------------------

         do 15 k=0,Nf+1 ! lower layer
            sig= dble(k)
            z= zf(i,j,k)
            zpd= z +dztop
#ifdef sigma_stretch
            temp= epm1*dnf/pfac/( epm1*(z-Df)+(D(i,j)+dzbot-Df) ) *(z-Df)/(D(i,j)+dzbot-Df) !!!
#else
            temp= -dnf/(Df-D(i,j)-dzbot) *(z-Df)/(D(i,j)+dzbot-Df)                       !!!
#endif
            wxk= Ddx(i,j)*temp
            wyk= Ddy(i,j)*temp
            gqk(i,j,k,1)= qpr*Jac(i,j,k)*(ux(i,j)*wxk +uy(i,j)*wyk)
            gqk(i,j,k,2)= qpr*Jac(i,j,k)*(vx(i,j)*wxk +vy(i,j)*wyk)
            gqk(i,j,k,3)= Jac(i,j,k)*(qpr*(wxk*wxk +wyk*wyk) +       &
   &           be2*wz(i,j,k)*wz(i,j,k))
   15    continue

      enddo ! end i
   enddo ! end j

! ------------------------------------------------------------------
!							Jacobians
! ------------------------------------------------------------------

   do 19 k=1,NK 
      do 21 i=0,NI 
         do 31 j=1,NJ 
            Jifc(i,j,k)= 0.5d0*(Jac(i,j,k)+ Jac(i+1,j,k)) 
            gi(i,j,k,1)= 0.5d0*(g11(i,j) +g11(i+1,j))*Jifc(i,j,k) 
            gi(i,j,k,2)= 0.5d0*(g12(i,j) +g12(i+1,j))*Jifc(i,j,k) 
            gqi(i,j,k,1)= qpr*gi(i,j,k,1) 
            gqi(i,j,k,2)= qpr*gi(i,j,k,2) 
   31    continue 
   21 continue 
   19 continue 
   
   do 28 k=1,NK 
      do 22 i=1,NI 
         do 32 j=0,NJ 
            Jjfc(i,j,k)= 0.5d0*(Jac(i,j,k)+ Jac(i,j+1,k)) 
            gj(i,j,k,1)= 0.5d0*(g12(i,j) +g12(i,j+1))*Jjfc(i,j,k) 
            gj(i,j,k,2)= 0.5d0*(g22(i,j) +g22(i,j+1))*Jjfc(i,j,k) 
            gqj(i,j,k,1)= qpr*gj(i,j,k,1) 
            gqj(i,j,k,2)= qpr*gj(i,j,k,2) 
   32    continue 
   22 continue 
   28 continue 

   do j=1,NJ 
      do i=0,NI 
         do k=1,NK 
            gi3(i,j,k)= 0.5d0*(g13(i,j,k) +g13(i+1,j,k))*Jifc(i,j,k) 
            gqi3(i,j,k)= qpr*gi3(i,j,k) 
         enddo
      enddo
   enddo

   do j=0,NJ 
      do i=1,NI 
         do k=1,NK 
            gj3(i,j,k)= 0.5d0*(g23(i,j,k) +g23(i,j+1,k))*Jjfc(i,j,k) 
            gqj3(i,j,k)= qpr*gj3(i,j,k) 
         enddo
      enddo
   enddo
 
!  ------------------------------------------------------------------
!   	  Bottom boundary layers - evenly spaced, no stretching
!  ------------------------------------------------------------------
#ifdef fixed_bottom_thickness

      do j=0,NJ+1
         do i=0,NI+1

            do k=0,1        !     at cell centers
               z= zc(i,j,k)
               wz(i,j,k)= 1./dzbot
               Jac(i,j,k)= J2d(i,j)/wz(i,j,k) !  z is  dimensionaless and is equal to  z*/DL
               sig= dble(k)-0.5d0
               wx(i,j,k)= -wz(i,j,k)*Ddx(i,j)
               wy(i,j,k)= -wz(i,j,k)*Ddy(i,j)
               g13(i,j,k)= ux(i,j)*wx(i,j,k) +uy(i,j)*wy(i,j,k)
               g23(i,j,k)= vx(i,j)*wx(i,j,k) +vy(i,j)*wy(i,j,k)
               ! g33(i,j,k)= wx(i,j,k)*wx(i,j,k) +wy(i,j,k)*wy(i,j,k) +
            end do

            do k=0,1-1      ! at cell faces
               sig= dble(k)
               wt(i,j,k)= 0.d0
               z= zc(i,j,k)
               zpd= z +dztop
               wzk(i,j,k)= 1./dzbot
            end do

            wzk(i,j,1)= 0.5*(wz(i,j,1) + wz(i,j,1+1))

            do k=0,1
               z= zf(i,j,k)
               sig= dble(k)
               temp= 1./dzbot
               wxk= -Ddx(i,j)*temp
               wyk= -Ddy(i,j)*temp
               gqk(i,j,k,1)= qpr*Jac(i,j,k)*(ux(i,j)*wxk +uy(i,j)*wyk)
               gqk(i,j,k,2)= qpr*Jac(i,j,k)*(vx(i,j)*wxk +vy(i,j)*wyk)
               gqk(i,j,k,3)= Jac(i,j,k)*(qpr*(wxk*wxk +wyk*wyk) + be2*wzk(i,j,k)*wzk(i,j,k))
            end do

         end do ! end i
      end do ! end j

      ! interpolation of Jac, g11, g12 and gi on the face x. 
      do k=1,1
         do i=0,NI
            do j=1,NJ
               Jifc(i,j,k)= 0.5d0*(Jac(i,j,k)+ Jac(i+1,j,k))
               gi(i,j,k,1)= 0.5d0*(g11(i,j) +g11(i+1,j))*Jifc(i,j,k)
               gi(i,j,k,2)= 0.5d0*(g12(i,j) +g12(i+1,j))*Jifc(i,j,k)
               gqi(i,j,k,1)= qpr*gi(i,j,k,1) 
               gqi(i,j,k,2)= qpr*gi(i,j,k,2) 
            end do
         end do
      end do

      ! interpolation of Jac, g11, g12 and gi on the face y. 
      do k=1,1
         do i=1,NI
            do j=0,NJ
               Jjfc(i,j,k)= 0.5d0*(Jac(i,j,k)+ Jac(i,j+1,k))
               gj(i,j,k,1)= 0.5d0*(g12(i,j) +g12(i,j+1))*Jjfc(i,j,k)
               gj(i,j,k,2)= 0.5d0*(g22(i,j) +g22(i,j+1))*Jjfc(i,j,k)
               gqj(i,j,k,1)= qpr*gj(i,j,k,1) 
               gqj(i,j,k,2)= qpr*gj(i,j,k,2) 
            end do
         end do
      end do

#endif ! end ifdef fixed_bottom_thickness

!  ------------------------------------------------------------------
!							Jacobian inv
!  ------------------------------------------------------------------

   do i=0,NI
      do j=0,NJ
         do k=0,NK
            Jacinv(i,j,k)=1./Jac(i,j,k)
         enddo
      enddo
   enddo

   
return 
END                                           
