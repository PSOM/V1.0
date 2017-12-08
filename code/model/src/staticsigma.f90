subroutine staticsigma 

#include "cppdefs.h"
  USE header
  !     ------------------------                                          
  !     FOR THE REGION BELOW THE TOPMOST LAYER                            
  !     modified a for periodicew bc                                      
  !     This subroutine updates all those quantities that are a function  
  !     of time and  (time and depth). It is called  every time step.     
  !     The surface bc  wfbct is updated here.                            
!  implicit REAL(kind=rc_kind) :: (a-h,o-z) 
  integer i,j,k 
  REAL(kind=rc_kind) :: hpd,hpdinv,hu,hv,hx,hy,hxpdx,hypdy,z,temp,       &
       be2,wxk,wyk,d2,dnk,dnkm1,sig,ep,epm1,dnkm1p,zpd,           &
       g13(0:NI+1,0:NJ+1,0:NK+1),g23(0:NI+1,0:NJ+1,0:NK+1)          
  !     Ddx and Ddy are known from init                                   
  !                                                                       
  !     pfac is defined in findzall                                       
  dnk= dble(NK) 

#ifdef fixed_bottom_thickness
  dnkm1= dble(NK-1-1) 
#else
  dnkm1= dble(NK-1) 
#endif

  dnkm1p= dnkm1/pfac 
  be2= beta*EPS*EPS 
  epm1= exp(pfac) -1.d0 
  ep = pfac*exp(1.d0) 
  !                                                                       
  do j=0,NJ+1 
     do i=0,NI+1 
        !     All these variables are functions of time                         
            hpd= dztop +D(i,j) 
            hpdinv= 1.d0/hpd 
!-            if (i.eq.0) then                                          
!-               hu= HDL*(h(i+1,j)-h(NI-1,j))                           
!-            else if  (i.eq.NI+1) then                                 
!-               hu= HDL*(h(2,j)-h(i-1,j))                              
!-            else                                                      
!-               hu= 0.5d0*HDL*(h(i+1,j)-h(i-1,j))                      
!-            endif                                                     
!-            if (j.eq.0) then                                          
!-               hv= HDL*(h(i,j+1)-h(i,j))                              
!-            else if (j.eq.NJ+1) then                                  
!-               hv= HDL*(h(i,j)-h(i,j-1))                              
!-            else                                                      
!-               hv= 0.5d0*HDL*(h(i,j+1)-h(i,j-1))                      
!-            endif                                                     
!-            hx= hu*ux(i,j) + hv*vx(i,j)                               
!-            hy= hu*uy(i,j) + hv*vy(i,j)                               
!-            hxpdx= hx + Ddx(i,j)                                      
!-            hypdy= hy + Ddy(i,j)                                      
!     wz is not a function of depth when the sigma lines are equally spa
!     Hence wz is wz(i,j,time). For a stretched grid wz would be w(i,j,k
!     Then Jac which is now Jac(i,j,time) would  become  Jac(i,j,k,time)
!                                                                       
!            write(100,*) J2d(i,j)                                      
!            write(200,*) Jac(i,j)                                      
!                                                                       
            do 20 k=0,NK-1 
               z= zc(i,j,k) 
               zpd= z +dztop 
               wz(i,j,k)= -epm1*dnkm1p/(epm1*zpd +hpd) 
!             if ((i.eq.10).and.(j.eq.10)) write(6,*) k,wz(i,j,k),zpd   
               Jac(i,j,k)= J2d(i,j)/wz(i,j,k) 
!     All these variables are functions of time and depth               
!     z is  dimensionaless and is equal to  z*/DL                       
!c               z= ((dble(k)-0.5d0)/dble(NK)) *hpd - D(i,j)        
!c               temp= ( z+D(i,j))*hpdinv                               
!               wt(i,j,k)= -temp*HDL*hdt(i,j)*dble(NK)*hpdinv         
!     now hdt computed in vhydro already contains HDL                   
               sig= dble(k)-0.5d0 
!              wt is at cell faces                                      
!-               wt(i,j,k)= -sig*hdt(i,j)*hpdinv                        
!c               wt(i,j,k)= -temp*hdt(i,j)*dble(NK)*hpdinv            
!c               wx(i,j,k)= (Ddx(i,j) - temp*hxpdx)*                    
!c     &              dble(NK)*hpdinv                                 
!c               wy(i,j,k)= (Ddy(i,j) - temp*hypdy)*                    
!c     &              dble(NK)*hpdinv                                 
               wx(i,j,k)= wz(i,j,k)*zpd*hpdinv*Ddx(i,j) 
               wy(i,j,k)= wz(i,j,k)*zpd*hpdinv*Ddy(i,j) 
               g13(i,j,k)= ux(i,j)*wx(i,j,k) +uy(i,j)*wy(i,j,k) 
               g23(i,j,k)= vx(i,j)*wx(i,j,k) +vy(i,j)*wy(i,j,k) 
!               g33(i,j,k)= wx(i,j,k)*wx(i,j,k) +wy(i,j,k)*wy(i,j,k) +  
!     &              be2*wz(i,j)*wz(i,j)                                
   20       continue 
            do k=0,NK-1 
               sig= dble(k) 
               wt(i,j,k)= 0.d0 
               z= zc(i,j,k) 
               zpd= z +dztop 
               wzk(i,j,k)= -epm1*dnkm1p/(epm1*zpd +hpd) 
            end do 
            wzk(i,j,NK-1)= 0.5*(wz(i,j,NK) + wz(i,j,NK-1)) 
!     temporary fix ****** for k=0  - use linear extrapolation          
!-            wx(i,j,0)= 2.d0*wx(i,j,1) -wx(i,j,2)                      
!-            wz(i,j,0)= 2.d0*wz(i,j,1) -wz(i,j,2)                      
!-            wy(i,j,0)= 2.d0*wy(i,j,1) -wy(i,j,2)                      
!-            Jac(i,j,0)= 2.d0*Jac(i,j,1) -Jac(i,j,2)                   
!-            g13(i,j,0)= 2.d0*g13(i,j,1) -g13(i,j,2)                   
!-            g23(i,j,0)= 2.d0*g23(i,j,1) -g23(i,j,2)                   
!                                                                       
!     We evaluate gk at i=0,NI+1,j=0,NJ+1 only because these are used   
!     at k=0,NK to fill the horizontal edges in mgpfill or hnbc.        
            do 15 k=0,NK 
!c               z= (dble(k)/dble(NK))*hpd - D(i,j)                 
               z= zf(i,j,k) 
               zpd= z +dztop 
               sig= dble(k) 
!c               temp= ( z+D(i,j))*hpdinv                               
!               wtk= -temp*HDL*hdt(i,j)*dble(NK)*hpdinv               
!     now hdt computed in vhydro already contains HDL                   
!c               wtk= -temp*hdt(i,j)*dble(NK)*hpdinv                  
!c               wxk= (Ddx(i,j) - temp*hxpdx)*                          
!c     &              dble(NK)*hpdinv                                 
!c               wyk= (Ddy(i,j) - temp*hypdy)*                          
!c     &              dble(NK)*hpdinv                                 
!-               wtk= -sig*hdt(i,j)*hpdinv                              
                                                                        
               temp= (epm1*dnkm1p/(epm1*zpd +hpd))*zpd*hpdinv 
               wxk= Ddx(i,j)*temp 
               wyk= Ddy(i,j)*temp 
               gqk(i,j,k,1)= qpr*Jac(i,j,k)*(ux(i,j)*wxk +uy(i,j)*wyk) 
               gqk(i,j,k,2)= qpr*Jac(i,j,k)*(vx(i,j)*wxk +vy(i,j)*wyk) 
               gqk(i,j,k,3)= Jac(i,j,k)*(qpr*(wxk*wxk +wyk*wyk) +       &
     &              be2*wz(i,j,k)*wz(i,j,k))                            
   15       continue 
      enddo
    enddo

!      do j=0,NJ+1                                                      
!         do i=0,NI+1                                                   
!            write(120,*) gk(i,j,0,1), gk(i,j,3,1)                      
!            write(130,*) gk(i,j,0,2), gk(i,j,3,2)                      
!            write(140,*) gk(i,j,0,3), gk(i,j,3,3)                      
!         end do                                                        
!      end do                                                           
!      stop                                                             
!                                                                       
!                                                                       
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
!                                                                       
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

!                                                                       
!     vent recomputed using vfbcs                                       
!     v= [(-vx/d2)*uf +(ux/d2)*vf]/Jac , but uf=0                       
!      j=0                                                              
!      do 50 i=ient1,ient2                                              
!         d2= ux(i,j)*vy(i,j) -vx(i,j)*uy(i,j)                          
!        do 55 k=1,NK                                                   
!            vent(i,k)= (ux(i,j)/d2)*vfbcs(i,k)/Jac(i,j)                
! 55      continue                                                      
! 50   continue                                                         
!                                                                       
!-      write(6,*) 'enter k'                                            
!-      read(5,*) k                                                     
!      do k=0,NK+1                                                      
!         write(500,*) 'k= ',k                                          
!         do j=0,NJ+1                                                   
!            do i=0,NI+1                                                
!c               write(500,*) wx(i,j,k),wz(i,j,k)                       
!               write(500,*) gqi3(i,j,k),gqj3(i,j,k)                    
!               write(600,*) wy(i,j,k),Jac(i,j,k)                       
!            end do                                                     
!         end do                                                        
!      end do                                                           
!      write(6,*) 'stopping in staticsigma'                             
!      stop                                                             
!                   i

#ifdef fixed_bottom_thickness
!==================================================================================
!     Bottom boundary layers - evenly spaced, no stretching
!     -----------------------------------------------------
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

         end do
      end do

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


#endif

!i=48;j=10;
! do k=1,32
   ! PRINT"('QW ',7(E18.11,' '))",wx(i,j,k),wy(i,j,k),wz(i,j,k),g13(i,j,k),g23(i,j,k),Jac(i,j,k),Jjfc(i,j,k)
! enddo

 do i=0,NI
   do j=0,NJ
     do k=0,NK
       Jacinv(i,j,k)=1./Jac(i,j,k)
     enddo
   enddo
 enddo

                                                    
       return 
      END                                           
