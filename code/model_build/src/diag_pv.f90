subroutine diag_pv(n)


  !     ---------------------------------------------------               
  USE header
  !     computes the vorticity (for now just the vertical component)      
  !     vor in (per second)                                               
  !     PV in (per second^3)                                              
  implicit NONE 
      integer i,j,k,n,kp1,km1 
      REAL(kind=rc_kind) :: vdx,udy,udz,vdz,wdz,wdx,wdy,LENinv,grav,   &
     &     rdx,rdy,rdz,flocal,dz,dzinv                                  
                                                                        
      grav= gpr*10.d0 
                                                                        
!     j=0,NJ                                                            
      do i=1,NI 
         do k=1,NK 
            u(i,0,k,n)= 2.0*u(i,1,k,n) -u(i,2,k,n) 
            v(i,0,k,n)= 2.0*v(i,1,k,n) -v(i,2,k,n) 
            w(i,0,k,n)= 2.0*w(i,1,k,n) -w(i,2,k,n) 
            rho(i,0,k)= 2.0*rho(i,1,k) -rho(i,2,k) 
!                                                                       
            u(i,NJ+1,k,n)= 2.0*u(i,NJ,k,n) -u(i,NJ-1,k,n) 
            v(i,NJ+1,k,n)= 2.0*v(i,NJ,k,n) -v(i,NJ-1,k,n) 
            w(i,NJ+1,k,n)= 2.0*w(i,NJ,k,n) -w(i,NJ-1,k,n) 
            rho(i,NJ+1,k)= 2.0*rho(i,NJ,k) -rho(i,NJ-1,k) 
         end do 
      end do 
!     periodic-ew boundaries                                            
      do k=0,NK+1 
         do j=0,NJ+1 
            u(0,j,k,n)= u(NI,j,k,n) 
            v(0,j,k,n)= v(NI,j,k,n) 
            w(0,j,k,n)= w(NI,j,k,n) 
            rho(0,j,k)= rho(NI,j,k) 
!                                                                       
            u(NI+1,j,k,n)= u(1,j,k,n) 
            v(NI+1,j,k,n)= v(1,j,k,n) 
            w(NI+1,j,k,n)= w(1,j,k,n) 
            rho(NI+1,j,k)= rho(1,j,k) 
         end do 
      end do 
                                                                        
      DLinv= 1.0/DL 
      LENinv=1.0/LEN 
!     WL= eps*delta*UL  ! defined in main                               
      do k=1,NK 
         kp1=k+1 
         km1= k-1 
         if (k.eq.NK) then 
            kp1= k 
         else if (k.eq.1) then 
            km1= k 
         end if 
         dzinv= 1.d0/dble(kp1-km1) 
         do j=1,NJ 
            do i=1,NI 
               dz = DL*(zc(i,j,kp1)-zc(i,j,km1)) 
                                                                        
               vdx= 0.5*((v(i+1,j,k,n)-v(i-1,j,k,n))*ux(i,j) +          &
     &              (v(i,j+1,k,n) -v(i,j-1,k,n))*vx(i,j) )+             &
     &              dzinv*(v(i,j,kp1,n) -v(i,j,km1,n))*wx(i,j,k)        
               udy= 0.5*((u(i+1,j,k,n)-u(i-1,j,k,n))*uy(i,j) +          &
     &              (u(i,j+1,k,n) -u(i,j-1,k,n))*vy(i,j) )+             &
     &              dzinv*(u(i,j,kp1,n) -u(i,j,km1,n))*wy(i,j,k)        
               vor3(i,j,k)= (vdx -udy)*UL*LENinv 
                                                                        
               wdy=0.5*((w(i+1,j,k,n)-w(i-1,j,k,n))*uy(i,j) +           &
     &              (w(i,j+1,k,n) -w(i,j-1,k,n))*vy(i,j)) +             &
     &              dzinv*(w(i,j,kp1,n) -w(i,j,km1,n))*wy(i,j,k)        
               vdz= UL*(v(i,j,kp1,n)-v(i,j,km1,n))/dz 
               vor1(i,j,k)= wdy*WL*LENinv - vdz 
                                                                        
               udz= UL*(u(i,j,kp1,n)-u(i,j,km1,n))/dz 
               wdx= 0.5*((w(i+1,j,k,n)-w(i-1,j,k,n))*ux(i,j) +          &
     &              (w(i,j+1,k,n) -w(i,j-1,k,n))*vx(i,j) )  +           &
     &              dzinv*(w(i,j,kp1,n) -w(i,j,km1,n))*wx(i,j,k)        
               vor2(i,j,k)= udz -wdx*WL*LENinv 
            end do 
         end do 
      end do 
                                                                        
!     CALCULATE PV                                                      
      do k=1,NK 
         kp1=k+1 
         km1= k-1 
         if (k.eq.NK) then 
            kp1= k 
         else if (k.eq.1) then 
            km1= k 
         end if 
         dzinv = 1.d0/dble(kp1-km1) 
         do j=1,NJ 
            do i=1,NI 
                                                                        
               rdx = 0.5*LENinv*                                        &
     &              ((rho(i+1,j,k)-rho(i-1,j,k))*ux(i,j) +              &
     &              (rho(i,j+1,k) -rho(i,j-1,k))*vx(i,j) )+             &
     &              LENINV*dzinv*                                       &
     &              (rho(i,j,kp1) -rho(i,j,km1))*wx(i,j,k)              
                                                                        
               rdy= 0.5*LENinv*                                         &
     &              ((rho(i+1,j,k)-rho(i-1,j,k))*uy(i,j) +              &
     &              (rho(i,j+1,k) -rho(i,j-1,k))*vy(i,j) ) +            &
     &              dzinv*LENINV*(rho(i,j,kp1) -rho(i,j,km1))*wy(i,j,k) 
                                                                        
               rdz= dzinv*(rho(i,j,kp1) -rho(i,j,km1))*wz(i,j,k)*DLinv 
                                                                        
               flocal = FPAR*ffc(i,j) 
               pv3(i,j,k) = -(vor3(i,j,k)+flocal)*rdz*grav/R0 
               pv2(i,j,k) = -vor2(i,j,k)*rdy*grav/R0 
               pv1(i,j,k) = -vor1(i,j,k)*rdx*grav/R0 
                                                                        
            end do 
         end do 
      end do 
                                                                        
!=      call outpv(vor1,vor2,vor3,pv1,pv2,pv3)                          
                                                                        
      END                                           
