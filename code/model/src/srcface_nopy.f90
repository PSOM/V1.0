subroutine srcface_nopy(n) 
!----------------------------------------------------                   
  USE header
  !USE rpgrads
!     FOR PERIODICEW boundaries                                         
!     --------------------------                                        
!     We interpolate the source terms onto the cell faces               
!     It is important that u,v,rp have been filled outside the boundarie
      implicit none
      integer i,j,k,n,in 
      REAL(kind=rc_kind) ::  uint,vint,wint,uxi(0:NI,NJ),uyi(0:NI,NJ),        &
     &     vxj(NI,0:NJ),vyj(NI,0:NJ),rpxi,rpeta,rpsig,fac,fac2,         &
     &     ainv,be2,vcif,vcjf, Jack                                     
      REAL(kind=rc_kind) ::  hxi, heta, gradh, hy, py,px 
      REAL(kind=rc_kind) ::  wxsk,wysk,wzsk
                                                                        
!     We are using the values at the ghost points.                      
                                                                        
      fac= EPS*delta 
      ainv= 1.d0/apr 
      fac2= EPS*lambda 
      be2= beta*EPS*EPS 
!c      do 10 i=0,NI                                                    
      do j=1,NJ 
         do i=1,NI-1 
            uxi(i,j)= 0.5*(ux(i+1,j)+ux(i,j)) 
            uyi(i,j)= 0.5*(uy(i+1,j)+uy(i,j)) 
         end do 
!     periodic e-w boundaries                                           
         uxi(NI,j)= 0.5*(ux(1,j)+ux(NI,j)) 
         uyi(NI,j)= 0.5*(uy(1,j)+uy(NI,j)) 
         uxi(0,j) = uxi(NI,j) 
         uyi(0,j) = uyi(NI,j) 
      end do 
                                                                        
      do k=1,NK 
            do j=1,NJ 
               do i=1,NI-1 
                  vint= 0.5*(v(i+1,j,k,n)+v(i,j,k,n)) 
                  uint= 0.5*(u(i+1,j,k,n)+u(i,j,k,n)) 
                  wint= 0.5*(w(i+1,j,k,n)+w(i,j,k,n)) 
                  vcif= 0.5*EPS*(uvis(i+1,j,k)+uvis(i,j,k)) 
                  vcjf= 0.5*EPS*(vvis(i+1,j,k)+vvis(i,j,k)) 
                  sifc(i,j,k)= ( (   uxi(i,j)*(-ffi(i,j)*vint+ fac*bbi(i,j)*wint - vcif)  &
                          &        + uyi(i,j)*( ffi(i,j)*uint                    - vcjf))*Jifc(i,j,k) &
                          &     + grpifc(i,j,k))                                    
            end do 
         end do 
      end do 

!                                                                       
!     periodic-ew boundaries                                            
      i=NI 
      do j=1,NJ 
         do k=1,NK 
            vint= 0.5*(v(1,j,k,n)+v(NI,j,k,n)) 
            uint= 0.5*(u(1,j,k,n)+u(NI,j,k,n)) 
            wint= 0.5*(w(1,j,k,n)+w(NI,j,k,n)) 
            vcif= 0.5*EPS*(uvis(1,j,k)+uvis(NI,j,k)) 
            vcjf= 0.5*EPS*(vvis(1,j,k)+vvis(NI,j,k)) 
            sifc(i,j,k)=( (  uxi(i,j)*(-ffi(i,j)*vint +fac*bbi(i,j)*wint - vcif)                   &
                       &   + uyi(i,j)*( ffi(i,j)*uint                    - vcjf) ) * Jifc(i,j,k)   &
                       & + grpifc(i,j,k))                                       
         end do 
      end do 
!                                                                       
      do k=1,NK 
         do j=1,NJ 
            sifc(0,j,k)= sifc(NI,j,k) 
         end do 
      end do 
!                       

! sifc is computed.
! ------------------
                                                
!                                                                       
      do j=0,NJ 
         do i=1,NI 
            vxj(i,j)= 0.5*(vx(i,j+1)+vx(i,j)) 
            vyj(i,j)= 0.5*(vy(i,j+1)+vy(i,j)) 
         end do 
      end do 
                                                                        
      do k=1,NK 
         do i=1,NI 
            do j=1,NJ-1 
               vint= 0.5*(v(i,j+1,k,n)+v(i,j,k,n)) 
               uint= 0.5*(u(i,j+1,k,n)+u(i,j,k,n)) 
               wint= 0.5*(w(i,j+1,k,n)+w(i,j,k,n)) 
               vcif= 0.5*EPS*(uvis(i,j+1,k)+uvis(i,j,k)) 
               vcjf= 0.5*EPS*(vvis(i,j+1,k)+vvis(i,j,k)) 
               sjfc(i,j,k)= ( (  vxj(i,j)*(-ffj(i,j)*vint +fac*bbj(i,j) *wint - vcif)   &
                           &   + vyj(i,j)*( ffc(i,j)*uint                     - vcjf) )  * Jjfc(i,j,k) &
                           & + grpjfc(i,j,k))                                    
            end do 

            do j=0,NJ,NJ 
               if (j.eq.0) in=1 
               if (j.eq.NJ) in= NJ 
               vint= v(i,in,k,n) 
               uint= u(i,in,k,n) 
               wint= w(i,in,k,n) 
               vcif= EPS*uvis(i,in,k) 
               vcjf= EPS*vvis(i,in,k) 
               sjfc(i,j,k)= ( (  vxj(i,j)*(-ffj(i,j)*vint +fac*bbj(i,j)*wint - vcif)  &
                            &  + vyj(i,j)*( ffc(i,j)*uint                    - vcjf) ) * Jjfc(i,j,k) &
                            & + grpjfc(i,j,k))                                    

            end do 
         end do 
      end do 


! sjfc is computed.
! ------------------


      return 



!     ================                                                  
!     skfc computed in separate routine                                 
!                                                                       
      do 70 i=1,NI 
         do 80 j=1,NJ 
            do 90 k=1,NK-1 
               wxsk= 0.25d0*(wx(i,j,k+1) +wx(i,j,k))*(si(i,j,k+1)+      &
     &              si(i,j,k) )                                         
               wysk= 0.25d0*(wy(i,j,k+1) +wy(i,j,k))*(sj(i,j,k+1)+      &
     &              sj(i,j,k) )                                         
               wzsk= 0.25d0*(wz(i,j,k+1) +wz(i,j,k))*(sk(i,j,k+1)+      &
     &              sk(i,j,k) )                                         
               Jack= 0.5d0*(Jac(i,j,k+1) + Jac(i,j,k)) 
               skfc(i,j,k)= ( be2*wzsk +wxsk +wysk )*Jack 
!               skfc(i,j,k)= 0.                                         
   90       continue 
            k=0 
!     linear extrapolation for wzsk                                     
            wzsk= wzk(i,j,k)*0.5*(3.0*sk(i,j,k+1)-sk(i,j,k+2)) 
!            wzsk= wzk(i,j,k)*sk(i,j,k+1)                               
            wxsk= wx(i,j,k+1)*si(i,j,k+1) 
            wysk= wy(i,j,k+1)*sj(i,j,k+1) 
            Jack= Jac(i,j,k+1) 
            skfc(i,j,k)= ( be2*wzsk +wxsk +wysk )*Jack 
!            skfc(i,j,k)= 0.0                                           
            k= NK 
!     linear extrap                                                     
            wzsk= wzk(i,j,k)*0.5*(3.0*sk(i,j,k)-sk(i,j,k-1)) 
!            wzsk= wzk(i,j,k)*sk(i,j,k)                                 
            wxsk= 0.5d0*(wx(i,j,k+1) +wx(i,j,k))*si(i,j,k) 
            wysk= 0.5d0*(wy(i,j,k+1) +wy(i,j,k))*sj(i,j,k) 
            Jack= Jac(i,j,k) 
            skfc(i,j,k)= ( be2*wzsk +wxsk +wysk )*Jack 
!            skfc(i,j,k)= 0.0                                           
                                                                        
! july 11, 2001;  compute skfc more accurately at top boundary          
!c            skfc(i,j,0)= skfc(i,j,1)                                  
!c            skfc(i,j,NK)= skfc(i,j,NK-1)                              
   80    continue 
   70 continue 
                                                                        
      return 
      END                                           
