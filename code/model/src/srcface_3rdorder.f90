subroutine srcface(n,step) 
!----------------------------------------------------                   
  USE header
!     FOR PERIODICEW boundaries                                         
!     --------------------------                                        
!     We interpolate the source terms onto the cell faces               
!     It is important that u,v,rp have been filled outside the boundarie
   implicit none
   integer i,j,k,n,in,inm1,step,iter 
   integer ip1,jp1,ip2,jp2,im1,jm1,i0,j0 
   REAL(kind=rc_kind) :: uint,vint,wint,uxi(0:NI,NJ),uyi(0:NI,NJ),        &
  &     vxj(NI,0:NJ),vyj(NI,0:NJ),rpxi,rpeta,rpsig,fac,fac2,         &
  &     ainv,be2,vcif,vcjf, Jack                                     
   REAL(kind=rc_kind) :: hxi, heta, gradh, hy, py,px 
! the following is not defined in original code, cause problem after changing implicit REAL(kind=rc_kind) :: to implicit none
   REAL(kind=rc_kind) ::  xa, xb, wxsk, wysk, wzsk, kaph1
                                                                        
!     We are using the values at the ghost points.                      
                                                                       

   fac= EPS*delta 
   ainv= 1.d0/apr 
   fac2= EPS*lambda 
   be2= beta*EPS*EPS 


    do j=1,NJ 
      do i=1,NI-1 
        uxi(i,j)= 0.5*(ux(i+1,j)+ux(i,j)) 
        uyi(i,j)= 0.5*(uy(i+1,j)+uy(i,j)) 
      end do 
      !  periodic e-w boundaries                                           
      uxi(NI,j)= 0.5*(ux(1,j)+ux(NI,j)) 
      uyi(NI,j)= 0.5*(uy(1,j)+uy(NI,j)) 
      uxi(0,j) = uxi(NI,j) 
      uyi(0,j) = uyi(NI,j) 
   end do 
                                                                        
   do i=1,NI 
     xa= 9./16. 
     xb= -1./16. 
     !     periodic ew boundaries                                            
     if (i.eq.1) then 
        i0= i 
        im1=NI 
        ip1= i+1 
        ip2= i+2 
     else if (i.eq.(NI-1)) then 
        i0= i 
        im1= i-1 
        ip1= i+1 
        ip2= 1 
     else if (i.eq.NI) then 
        i0= i 
        im1= i-1 
        ip1= 1 
        ip2= 2 
     else 
        i0= i 
        im1= i-1 
        ip1= i+1 
        ip2= i+2 
     end if 
                                                                    
                                                                    
     do k=1,NK 
       do j=1,NJ 
         uint= xa*(u(i0,j,k,n) + u(ip1,j,k,n)) + xb*(u(im1,j,k,n) + u(ip2,j,k,n))                  
         vint= xa*(v(i0,j,k,n) + v(ip1,j,k,n)) + xb*(v(im1,j,k,n) + v(ip2,j,k,n))                  
         wint= xa*(w(i0,j,k,n) + w(ip1,j,k,n)) + xb*(w(im1,j,k,n) + w(ip2,j,k,n))                  
         vcif= EPS*(xa*(uvis(i0,j,k) + uvis(ip1,j,k)) + xb*(uvis(im1,j,k) + uvis(ip2,j,k)))               
         vcjf= EPS*(xa*(vvis(i0,j,k) + vvis(ip1,j,k)) + xb*(vvis(im1,j,k) + vvis(ip2,j,k)))               
                                                                  
         px=         (p(ip1,j,k) -p(i,j,k))                            * gqi(i,j,k,1)    &
            & + 0.25*(p(ip1,j+1,k)+p(i,j+1,k)-p(ip1,j-1,k)-p(i,j-1,k)) * gqi(i,j,k,2)  &
            & + 0.25*(p(ip1,j,k+1)+p(i,j,k+1)-p(ip1,j,k-1)-p(i,j,k-1)) * gqi3(i,j,k)                
                                                                  
         sifc(i,j,k)= ((uxi(i,j)*(-ffi(i,j)*vint + fac*bbi(i,j)*wint - vcif)                  &
                    & + uyi(i,j)*( ffi(i,j)*uint                     - vcjf))*Jifc(i,j,k) + grpifc(i,j,k))  + px 
       end do 
     end do 
   end do 
!                                                                       
!     periodic-ew boundaries                                            
!                                                                       

   do k=1,NK 
     do j=1,NJ 
       sifc(0,j,k)= sifc(NI,j,k) 
     end do 
   end do 
                                                                       
!                                                                       
    do j=0,NJ 
      do i=1,NI 
        vxj(i,j)= 0.5*(vx(i,j+1)+vx(i,j)) 
        vyj(i,j)= 0.5*(vy(i,j+1)+vy(i,j)) 
      end do 
    end do 
                                                                        
    do j=1,NJ-1 
      if (j.eq.1) then 
        j0= j 
        jm1=NJ 
        jp1= j+1 
        jp2= j+2 
        xa = 0.5 
        xb = 0.0 
      else if (j.eq.NJ-1) then 
        j0= j 
        jm1= j-1 
        jp1= NJ 
        jp2= 1 
        xa = 0.5 
        xb = 0.0 
      else 
        j0= j 
        jm1= j-1 
        jp1= j+1 
        jp2= j+2 
        xa= 9./16. 
        xb= -1./16. 
      end if 
      do k=1,NK 
        do i=1,NI 
           vint= xa * (v(i,j0,k,n)+v(i,jp1,k,n)) + xb * (v(i,jm1,k,n)+v(i,jp2,k,n))                      
           uint= xa * (u(i,j0,k,n)+u(i,jp1,k,n)) + xb * (u(i,jm1,k,n)+u(i,jp2,k,n))                      
           wint= xa * (w(i,j0,k,n)+w(i,jp1,k,n)) + xb * (w(i,jm1,k,n)+w(i,jp2,k,n))                      
           vcif= EPS*(xa*(uvis(i,j0,k)+uvis(i,jp1,k)) + xb*(uvis(i,jm1,k)+uvis(i,jp2,k)) )                  
           vcjf= EPS*(xa*(vvis(i,j0,k)+vvis(i,jp1,k)) + xb*(vvis(i,jm1,k)+vvis(i,jp2,k)) )                  

           py=              (p(i,j+1,k)-p(i,j,k))                        * gqj(i,j,k,2) &
              & + 0.25*(p(i+1,j+1,k)+p(i+1,j,k)-p(i-1,j+1,k)-p(i-1,j,k)) * gqj(i,j,k,1) &
              & + 0.25*(p(i,j+1,k+1)+p(i,j,k+1)-p(i,j+1,k-1)-p(i,j,k-1)) * gqj3(i,j,k) 

           sjfc(i,j,k)= ((vxj(i,j)*(-ffj(i,j)*vint + fac*bbj(i,j)*wint - vcif)                 &
                      & + vyj(i,j)*( ffc(i,j)*uint                     - vcjf) ) * Jjfc(i,j,k) &
                      & + grpjfc(i,j,k)) + py                               

        end do 
      end do 
    end do


     do k=1,NK 
       do i=1,NI 
         do j=0,NJ,NJ 
           !     Use linear extrapolation over 2pts.  to get uint                  
           if (j.eq.0) then 
              in=1 
              inm1= 2 
           else if (j.eq.NJ) then 
              in= NJ 
              inm1 = NJ-1 
           end if 

           vint= 0.5d0*(3.d0*v(i,in,k,n)-v(i,inm1,k,n)) 
           uint= 0.5d0*(3.d0*u(i,in,k,n)-u(i,inm1,k,n)) 
           wint= 0.5d0*(3.d0*w(i,in,k,n)-w(i,inm1,k,n)) 
           vcif= EPS*uvis(i,in,k) 
           vcjf= EPS*vvis(i,in,k) 

           !     mgpfill should have been called before                            
           py=         (p(i,j+1,k)-p(i,j,k))                             * gqj(i,j,k,2)  & 
              & + 0.25*(p(i+1,j+1,k)+p(i+1,j,k)-p(i-1,j+1,k)-p(i-1,j,k)) * gqj(i,j,k,1)  &
              & + 0.25*(p(i,j+1,k+1)+p(i,j,k+1)-p(i,j+1,k-1)-p(i,j,k-1)) * gqj3(i,j,k) 
                                                                   
           sjfc(i,j,k)= ((  vxj(i,j)*(-ffj(i,j)*vint + fac*bbj(i,j)*wint - vcif)                  &
                   &      + vyj(i,j)*( ffc(i,j)*uint                     - vcjf) )*Jjfc(i,j,k)    &
                   &      + grpjfc(i,j,k)) +py                                
         end do 
       end do 
     end do 


 return 
                                                                        
 END                                          
 
