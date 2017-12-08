subroutine intpol 
!----------------------------------------------------                   
!     Use 3rd Order accurate interpolation                              
!     interpolate cx,cy,cz onto the cell faces : cxf                    
!     doesn't use 2nd order QUICK scheme of interpolation now           
!     n is the previous time level                                      
  USE header
      integer i,j,k,i0,im1,ip1,ip2,j0,jm1,jp1,jp2 
      REAL(kind=rc_kind) :: cx1,cx2,cy1,cy2,Jack,xa,xb 
!                                                                       
!     we need to impose boundary conditions like uu_x+vu_y+wu_z=0       
!     at i=0 and NI solid boundaries. If we use just interpolation      
!     to get these at the faces, we will not get the right value.       
!                                                                       
      if ((NI.lt.4).or.(NJ.lt.4)) then 
         write(6,*) 'prob in intpol,insuff i pts' 
         stop 
      end if 
!                                                                       
      xa= 9./16. 
      xb= -1./16. 
!     periodic ew boundaries                                            
      do i=0,NI 
         if (i.eq.0) then 
            i0= NI 
            im1=NI-1 
            ip1= i+1 
            ip2= i+2 
         else if (i.eq.1) then 
            i0= i 
            im1=NI 
            ip1= i+1 
            ip2= i+2 
         else if (i.eq.NI) then 
            i0= i 
            im1= i-1 
            ip1= 1 
            ip2= 2 
         else if (i.eq.(NI-1)) then 
            i0= i 
            im1= i-1 
            ip1= i+1 
            ip2= 1 
         else 
            i0= i 
            im1= i-1 
            ip1= i+1 
            ip2= i+2 
         end if 
                                                                        
         do k=1,NK 
            do j=1,NJ 
              !     do for face to the east of the i-th cell                          
              cxf(i,j,k)= ( xa*(ux(i0,j)*cx(i0,j,k)   + ux(ip1,j)*cx(ip1,j,k))    &
                       & +  xb*(ux(im1,j)*cx(im1,j,k) + ux(ip2,j)*cx(ip2,j,k))    &
                       & + 0.5*(uy(i0,j)*cy(i0,j,k)   + uy(ip1,j)*cy(ip1,j,k)) )  &
                       & *Jifc(i,j,k)                                        
          end do 
        end do 
      end do 
!                                                                       
      if (periodicew) then 
         continue 
      else 
         do 40 k=1,NK 
            do 50 j=1,NJ 
               cxf(0,j,k)= ufbcw(j,k) 
               cxf(NI,j,k)= ufbce(j,k) 
   50       continue 
   40    continue 
      endif 
!                                                                       
!                                                                       
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
              !     do for face to the north of the j-th cell                         
              cyf(i,j,k)= Jjfc(i,j,k)* (                                      &
                    &   0.5*(vx(i,j0)*cx(i,j0,k)   + vx(i,jp1)*cx(i,jp1,k))   &
                    &  + xa*(vy(i,j0)*cy(i,j0,k)   + vy(i,jp1)*cy(i,jp1,k))   &
                    &  + xb*(vy(i,jm1)*cy(i,jm1,k) + vy(i,jp2)*cy(i,jp2,k)) )                            
          end do 
        end do 
      end do 
!                                                                       
      do 140 k=1,NK 
         do 150 i=1,NI 
            cyf(i,NJ,k)= vfbcn(i,k) 
            cyf(i,0,k)= vfbcs(i,k) 
  150    continue 
  140 continue 
!                                                                       
!                                                                       
      do 210 k=1,NK-1 
         do 220 j=1,NJ 
            do 230 i=1,NI 
               Jack= 0.5d0*(Jac(i,j,k) + Jac(i,j,k+1)) 
               czf(i,j,k)= 0.5*Jack*(wx(i,j,k)*cx(i,j,k) + wx(i,j,k+1)*cx(i,j,k+1)    &
                          &  +       wy(i,j,k)*cy(i,j,k) +wy(i,j,k+1)*cy(i,j,k+1)     & 
                          &  + EPS*( wz(i,j,k)*cz(i,j,k)+wz(i,j,k+1)*cz(i,j,k+1)) )                         
!     &              EPS*wz(i,j)*( cz(i,j,k) +cz(i,j,k+1)) + wt(i,j,k) +
!     &              wt(i,j,k+1) )                                      
  230       continue 
  220    continue 
  210 continue 
!                                                                       
!     czf would be zero at k=0,NK impermeable boundaries                
      do 240 j=1,NJ 
         do 250 i=1,NI 
            czf(i,j,0)= wfbcb(i,j) 
            czf(i,j,NK)= Jac(i,j,NK)*( wx(i,j,NK)*cx(i,j,NK) +          &
     &           wy(i,j,NK)*cy(i,j,NK) +                                &
     &           EPS*wz(i,j,NK)*cz(i,j,NK) )                            
!cc     &           EPS*wz(i,j)*cz(i,j,NK) + wt(i,j,NK) )               
!c            czf(i,j,NK)= wfbct(i,j)                                   
  250    continue 
  240 continue 
!                                                                       
      return 
      END                                           
