subroutine uv_geostrophic(ugeo,vgeo)
!     ---------------------------------------------                     
  !   Finds the geostrophic velocities ugeo,vgeo
  USE header
!  USE rpgrads
!     modified for periodicew bc                                        
!     ---------------------------                                       
!     Sets up the initial velocities so that they are in geostrophic bal
!     with the initial density field.                                   
!      implicit logical (a-z)                                           
      implicit none 
      integer i,j,k,n,imax,jmax,kmax,m 
      real(kind=rc_kind):: uxi,vyj,hxi,heta,hx,hy,px,py,ujfc 
      real(kind=rc_kind):: res,resmax,kaph1
      real(kind=rc_kind):: ainv,be2,fac2,wzsk,wxsk,wysk,Jack,pxi,peta,      &
     &     psig,pgrad,con,pz                                            
      real(kind=rc_kind), dimension (0:NI+1,0:NJ+1,0:NK+1) :: ugeo, vgeo
                                                                        
      call rpevalgrad_Song(0) 
                                                                        
      kaph1= 1.d0-kappah 
      con  = 1.0 -qpr 
                                                                        
      do j=1,NJ 
         do i=1,NI 
            hxi  = 0.5d0*( h(i+1,j)-h(i-1,j) ) 
            heta = 0.5d0*( h(i,j+1)-h(i,j-1) ) 
            hx   = ux(i,j)*hxi +vx(i,j)*heta 
            hy   = uy(i,j)*hxi +vy(i,j)*heta 
            do k=1,NK 
               pxi  = 0.5d0*(p(i+1,j,k)-p(i-1,j,k)) 
               peta = 0.5d0*(p(i,j+1,k)-p(i,j-1,k)) 
               psig = 0.5d0*(p(i,j,k+1)-p(i,j,k-1)) 
               px   = (ux(i,j)*pxi +vx(i,j)*peta +wx(i,j,k)*psig) 
               py   = (uy(i,j)*pxi +vy(i,j)*peta +wy(i,j,k)*psig) 
                                                                        
               ugeo(i,j,k) = - (qpr*py +drpy(i,j,k)+gpr*hy)/             &
                    (ffc(i,j))                                          
               vgeo(i,j,k) = (qpr*px +drpx(i,j,k) +gpr*hx)/              &
                    (ffc(i,j))                                          
            end do
         end do 
         do k=1,NK 
            ugeo(0,j,k)= ugeo(NI,j,k) 
            vgeo(0,j,k)= vgeo(NI,j,k) 
            ugeo(NI+1,j,k) = ugeo(1,j,k) 
            vgeo(NI+1,j,k) = vgeo(1,j,k) 
         end do 
      end do 

      return
      end


                                                                        
