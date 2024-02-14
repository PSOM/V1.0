subroutine ini_topog 
#include "cppdefs.h"
  !   DEFINES Wavy topography -- sin(2pi x/ L) 
  !     ----------------                                                  
  USE header
  !     dy is the y grid spacing in m                                     
  !     creates the bottom topography                                     
  !     fits a spline to the described topog and computed Ddx(i)          
      integer i,j, nwaves
      REAL(kind=rc_kind) :: dep,LbyD,xlen,xmeter,wavelen,deriv

      xlen= dx*dble(NI)   ! length of channel (x-dimension) in meters
      nwaves = xlen / topo_wavelength ! nwaves is integer = the number of wavelengths that fit in the channel 
      wavelen = xlen / nwaves 

      LbyD= LEN/DL 
      DLinv=1.d0/DL 


      do i=0,NI+1 
         xmeter= xc(i)*1000.0   ! x-coord in meters
!         dep= depmean_dim*(1. + 0.1 *dsin(10.*2.*pi*xc(i)/xlen))
!         deriv =  0.1*depmean_dim * dcos(10.*2.*pi*xc(i)/xlen)*( 10.*2.*pi/xlen)
         dep  = -depmean_dim - topo_amplitude* dsin(2.*pi*xmeter/wavelen)
         deriv= -topo_amplitude* (2.*pi/wavelen)* dcos(2.*pi*xmeter/wavelen)
         ! write(6,*) j,yc(j),dep                                       

         do j=0,NJ+1 
            D(i,j)= dep*DLinv 
!            Ddx(i,j)= deriv*LbyD*1.d-3 ! convert from km to m and non-dimensionalizing???
            Ddx(i,j)= deriv*LbyD
            Ddy(i,j)= 0.d0                                      
         end do
         j=1
!          write(6,*) i,xc(i), D(i,j),Ddx(i,j)
  101 end do 


#ifdef file_output
      open (unit=70, file=TRIM(dirout)//'Ddxdy.dat')
      do j=0,NJ+1
         write(70,*) 'j= ',j
         do i=0,NI+1
            write(70,*) D(i,j),Ddx(i,j),Ddy(i,j)
        end do
      end do
      close(70)
#endif             
      return 
      END                                           