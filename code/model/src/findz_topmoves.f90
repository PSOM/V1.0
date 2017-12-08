subroutine findzall 
  !     -------------------------------------------------------           
#include "cppdefs.h"
  USE header
  !     finds the value of z (non-dim by DL) given the value of sigma.    
  !     ht and dep are non-dim by DL.                                     
  !     At every (i,j), the column is divided into NK equal-depth cells.  
  integer i,j,k 
  REAL(kind=rc_kind) :: sigma,dep,ht,epm1,dnkm1,epm1inv,hpd, dnkm1inv,xfac                                                


 !---------------------------

  ! around NK...
  dnkm1= dble(NK-1) 
  dnkm1inv= 1.d0/dnkm1 

  ! stretching of the vertical grid
  !    pfac is the stretching in z. higher pfac gives more points near surf.
  pfac= 2.0d0 
  !      pfac= 5.d0  !c=NK32,dztop=0.5m                                   
  !==      pfac= 4.d0  !c=NK32,dztop=0.5m USED FOR SEVERAL MLI runs       
  !=      pfac=3.d0  !c used with NK=32 for first set of runs             
  epm1= exp(pfac) -1.d0 
  epm1inv= 1.d0/(exp(pfac) -1.d0) 


 !---------------------------

  do j=0,NJ+1 
    do i=0,NI+1 

      !  In the surface layer                                              

      hpd= h(i,j)*HDL +dztop 
      do k=NK,NK+1 
        sigma= dble(k)-0.5 
        zc(i,j,k)= (sigma -dnkm1)*hpd -dztop 
      end do
      do k=NK-1,NK+1 
        sigma= dble(k) 
        zf(i,j,k)= (sigma -dnkm1)*hpd -dztop 
      end do
   enddo
  enddo


#ifdef fixed_bottom_thickness

  dnkm1= dble(NK-1-1) 
  dnkm1inv= 1.d0/dnkm1 
#endif


  do j=0,NJ+1 
    do i=0,NI+1 
        
     !  Below the surface layer                                           

      do k=0,NK-1 
        sigma= dble(k)-0.5
#ifdef fixed_bottom_thickness
        sigma= dble(k-1)-0.5
#endif 
        xfac = (dnkm1 -sigma)*dnkm1inv 
 
#ifdef fixed_bottom_thickness
       zc(i,j,k)= (exp(pfac*xfac)-1.d0)*epm1inv*(D(i,j)+dzbot+dztop) -dztop                               
#else
       zc(i,j,k)= (exp(pfac*xfac)-1.d0)*epm1inv*(D(i,j)+dztop) -dztop                               
#endif 

      end do
      do k=-1,NK-2 
        sigma= dble(k) 
#ifdef fixed_bottom_thickness
        sigma= dble(k-1)
#endif 
        xfac = (dnkm1 -sigma)*dnkm1inv 

#ifdef fixed_bottom_thickness
        zf(i,j,k)= (exp(pfac*xfac)-1.d0)*epm1inv*(D(i,j)+dzbot+dztop) -dztop                               
#else
        zf(i,j,k)= (exp(pfac*xfac)-1.d0)*epm1inv*(D(i,j)+dztop) -dztop                               

#endif 
      end do


#ifdef fixed_bottom_thickness

!     For the bottom boundary layer
!     -------------------------------
            do k=0,1
               sigma= dble(k)-0.5
               zc(i,j,k)= ((sigma -0.)/dble(1))*dzbot+D(i,j)
            end do
            do k=-1,1
               sigma= dble(k)
               zf(i,j,k)= ((sigma-0.)/dble(1))*dzbot+D(i,j)
            end do
#endif
         end do
      end do



 return 

end subroutine findzall
