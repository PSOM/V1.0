subroutine findzall 

#include "cppdefs.h"
USE header
! ------------------------------------------------------------------
!     finds the value of z (non-dim by DL) given the value of sigma.    
!     ht and dep are non-dim by DL.                                     
!     At every (i,j), the column is divided into NK equal-depth cells.  
! ------------------------------------------------------------------
   implicit none
   integer i,j,k
   REAL(kind=rc_kind) :: sigma,epm1,epm1inv,dnkm1,dnf,dnfinv,hpd,xfac,Df
 
#ifdef sigma_stretch 
   epm1= exp(pfac)-1.d0 
   epm1inv= 1.d0/epm1
#endif
   dnkm1= dble(NK-1) 
   dnf  = dble(Nf)
   Df=(-depmean_dim+distance)/DL
    
! ---------------------------------------------------------------------------------
!  								In the surface layers
! ---------------------------------------------------------------------------------

   do j=0,NJ+1 
      do i=0,NI+1 

         ! ---- CENTERS ----
         hpd= h(i,j)*HDL +dztop 
         do k=NK,NK+1 
            sigma= dble(k)-0.5 
            zc(i,j,k)= (sigma -dnkm1)*hpd -dztop 
         end do

         ! -----  FACES -----
         do k=NK-1,NK+1 
            sigma= dble(k) 
            zf(i,j,k)= (sigma -dnkm1)*hpd -dztop 
         end do

      enddo
   enddo

! ---------------------------------------------------------------------------------
!  							In the upper layers (rectangular)
! ---------------------------------------------------------------------------------

   do j=0,NJ+1 
      do i=0,NI+1 

         ! ---- CENTERS ----
         do k=Nf+1,NK-1
            sigma= dble(k)-0.5
            xfac = (dnkm1-sigma)/(dnkm1-Nf)
#ifdef sigma_stretch
            zc(i,j,k) = (exp(pfac*xfac)-1.d0)*epm1inv*(Df+dztop) - dztop
#else            
            zc(i,j,k) = xfac*(Df+dztop) - dztop
#endif
         end do

         ! ----- FACES -----
         do k=Nf,NK-2
            sigma= dble(k)
            xfac = (dnkm1-sigma)/(dnkm1-Nf)
#ifdef sigma_stretch            
            zf(i,j,k) = (exp(pfac*xfac)-1.d0)*epm1inv*(Df+dztop) - dztop
#else            
            zf(i,j,k) = xfac*(Df+dztop) - dztop
#endif
         end do

! ---------------------------------------------------------------------------------
!  							In the lower layers (corrugated)
! ---------------------------------------------------------------------------------	  
#ifdef fixed_bottom_thickness
   dnf= dble(Nf-1) 
   dnfinv= 1.d0/dnf 
#else
   dzbot = 0.
#endif

         ! ---- CENTERS ----
         do k=0,Nf
            sigma= dble(k)-0.5
#ifdef fixed_bottom_thickness
            sigma= dble(k-1)-0.5
#endif 
            xfac = (dnf-sigma)*dnfinv 
#ifdef sigma_stretch
            zc(i,j,k)= (exp(pfac*xfac)-1.d0)*epm1inv*(D(i,j)+dzbot-Df) +Df                               
#else
            zc(i,j,k)= xfac*(D(i,j)+dzbot-Df) +Df                               
#endif
         end do

         ! ----- FACES -----
         do k=-1,Nf-1 
            sigma= dble(k) 
#ifdef fixed_bottom_thickness
            sigma= dble(k-1)
#endif 
            xfac = (dnf-sigma)*dnfinv 
#ifdef sigma_stretch
            zf(i,j,k)= (exp(pfac*xfac)-1.d0)*epm1inv*(D(i,j)+dzbot-Df) +Df                               
#else
            zf(i,j,k)= xfac*(D(i,j)+dzbot-Df) +Df                               
#endif
         end do

! ---------------------------------------------------------------------------------
!  								In the bottom layers
! ---------------------------------------------------------------------------------

#ifdef fixed_bottom_thickness
         ! ---- CENTERS ----
         do k=0,1
            sigma= dble(k)-0.5
            zc(i,j,k)= sigma*dzbot+D(i,j)
         end do

         ! ----- FACES -----
         do k=-1,1
            sigma= dble(k)
            zf(i,j,k)= sigma*dzbot+D(i,j)
         end do
#endif

      end do ! end i
   end do ! end j


return 

end subroutine findzall
