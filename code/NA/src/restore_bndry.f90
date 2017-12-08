      subroutine restore_bndry(n)
!     ---------------------------------------------
!     restores the density below  MLD at the N-S boundaries
!     so as to counter the effect of wind stress curl - 
!     added on 4/21/11

      USE header
      implicit none
      integer i,j,k,n,nw1,nw2
      real(KIND=rc_kind) rhonorth(NK),rhosouth(NK),rhomean,zdepth,restoretime,fac

!     Wind stress curl is in the region 1 to nw1, and nw2 to NJ
!      nw1=njf1/2                 
!      nw2=njf4+(nj-njf4)/2
      nw1=0.1*nj
      nw2=0.9*nj


      do k=1,NK
         rhosouth(k)=0.d0
         rhonorth(k)= 0.d0
         do i=1,NI
            rhosouth(k)=rhosouth(k) + s(i,1,k,n)
            rhonorth(k)=rhonorth(k) + s(i,NJ,k,n)
         end do
         rhosouth(k)= rhosouth(k)/dble(NI)
         rhonorth(k)= rhonorth(k)/dble(NI)
      end do

!     Restore at southern boundary, below MLD
!     restore time = 3 days = 3*86400/TL
      restoretime=3.d0*86400.d0/TL
      fac= dtime/restoretime
      do k=1,NK
         do j=2,nw1
            rhomean=0.d0
            do i=1,NI
               rhomean= s(i,j,k,n)+ rhomean
            end do
            rhomean= rhomean/dble(NI)
            do i=1,NI
               zdepth=-zf(i,j,k)*DL
               if (zdepth.gt.mld(i,j)) then   !resotre the mean to rhosouth
                  s(i,j,k,n)= s(i,j,k,0)- fac*(rhomean-rhosouth(k))
               end if
            end do
         end do

         do j=nw2+1,NJ-1
            rhomean=0.d0
            do i=1,NI
               rhomean= s(i,j,k,n)+ rhomean
            end do
            rhomean= rhomean/dble(NI)
            do i=1,NI
               zdepth=-zf(i,j,k)*DL
               if (zdepth.gt.mld(i,j)) then   !resotre the mean to rhonorth
                  s(i,j,k,n)= s(i,j,k,0)- fac*(rhomean-rhonorth(k))
               end if
            end do
         end do

      end do
                  
      return
      end
