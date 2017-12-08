      subroutine conadjust(stepl,n)
!----------------------------------------------------
!     T(i,j,k) is set to 0, so Treal is T0
!     Performs convective adjustment
      USE header
      implicit none
      integer i,j,k,kk,n,stepl,iclean,icount
      double precision rh,rhtop,zt,zb,dz,zinv
      double precision dzupper
      double precision sreal,Treal,rhreal,pbar
!      double precision rho(0:NI+1,0:NJ+1,0:NK+1)
!
!-      call evalrho(rho,n)

      do j=1,NJ
         do i=1,NI
            icount=1
!     
 101        iclean=1
            do k=NK-1,1,-1
!               conv(i,j,k)= 0
               rhtop= s(i,j,k+1,n)           !top 
!     this level
               rh =  s(i,j,k,n)

               if (rh.lt.rhtop) then
                  iclean=0
                  dz = zf(i,j,k) -zf(i,j,k-1)  !this level
                  dzupper = zf(i,j,k+1) -zf(i,j,k)   ! upper level
!     mix
                  conv(i,j,k)= 1
                  zinv= 1.d0/(dzupper +dz)
                  s(i,j,k,n)=(dzupper*s(i,j,k+1,n) +dz*s(i,j,k,n))*zinv
                  s(i,j,k+1,n)= s(i,j,k,n)
                  do it=1,ntr
                     Tr(it,i,j,k,n)=(dzupper*Tr(it,i,j,k+1,n) & 
     &                    +dz*Tr(it,i,j,k,n))*zinv
                     Tr(it,i,j,k+1,n)= Tr(it,i,j,k,n)
                  end do
               end if
            end do
            if ((iclean.eq.0).and.(icount.lt.4)) then 
!did not restrat            if ((iclean.eq.0).and.(icount.lt.7)) then 
!               write(6,*) 'icount', icount
               icount=icount+1
               goto 101
            end if
         end do
      end do

      if (mod((stepl-1),100).eq. 0) then
         do k=0,NK+1
            do j=0,NJ+1
               do i=0,NI+1
                  con100(i,j,k) = 0
               end do
            end do
         end do
      else
         do k=0,NK+1
            do j=0,NJ+1
               do i=0,NI+1
                  con100(i,j,k) = con100(i,j,k) + conv(i,j,k)
               end do
            end do
         end do
      end if
!     re-evaluate density (if needed)
!      call evalrho(rho,n)

!-      call sTbc_periodicew(n)
      return
      end

