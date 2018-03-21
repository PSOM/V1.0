!! This subroutine calculates shear and buoyancy frequencies at cell faces, by a linear
!! interpolation from the cell centers

subroutine shearn2(ii,jj)

use header
implicit none

!!! Local Variables
integer :: ii,jj,k

REAL(kind=rc_kind):: temp, dtr0 

dtr0= -0.17

!!!! Calculate shear separately
do k=0, (NK-1)   
   ush(ii,jj,k) = (u(ii,jj,k+1,1) - u(ii,jj,k,1))*UL / ((zc(ii,jj,k+1) - zc(ii,jj,k))*DL)
   vsh(ii,jj,k) = (v(ii,jj,k+1,1) - v(ii,jj,k,1))*UL / ((zc(ii,jj,k+1) - zc(ii,jj,k))*DL)
   shq(ii,jj,k) = ush(ii,jj,k)*ush(ii,jj,k) + vsh(ii,jj,k)*vsh(ii,jj,k)
end do   

  ush(ii,jj,NK) = ush(ii,jj,NK-1)
  vsh(ii,jj,NK) = vsh(ii,jj,NK-1)

!  ush(ii,jj,NK) = stressx(jj)/(KzMom(ii,jj,NK) * 1027.d0 )
!  vsh(ii,jj,NK) = stressy(jj)/(KzMom(ii,jj,NK) * 1027.d0 )
  shq(ii,jj,NK) = ush(ii,jj,NK)*ush(ii,jj,NK) + vsh(ii,jj,NK)*vsh(ii,jj,NK)
!!!-----------------------------------------------------------------------------------------------------------
             
!! Calculate buoyancy frequency directly from the density
do k=1, (NK-1)   
  nnface(ii,jj,k) = (1/DL)*(-9.81d0/1027.d0) * (rho_old(ii,jj,k+1) - rho_old(ii,jj,k)) /  (zc(ii,jj,k+1) - zc(ii,jj,k))

if (k==1) then
  nnface(ii,jj,k-1) = nnface(ii,jj,k)
end if
end do

k=NK   
  nnface(ii,jj,k) = (-9.81d0/1027.d0)*( dtr0 * ( (1/KzTr(ii,jj,k))*( swr(jj) - qloss(jj) )/(1027.d0*4187.d0) ) )

do k=0,NK
  Rig(ii,jj,k) = nnface(ii,jj,k)/shq(ii,jj,k)
end do


return
end subroutine shearn2
