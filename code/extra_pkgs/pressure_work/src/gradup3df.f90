subroutine gradup3df

USE header

implicit none

integer :: i,j,k

REAL(kind=rc_kind)  :: dufpdx, dvfpdy, dwfpdz
REAL(kind=rc_kind)  :: dufdx , dvfdy , dwfdz
REAL(kind=rc_kind)  :: udpdx , vdpdy , wdpdz

REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1)  :: ptc
REAL(kind=rc_kind), dimension(    0:NI,    NJ  ,   NK  )  :: ptfx
REAL(kind=rc_kind), dimension(      NI,  0:NJ,     NK  )  :: ptfy
REAL(kind=rc_kind), dimension(      NI,    NJ,   0:NK  )  :: ptfz

call diag_ptc  (ptc)
call intpol_pt (ptc,ptfx,ptfy,ptfz)

do k=1,NK
do j=1,NJ
do i=1,NI

  dufpdx= uf(i,j,k)*ptfx(i,j,k)-uf(i-1,j  ,k  )*ptfx(i-1,j  ,k  )
  dvfpdy= vf(i,j,k)*ptfy(i,j,k)-vf(i  ,j-1,k  )*ptfy(i  ,j-1,k  )
  dwfpdz= wf(i,j,k)*ptfz(i,j,k)-wf(i  ,j  ,k-1)*ptfz(i  ,j  ,k-1)

  !dum1(i,j,k)=(dufpdx+dvfpdy+dwfpdz)/Jac(i,j,k)/TL

  dufdx= (uf(i,j,k)-uf(i-1,j  ,k  )) 
  dvfdy= (vf(i,j,k)-vf(i  ,j-1,k  )) 
  dwfdz= (wf(i,j,k)-wf(i  ,j  ,k-1)) 

  !dum2(i,j,k)=(dufdx+dvfdy+dwfdz)/Jac(i,j,k)/TL

  udpdx=0.5*(uf(i,j,k)+uf(i-1,j  ,k  ))*(ptfx(i,j,k)-ptfx(i-1,j  ,k  ))
  vdpdy=0.5*(vf(i,j,k)+vf(i  ,j-1,k  ))*(ptfy(i,j,k)-ptfy(i  ,j-1,k  ))
  wdpdz=0.5*(wf(i,j,k)+wf(i  ,j  ,k-1))*(ptfz(i,j,k)-ptfz(i  ,j  ,k-1))

  !dum3(i,j,k)=(udpdx+vdpdy+wdpdz)/Jac(i,j,k)/TL
  
enddo
enddo
enddo

return


end
