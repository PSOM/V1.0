  subroutine write_cdf(step,n)

  USE header,ONLY : NI,NJ,NK,ntr,nconsume,dirout,out1d_int,out2d_int,out3d_int,rc_kind,pvt,pv1,pv2,pv3
  IMPLICIT NONE

  INTEGER :: step,counter_2d,counter_3d,counter_1d,ksurf,islice,jslice,imooring,jmooring,n
  REAL :: sigma,z

! out1d_int  - frequency of 1d output
! out2d_int  - frequency of 2d output
! out3d_int  - frequency of 3d output

  ! 1D output
!  if (mod(step,out1d_int).eq.0) then
!    counter_1d= step/out1d_int +1 
!    imooring=NI/2;jmooring=NJ/2;    call write_cdf_1D_mooring(imooring,jmooring,counter_1d)
!  end if

  ! 2D output
!  if (mod(step,out2d_int).eq.0) then
!    counter_2d= step/out2d_int +1 
!    call diag_pv(n);pvt=pv1+pv2+pv3;
!    ksurf=NK;    call write_cdf_2D_sigma(ksurf,counter_2d,n)
!    islice=  5;call write_cdf_2D_x(islice,counter_2d,n)
!  end if

  ! 3D output
  if (mod(step,out3d_int).eq.0) then
    counter_3d= step/out3d_int +1 
   call write_cdf_3D(step,n)
!   call write_cdf_3D_strain(step,n)
  end if
     
  END 
