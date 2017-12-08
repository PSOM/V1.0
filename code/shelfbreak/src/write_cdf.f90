  subroutine write_cdf(step,n)

  USE header,ONLY : NI,NJ,NK,ntr,nconsume,dirout,out1d_int,out2d_int,out3d_int,rc_kind,pvt,pv1,pv2,pv3
  IMPLICIT NONE

  INTEGER :: step,counter_2d,counter_3d,counter_1d,ksurf,islice,jslice,imooring,jmooring,n
  REAL :: sigma,z

! out1d_int  - frequency of 1d output
! out2d_int  - frequency of 2d output
! out3d_int  - frequency of 3d output

  ! 1D output
  if (mod(step,out1d_int).eq.0) then
    counter_1d= step/out1d_int +1 
imooring=NI/2;jmooring=NJ/2;    call write_cdf_1D_mooring(imooring,jmooring,counter_1d)
  end if

  ! 2D output
  if (mod(step,out2d_int).eq.0) then
    counter_2d= step/out2d_int +1 

    call diag_pv(n);pvt=pv1+pv2+pv3;

!   ksurf=NK/2;  call write_cdf_2D_sigma(ksurf,counter_2d,n)
 !   ksurf=9;  call write_cdf_2D_sigma(ksurf,counter_2d,n)
 !   ksurf=4;  call write_cdf_2D_sigma(ksurf,counter_2d,n)
!       ksurf=1;  call write_cdf_2D_sigma(ksurf,counter_2d,n)
!    ksurf=NK;    call write_cdf_2D_sigma(ksurf,counter_2d,n)

!   islice=  5;call write_cdf_2D_x(islice,counter_2d,n)
!   islice= 10;call write_cdf_2D_x(islice,counter_2d,n)
!   islice= 15;call write_cdf_2D_x(islice,counter_2d,n)
!   islice= 20;call write_cdf_2D_x(islice,counter_2d,n)
!   islice= 25;call write_cdf_2D_x(islice,counter_2d,n)
!   islice= 30;call write_cdf_2D_x(islice,counter_2d,n)
!   islice= 35;call write_cdf_2D_x(islice,counter_2d,n)
!   islice= 40;call write_cdf_2D_x(islice,counter_2d,n)
    islice= 48;call write_cdf_2D_x(islice,counter_2d,n)
!   islice= 60;call write_cdf_2D_x(islice,counter_2d,n)
!   islice= 70;call write_cdf_2D_x(islice,counter_2d,n)
!   islice= 80;call write_cdf_2D_x(islice,counter_2d,n)
!   islice= 90;call write_cdf_2D_x(islice,counter_2d,n)
 !  islice= 96;call write_cdf_2D_x(islice,counter_2d,n)
  jslice= NJ/2;call write_cdf_2D_y(jslice,counter_2d,n)
  jslice= 110;call write_cdf_2D_y(jslice,counter_2d,n)

    sigma=1025.00d0;  call write_cdf_2D_isopycnal(sigma,counter_2d,n) ! writes the solution on the sigma isopycnal
    sigma=1025.20d0;  call write_cdf_2D_isopycnal(sigma,counter_2d,n) ! writes the solution on the sigma isopycnal
    sigma=1025.60d0;  call write_cdf_2D_isopycnal(sigma,counter_2d,n) ! writes the solution on the sigma isopycnal


    z=-5.;    call write_cdf_2D_geopotential(z,counter_2d,n) ! writes the solution on the geopotential
    z=-25.;   call write_cdf_2D_geopotential(z,counter_2d,n) ! writes the solution on the geopotential
    z=-50.;   call write_cdf_2D_geopotential(z,counter_2d,n) ! writes the solution on the geopotential
    z=-75.;   call write_cdf_2D_geopotential(z,counter_2d,n) ! writes the solution on the geopotential
    z=-90.;   call write_cdf_2D_geopotential(z,counter_2d,n) ! writes the solution on the geopotential
 !   z=-100.;  call write_cdf_2D_geopotential(z,counter_2d,n) ! writes the solution on the geopotential
 !   z=-150.;  call write_cdf_2D_geopotential(z,counter_2d,n) ! writes the solution on the geopotential
 !   z=-200.;  call write_cdf_2D_geopotential(z,counter_2d,n) ! writes the solution on the geopotential
 !   z=-250.;  call write_cdf_2D_geopotential(z,counter_2d,n) ! writes the solution on the geopotential
 !   z=-300.;  call write_cdf_2D_geopotential(z,counter_2d,n) ! writes the solution on the geopotential

  end if

  ! 3D output
  if (mod(step,out3d_int).eq.0) then
    counter_3d= step/out3d_int +1 
   call write_cdf_3D(step,n)
!   call write_cdf_3D_strain(step,n)

  end if
     
! ksurf=NK;call write_cdf_2D_sigma(frame_int,step,ksurf,h,consump,Tr,s,T,rho,u,v,w,p,vor,strain,freqN2,xc,yc,zc,DL,LEN,Jac,dtf*TL)
! ksurf= INT(NK/3); call write_cdf_2D_sigma(frame_int,step,ksurf,h,consump,Tr,s,T,rho,u,v,w,p,vor,strain,freqN2,xc,yc,zc,DL,LEN,Jac,dtf*TL)
! ksurf= 1; call write_cdf_2D_sigma(frame_int,step,ksurf,h,consump,Tr,s,T,rho,u,v,w,p,vor,strain,freqN2,xc,yc,zc,DL,LEN,Jac,dtf*TL)

  END 
