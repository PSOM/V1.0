  subroutine write_cdf(step,n)

  USE header,ONLY : NI,NJ,NK,ntr,nconsume,nsteps,dirout,out1d_int,out2d_int,out3d_int,rc_kind,pvt,pv1,pv2,pv3,time_days
  IMPLICIT NONE

  INTEGER :: step,counter_1d,counter_2d,counter_3d,ksurf,islice,jslice,imooring,jmooring,n
  REAL :: sigma,z
!   REAL(kind=rc_kind) :: ptc(0:NI+1,0:NJ+1,0:NK+1)
  
! out1d_int - frequency of 1d output
! out2d_int - frequency of 2d output
! out3d_int - frequency of 3d output

! -------------------  1D output  -------------------
!    if (mod(step,out1d_int).eq.0) then
! !      write(6,*) 'write_cdf_1D_mooring'
!      counter_1d= step/out1d_int +1
!      imooring=NI/2;jmooring=NJ/2; call write_cdf_1D_mooring(imooring,jmooring,counter_1d)
!    end if

! ! -------------------  2D output  -------------------
!   if (mod(step,out2d_int).eq.0) then
! !     write(6,*) 'write_cdf_2D'
!     counter_2d= step/out2d_int +1
!
! !     call diag_pv(n);pvt=pv1+pv2+pv3;
! !     call diag_ptc(ptc)
!
! !     write(6,*) 'write_cdf_2D_x'
!     islice= NI/2; call write_cdf_2D_x(islice,counter_2d,n)
!
! !     write(6,*) 'write_cdf_2D_y'
!     jslice= NJ/2; call write_cdf_2D_y(jslice,counter_2d,n)
!
! !   write(6,*) 'write_cdf_2D_sigma'
! !   ksurf=NK/2;  call write_cdf_2D_sigma(ksurf,counter_2d,n)
! !   call write_cdf_2D_sigma(frame_int,step,ksurf,h,consump,Tr, &
! !   & s,T,rho,u,v,w,p,vor,strain,freqN2,xc,yc,zc,DL,LEN,Jac,dtf*TL)
!
! !   write(6,*) 'write_cdf_2D_x_face'
! !   islice= NI/2; call write_cdf_2D_x_face(islice,counter_2d)
!
! !   write(6,*) 'write_cdf_2D_isopycnal' ! writes the solution on the sigma isopycnal
! !   sigma=25.d0;  call write_cdf_2D_isopycnal(sigma,counter_2d,n)
!
! !   write(6,*) 'write_cdf_2D_geopotential' ! writes the solution on the geopotential
! !   z=-300.;  call write_cdf_2D_geopotential(z,counter_2d,n)
!
!   end if

! -------------------  3D output  -------------------

if ((time_days<5).OR.(time_days>20)) then ! Day<7 or Day>14
!    if (mod(step,out3d_int).eq.0) then ! output every day
!       counter_3d= step/out3d_int +1
!       call write_cdf_3D(step,n)
!    endif
   
else ! more frequent
  if (mod(step,out2d_int).eq.0) then
    counter_3d= step/out3d_int +1
    call write_cdf_3D(step,n)
  endif
endif

  END 
