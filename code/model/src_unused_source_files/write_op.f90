SUBROUTINE write_op(step)
#include "cppdefs.h"
  USE header
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: step
  CHARACTER(LEN=10) :: stepchar 



! ----------------------------
! Output every ngraph2d steps

if (mod(step,ngraph2d) .eq. 0) then

#ifdef file_output_bin 
  call save2d(NI+2,NJ+2,rho(:,:,NK),TRIM(dirout)//'op.rho.k01.'//stepchar(step)//'.bin')
  call save2d(NI+2,NJ+2,h,TRIM(dirout)//'op.h.'//stepchar(step)//'.bin')
#endif

#ifdef file_output_cdf 
ksurf=NK;call writeksurf(frame_int,step,ksurf,h,consump,Tr,s,T,rho,u,v,w,p,vor,strain,freqN2,xc,yc,zc,DL,LEN,Jac,dtf*TL)
ksurf= INT(NK/3); call writeksurf(frame_int,step,ksurf,h,consump,Tr,s,T,rho,u,v,w,p,vor,strain,freqN2,xc,yc,zc,DL,LEN,Jac,dtf*TL)
#endif

endif
! ----------------------------



! ----------------------------
! Output every ngraph2d steps
!   save the whole model state every save_steps steps

if (mod(step,save_steps) .eq. 0) then

#ifdef file_output_bin 
  call save3d(NI+2,NJ+2,NK+2,u(:,:,:,1),TRIM(dirout)//'op.u.'//stepchar(step)//'.bin')
  call save3d(NI+2,NJ+2,NK+2,v(:,:,:,1),TRIM(dirout)//'op.v.'//stepchar(step)//'.bin')
  call save3d(NI+2,NJ+2,NK+2,w(:,:,:,1),TRIM(dirout)//'op.w.'//stepchar(step)//'.bin')
  call save3d(NI+2,NJ+2,NK+2,rho,TRIM(dirout)//'op.rho.'//stepchar(step)//'.bin')
  CALL vort(0)
  call save3d(NI+2,NJ+2,NK+2,vor,TRIM(dirout)//'op.vor.'//stepchar(step)//'.bin')

#endif

endif
! ----------------------------



! ----------------------------
! Output every pickup_int
!  save a pickup file

if (mod(step,pickup_int) .eq. 0) then

#ifdef file_output_bin 
 call w_pickup(TRIM(dirout)//'op.pickup.'//stepchar(step)//'.bin')
#endif

endif
! ----------------------------





END SUBROUTINE
