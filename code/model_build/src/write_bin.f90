SUBROUTINE write_bin(step)
#include "cppdefs.h"
  USE header
  IMPLICIT NONE

  INTEGER :: step
  CHARACTER(LEN=10) :: stepchar 
  


!-------
! OUTPUT


if (mod(step,out2d_int) .eq. 0) then
  call save2d(NI+2,NJ+2,rho(:,:,NK),TRIM(dirout)//'op.rho.k01.'//stepchar(step)//'.bin')
  call save2d(NI+2,NJ+2,h,TRIM(dirout)//'op.h.'//stepchar(step)//'.bin')
endif


if (mod(step,out3d_int) .eq. 0) then
  call save3d(NI+2,NJ+2,NK+2,u(:,:,:,1),TRIM(dirout)//'op.u.'//stepchar(step)//'.bin')
  call save3d(NI+2,NJ+2,NK+2,v(:,:,:,1),TRIM(dirout)//'op.v.'//stepchar(step)//'.bin')
  call save3d(NI+2,NJ+2,NK+2,w(:,:,:,1),TRIM(dirout)//'op.w.'//stepchar(step)//'.bin')
  call save3d(NI+2,NJ+2,NK+2,rho,TRIM(dirout)//'op.rho.'//stepchar(step)//'.bin')
  CALL vort(0)
  call save3d(NI+2,NJ+2,NK+2,vor,TRIM(dirout)//'op.vor.'//stepchar(step)//'.bin')
endif




!-------
! PICKUP

!Reading pickup files, only once, when step is equal to pickup_step
IF(pickup_step>-0.5) then
  if(step==pickup_step) then
    call r_pickup(step)
  endif
ENDIF


! Writing pickup files at some frequency
if(mod(step,pickup_int) .eq. 0) then
  call w_pickup(TRIM(dirout)//'op.pickup.'//stepchar(step)//'.bin')
endif

END SUBROUTINE
