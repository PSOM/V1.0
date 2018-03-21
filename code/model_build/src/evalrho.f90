subroutine evalrho(rhonew,nl) 
!     ---------------------------------------------                     
#include "cppdefs.h"
  USE header,only:NI,NJ,NK,rc_kind
  IMPLICIT NONE

  REAL(kind=rc_kind) :: rhonew(0:NI+1,0:NJ+1,0:NK+1)
  INTEGER :: nl

#ifdef rhoonly
  call evalrho_rho(rhonew,nl) 
#else
  call evalrho_sT(rhonew,nl) 
#endif

  return 
END subroutine evalrho
