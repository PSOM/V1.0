module relaxation

  use header, only : NI,NJ,NK,UL,LEN,T,T_ref,r_sponge,rc_kind
!  use grids, only : et_c, yc

 REAL(kind=rc_kind) :: r_T(0:NJ+1,0:NK+1)

contains

subroutine set_coef()
 implicit none
 integer :: i,j,k
 REAL(kind=rc_kind) ::  gamma_T,tmp(0:NJ+1)
! specific restoring profile

 gamma_T = 1d0 / (5d0 * 86400d0 * UL/LEN)
 do k = 0, NK+1
    r_T(:,k) = gamma_T * r_sponge
    !r_T(:,k) = gamma_T 
 enddo

end subroutine set_coef

subroutine sponge(n,dtime)
 implicit none

 REAL(kind=rc_kind) ::  varbar(0:NJ+1,0:NK+1),dtime
 integer :: i,j,k,n

! calculate the bar{var} in x

 do k = 0, NK+1
    do j = 0, NJ+1
       varbar(j,k) = sum(T(1:NI,j,k,n))/real(NI)
    enddo
 enddo

 do j = 0, NJ+1
    do k = 0, NK+1
       do i = 0, NI+1
       !T(i,j,k,n) = T(i,j,k,n) - dtime * r_T(j,k) * (varbar(j,k) - T_ref(i,j,k))
       T(i,j,k,n) = T(i,j,k,n) - dtime * r_T(j,k) * (T(i,j,k,n) - T_ref(i,j,k))
       enddo
    enddo
 enddo


end subroutine sponge


end module relaxation
