subroutine sponge(n,dum_time)

  USE header

  implicit none

  REAL(kind=rc_kind) ::  dum_time, gamma_T, dummy
  REAL(kind=rc_kind) ::  rinn, rout, rrr, rdum, s_relax, u_relax, v_relax
  INTEGER, PARAMETER ::  spl=15, splz=60
  
  INTEGER :: i,j,k,n

  gamma_T = 1.d0 / (dum_time)

  do j=0,spl+1

    dummy=gamma_T*cos((j-1)/spl*PI/2) 
    if (j==0) dummy=gamma_T

    do k = splz+1, NK+1
    do i = 0 , NI+1
      s(i,j,k,n) = s(i,j,k,n) - dum_time * dummy * (s(i,j,k,n) - rho_refS(k))
      u(i,j,k,n) = u(i,j,k,n) - dum_time * dummy *  u(i,j,k,n)
      v(i,j,k,n) = v(i,j,k,n) - dum_time * dummy *  v(i,j,k,n)
      w(i,j,k,n) = w(i,j,k,n) - dum_time * dummy *  w(i,j,k,n)
    enddo
    enddo

  enddo

  do j=NJ-spl,NJ+1

    dummy=gamma_T*sin((j-(NJ-spl))/spl*PI/2) 
    if (j==NJ+1) dummy=gamma_T

    do k = splz+1, NK+1
    do i = 0 , NI+1
      s(i,j,k,n) = s(i,j,k,n) - dum_time * dummy * (s(i,j,k,n) - rho_refN(k))
      u(i,j,k,n) = u(i,j,k,n) - dum_time * dummy *  u(i,j,k,n)
      v(i,j,k,n) = v(i,j,k,n) - dum_time * dummy *  v(i,j,k,n)
      w(i,j,k,n) = w(i,j,k,n) - dum_time * dummy *  w(i,j,k,n)
    enddo
    enddo

  enddo

  !rinn=dum_time*240.d0*6.d0
  !print*, ffc(NI/2,NJ/2)*FPAR, 1/(ffc(NI/2,NJ/2)*FPAR), 2.d0*PI/(ffc(NI/2,NJ/2)*FPAR)
  rinn=2.d0*PI/(ffc(NI/2,NJ/2)*FPAR)
  rout=rinn/1.d3

  do k=0,splz

    rdum =(dble(splz-k)*rout+dble(k)*rinn)/dble(splz)
    if (rdum .ne. 0) then
      rrr=1.d0/rdum
    else
      rrr=0.d0
    endif

    do j = 0, NJ+1
    do i = 0, NI+1

      s_relax=( dble(splz-k)*rho_refB(j)+dble(k)*s(i,j,k,n) )/dble(splz)
      u_relax=(                          dble(k)*u(i,j,k,n) )/dble(splz)
      v_relax=(                          dble(k)*v(i,j,k,n) )/dble(splz)

      s(i,j,k,n) = s(i,j,k,n) - dum_time*rrr*(s(i,j,k,n)-s_relax)
      u(i,j,k,n) = u(i,j,k,n) - dum_time*rrr*(u(i,j,k,n)-u_relax)
      v(i,j,k,n) = v(i,j,k,n) - dum_time*rrr*(v(i,j,k,n)-v_relax)

    !  w(i,j,k,n) = w(i,j,k,n) - dum_time * dummy *  w(i,j,k,n)
    enddo
    enddo

  enddo

end subroutine sponge  
