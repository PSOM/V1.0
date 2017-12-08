subroutine diag_ptc(ptc)

USE header

implicit none

integer i,j,k
real(kind=rc_kind) :: const

REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1)  :: ptc

const=10.d0*gpr

call rpevalgrad_Song(0)

do k=1, NK
  do j=1, NJ
  do i=1, NI
    ptc(i,j,k)=qpr*P1*p (i,j,k)+rp(i,j,k)*P1+h (i,j)*gpr*10.*R0
  end do
  end do
  ptc(     0,1:NJ,k)=ptc(     1,1:NJ,k)
  ptc(  NI+1,1:NJ,k)=ptc(    NI,1:NJ,k)
  ptc(0:NI+1,   0,k)=ptc(0:NI+1,   1,k)
  ptc(0:NI+1,NJ+1,k)=ptc(0:NI+1,  NJ,k)
end do
ptc(0:NI+1,0:NJ+1,   0)=ptc(0:NI+1,0:NJ+1, 1)
ptc(0:NI+1,0:NJ+1,NK+1)=ptc(0:NI+1,0:NJ+1,NK)

return
end
