subroutine parti_forward()
  use header, only : NP,dtf
  use particles, only : parti
  integer :: i
print*, "dtf = ", dtf
do i = 1, NP
  parti(i)%i = parti(i)%i + 0.5d0 * dtf * (3d0 * parti(i)%u - parti(i)%u0)
  parti(i)%j = parti(i)%j + 0.5d0 * dtf * (3d0 * parti(i)%v - parti(i)%v0)
  parti(i)%k = parti(i)%k + 0.5d0 * dtf * (3d0 * parti(i)%w - parti(i)%w0)

  parti(i)%u0 = parti(i)%u
  parti(i)%v0 = parti(i)%v
  parti(i)%w0 = parti(i)%w
enddo

end subroutine parti_forward
