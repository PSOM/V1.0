function stepchar(step)
    implicit none
    integer, intent(in) :: step
    character(len=10) :: stepchar
    write(stepchar,'(I10.10)') step

end
