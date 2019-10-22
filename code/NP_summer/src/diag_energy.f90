subroutine diag_energy(step) 
!---------------------------------------------------                    
  USE header
                                                                       
  integer :: step 
  REAL(kind=rc_kind) ::  const,enkin,enkin2
  logical :: exist

  const= EPS*EPS*delta*delta 

  enkin = SUM( Jac(1:NI,1:NJ,1:NK) * ( u(1:NI,1:NJ,1:NK,0)**2 + v(1:NI,1:NJ,1:NK,0)**2 + const*w(1:NI,1:NJ,1:NK,0)**2  ))


  ! Print the TKE into a file TKE.out
  inquire(file=TRIM(dirout)//'TKE.out', exist=exist)
  if (exist) then
    open(12, file=TRIM(dirout)//'TKE.out', status="old", position="append", action="write")
  else
     open(12, file=TRIM(dirout)//'TKE.out', status="new", action="write")
  end if
  WRITE(12,*) step,", ", enkin
  close(12)

  ! Print the TKE to the command line
  print*, "#total kinetic energy = ", step, enkin
  return 
END subroutine diag_energy                                       
