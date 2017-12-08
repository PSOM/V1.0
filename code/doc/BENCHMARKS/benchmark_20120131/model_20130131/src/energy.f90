subroutine energy(step) 
!---------------------------------------------------                    
  USE header
                                                                       
  integer :: step 
  REAL(kind=rc_kind) ::  const,enkin,enkin2
                                                                   
  const= EPS*EPS*delta*delta 

  enkin = SUM( Jac(1:NI,1:NJ,1:NK) * ( u(1:NI,1:NJ,1:NK,0)**2 + v(1:NI,1:NJ,1:NK,0)**2 + const*w(1:NI,1:NJ,1:NK,0)**2  ))
              
  print*, "#total kinetic energy = ", step, enkin
  return 
END subroutine energy                                       
