subroutine wind_stress(udif,vdif,step) 
  !     ---------------------------------------------                     
  USE header

  implicit none 
  integer i,j,k,m,step
                                                                        
  REAL(kind=rc_kind) :: udif(NI,NJ,NK),vdif(NI,NJ,NK)

  REAL(kind=rc_kind) :: Kdudzt,Kdvdzt,fact,fac,rhoinv

  real(kind=rc_kind) :: yw , yset, y0
  integer            :: iyw, iy0

  real :: stressmax

  udif=0.;vdif=0.;


!*********************************
! COMPUTATION OF THE WIND STRESS

  stress_top_x=0.;stress_top_y=0; 
  !stress_top_x=0.;stress_top_y=0.26; 
  stressmax= 0.d0  
  
  !************************************
! COMPUTATION OF THE SOURCE TERM


  fac= 1.d0/(UL*DL*delta) 
  fact = DL/UL 

  do j=1,NJ 
    do i=1,NI 
      stress_top(i,j) = sqrt(stress_top_x(i,j)*stress_top_x(i,j)+ stress_top_y(i,j)*stress_top_y(i,j))
      rhoinv = 1.d0/rho(i,j,NK) 
      !rhoinv = 1.d0/R0 
      Kdudzt= stress_top_x(i,j)*rhoinv*fact 
      Kdvdzt= stress_top_y(i,j)*rhoinv*fact 

      udif(i,j,NK)= fac*Jac(i,j,NK)*wz(i,j,NK)*Kdudzt   
      vdif(i,j,NK)= fac*Jac(i,j,NK)*wz(i,j,NK)*Kdvdzt   
    
    end do ! i
  end do ! j
                                                                        
return 
END                                           
                                                                        
                                                                        
