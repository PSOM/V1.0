subroutine wind_stress(udif,vdif,step) 
  !     ---------------------------------------------                     
  USE header

  implicit none 
  integer i,j,k,m,step
                                                                        
  REAL(kind=rc_kind) :: udif(NI,NJ,NK),vdif(NI,NJ,NK)

  REAL(kind=rc_kind) :: Kdudzt,Kdvdzt,fact,fac,rhoinv

  INTEGER :: yw, ymid, yset

  real :: stressmax, yyy

  udif=0.;vdif=0.;


!*********************************
! COMPUTATION OF THE WIND STRESS

  stress_top_x=0.;stress_top_y=0; 
  !stress_top_x=0.;stress_top_y=0.26; 
  stressmax= 0.26d0  
  
  !stressx(NJ/4:(3*NJ/4))= 0.02d0 
  !stress_top_x(:,NJ/4:(3*NJ/4))= 0.02d0 

  yw= 50.d0
  ymid= 0.5* dble(NJ)
  yset =30.d0

  do j=1,NJ/2
     yyy= dble(j)-0.5
     stress_top_y(1,j) = stressmax* 0.5* &
           & (1.d0 +tanh((( yyy-yset)/yw)*PI ))
  end do

  yset= NJ-30
  do j=NJ/2+1,NJ
    yyy= dble(j)-0.5
    stress_top_y(1,j) = stressmax*0.5* &
        & (1.d0 +tanh((-( yyy-yset)/yw)*PI ))
  end do

   do j=1,NJ
     if (abs(stress_top_y(1,j)-stressmax).lt.0.001) stress_top_y(1,j)=stressmax
     if (abs(stress_top_y(1,j)).lt.0.001) stress_top_y(1,j)=0.d0
   end do

   do i=1,NI
   do j=1,NJ
     stress_top_y(i,j)=stress_top_y(1,j)
   end do
   end do

   if (step>240) stress_top_y=0.d0

  !************************************
! COMPUTATION OF THE SOURCE TERM


  fac= 1.d0/(UL*DL*delta) 
  fact = DL/UL 

  do j=1,NJ 
    do i=1,NI 
      stress_top(i,j) = sqrt(stress_top_x(i,j)*stress_top_x(i,j)+ stress_top_y(i,j)*stress_top_y(i,j))
      !stress_top(i,j)=0
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
                                                                        
                                                                        
