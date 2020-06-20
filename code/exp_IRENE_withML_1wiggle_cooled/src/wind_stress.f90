subroutine wind_stress(udif,vdif,step) 
  !     ---------------------------------------------                     
  USE header

  IMPLICIT NONE

  INTEGER           , INTENT(in) :: step
  REAL(kind=rc_kind), INTENT(out) :: udif(NI,NJ,NK),vdif(NI,NJ,NK)

  INTEGER :: i,j,k
  INTEGER :: iw0, iw1
  INTEGER, PARAMETER :: npw=825

  REAL(kind=rc_kind) :: Kdudzt,Kdvdzt,fact,fac,rhoinv
  REAL(kind=rc_kind) :: stressmax, xstressmax, ystressmax 
  REAL(kind=rc_kind) :: ttt, iw, WT0, WT1
  REAL(kind=rc_kind) :: yw, yset, yyy
!  REAL(kind=rc_kind) :: stress(2,npw)
  REAL(kind=rc_kind), PARAMETER :: dtw=6.944444e-03

!*********************************
! COMPUTATION OF THE WIND STRESS

  !---------initialize arrays and vars
  udif(:,:,:)=0.d0; vdif(:,:,:)=0.d0;

  stress_top_x=0.d0;stress_top_y=0.d0; 
  stressmax=0.d0  

!  read in ini_st to avoid the loop three times every timestep
!  IF ( step .EQ. 1 ) THEN
!    OPEN(UNIT=40,FILE='../windsALBOREX.in')
!    DO i=1,npw
!      READ(40,*) iw, stress(1,i),stress(2,i)
!    ENDDO
!    CLOSE(40)
!  ENDIF

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
                                                                        
                                                                        
