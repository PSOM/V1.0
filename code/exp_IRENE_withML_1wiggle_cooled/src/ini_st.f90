subroutine ini_st

!  initializes TS profiles from ALBOREX data

  USE header

  IMPLICIT NONE
  INTEGER, PARAMETER                 :: nwork =50000
  INTEGER i,j,k,iseed
  real(kind=rc_kind) :: sd,am,da(NI)
  real(kind=rc_kind) :: dwork(nwork)

!----------------------------------------------------

  ! read data
  OPEN(UNIT=40,FILE='exp_IRENE_withML_1wiggle_cooled/init_T.in')
  OPEN(UNIT=41,FILE='exp_IRENE_withML_1wiggle_cooled/init_s.in')

  DO k=0,NK+1
    DO j=0,NJ+1
		DO i=0,NI+1
      	  READ(40,*) T(i,j,k,0)
      	  READ(41,*) s(i,j,k,0)
	  ENDDO
    ENDDO
  ENDDO

  CLOSE(40)
  CLOSE(41)

  ! copy vertical section along x
  !DO k=0,NK+1
 ! 	DO j=0,NJ+1
 ! 	  DO i=0,NI+1
 !   	  T(i,j,k,0)=T(0,j,k,0) 
 !   	  s(i,j,k,0)=s(0,j,k,0)
 ! 	  ENDDO  
 ! 	ENDDO  
 ! ENDDO  
  

  ! read input winds
!  OPEN(UNIT=40,FILE='./ALBOREX_WMOP_gld_winds_consump/windsALBOREX.in')
!  DO i=1,825
!    READ(40,*) dum, stress(1,i),stress(2,i)
!  ENDDO
!  CLOSE(40)


  RETURN

END SUBROUTINE ini_st


  
