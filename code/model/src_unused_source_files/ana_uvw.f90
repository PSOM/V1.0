SUBROUTINE ana_uvw(step)
  !     --------------------       
  !  specify idealized velocity field to test the particle module                                       
  USE header, ONLY: rc_kind,NI,NJ,NK,uf,vf,wf,dtf,pi
!  USE grids
  IMPLICIT NONE

  INTEGER i,j,step
  REAL(kind=rc_kind) ::  y,x,T0,l,r,L0
  INTEGER,PARAMETER :: seed = 86456  
  CALL RANDOM_SEED()

!  CALL ini_grids()

  ! the famous sine function
  ! ===========================
  ! cases 
  ! 1: the famous sine function
  !     u = (A sin(2pi k y /L),0) 0<=t<T/2
  !       = (0, A sin(2pi k x /L) T/2<=t<T
  ! 2:
  SELECT CASE (2)
  CASE(1)
     T0 = 0.25d0
     L0 = 1
     l = 1d0
     IF (MOD(step*dtf,T0)>=T0/2d0) THEN
        uf = 0d0
        DO i = 1, NJ
           CALL RANDOM_NUMBER(r)
           x = REAL(i)/REAL(NI)
           vf(i,:,:) = 1d-3*SIN(2d0*pi*l*x/L0) 
        ENDDO
     ELSE
        vf=0d0
        DO j = 1, NJ
           CALL RANDOM_NUMBER(r)
           y = REAL(j)/REAL(NJ)
           uf(:,j,:) = 1d-3*SIN(2d0*pi*l*y/L0)
        ENDDO
     ENDIF
     wf = 1d-3
  CASE(2)
     DO j = 0, NJ
        do i = 1, NI
           y = REAL(j)/REAL(NJ+1)
           x = REAL(i)/REAL(NI)
           vf(i,j,:) = 5d-3*COS(pi*x)*sin(pi*y)
        ENDDO
     ENDDO

     DO j = 1, NJ
        DO i=0, NI
           y = REAL(j)/REAL(NJ)
           x = REAL(i)/REAL(NI+1)
           uf(i,j,:) = -5d-3*sin(pi*x)*cos(pi*y)
        ENDDO
     enddo
     END SELECT


   END SUBROUTINE ana_uvw
