SUBROUTINE ini_rho()
  !     --------------------                                              
  USE header, ONLY: NI,NJ,NK,rho,s,T,pi,zc,DL,LEN,yc,xc,dx,dy,total_depth,T_ref
!  USE grids
  IMPLICIT NONE

  ! Initializes s,T from Andrey's section gridded for model grid
  ! s_init.dat and T_init.dat  written in ascii - NK+2 x NJx2

  INTEGER n,i,j,k,Nu,k0
  REAL(kind=rc_kind) ::  A,B,C,D,G,potdens,x,y(0:NJ+1),z,tmp,yfront,mldepth,width,slfac
  REAL(kind=rc_kind) ::  tmp1, tmp2, tmp3, tmp4
  INTEGER,PARAMETER :: seed = 86456  
  print*, 'ini_grids'
  CALL RANDOM_SEED()

  ! assign the density or temperature and salinity by either analytic functions or 
  ! any particular hydrographic section.
!  CALL ini_grids()

  ! t=A/B exp(A*z)erf(B*y) + G*(1+z)
  ! z=-1:0, y = -0.5,0.5, 

  ! ===========================
  ! cases 
  !-2: eddy
  ! 0: jet
  ! 1: idealized surface jet
  ! 2: krushio
  SELECT CASE (1)
  CASE(1)
     A = 4.0d0
     B = 20d0 !jet width in km as yc is in km
     C = 2d0
     G = 19d0
     D = 0.5d0  !linear temperature difference from south to north 
     !Nu =NK-9
     mldepth = 100d0
     s(:,:,:,0) = 34d0
     y = yc -yc(0) - (yc(NJ+1)-yc(0))/2d0
     print*, "y=",y
     print*, "yc=",yc
     DO k = NK, 0, -1
        DO j = 0, NJ+1
           DO i = 0, NI+1
              CALL RANDOM_NUMBER(tmp)
              !y = -(REAL(j)-REAL(NJ+1)/2d0)/REAL(NJ+1)
              !y = -(dble(j)-dble(NJ+1)/2d0)*dy
              ! y = - (yc(j) - (yc(NJ+1)-yc(0))/2d0)
              !IF (k>Nu) THEN
              IF (zc(1,1,k)*DL>-mldepth) THEN
                 z=0d0
                 Nu=k
              ELSE
                 z = dble(k-Nu)/dble(Nu)
              ENDIF
              if (k<=1) then
                  T(i,j,k,0) = C * EXP(A*z) + G*(1+z) + 0.02*tmp ! + 1d-2*SIN(DBLE(i)/DBLE(NI+1)*2d0*pi*1d0)
              else
              T(i,j,k,0) = C * EXP(A*z) * (erf(-y(j)/B) )  + G*(1+z) + 0.02*tmp ! + 1d-2*SIN(DBLE(i)/DBLE(NI+1)*2d0*pi*1d0)
              endif
              !T(i,j,k,0) = C * EXP(A*z) * erf(B*y) + G*(1+z) + 1d-3*SIN(DBLE(i)/DBLE(NI+1)*2d0*pi*3d0)
              rho(i,j,k) = potdens(s(0,j,k,0),T(0,j,k,0))
              !print*, "t=",t(1,j,k,0),"rho=",rho(1,j,k)
           ENDDO
        ENDDO
     ENDDO
     T(:,0,:,0)= T(:,1,:,0)
     T(:,NJ+1,:,0)= T(:,NJ,:,0)
     s(:,0,:,0)= s(:,1,:,0)
     s(:,NJ+1,:,0)= s(:,NJ,:,0)
     rho(:,0,:)= rho(:,1,:)
     rho(:,NJ+1,:)= rho(:,NJ,:)
     T(:,:,NK+1,:)=T(:,:,NK,:)
     s(:,:,NK+1,:)=s(:,:,NK,:)
     rho(:,:,NK+1)=rho(:,:,NK)
     DO k = 1, NK
        rho(:,:,k) = 0.25d0*(rho(:,:,k-1)+2d0*rho(:,:,k)+rho(:,:,k+1))
        s(:,:,k,0) = 0.25d0*(s(:,:,k-1,0)+2d0*s(:,:,k,0)+s(:,:,k+1,0))
        T(:,:,k,0) = 0.25d0*(T(:,:,k-1,0)+2d0*T(:,:,k,0)+T(:,:,k+1,0))
     ENDDO
  case(-4)
     s(:,:,:,:) = 34d0
     DO k = 0, NK+1
        z = zc(1,1,k)*DL
        T(:,:,k,:) = 28d0 + 5d0 * z/800d0
        DO i = 0, NI+1
           DO j = 0, NJ+1
              rho(i,j,k) = potdens(s(i,j,k,0),T(i,j,k,0))
           enddo
        enddo
      enddo

  case(-3)
     call stprofile
  end select
  T_ref = T(:,:,:,0)

END SUBROUTINE ini_rho
