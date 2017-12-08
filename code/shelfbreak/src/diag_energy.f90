subroutine diag_energy(step) 
!---------------------------------------------------                    
  USE header
           
                                                            
  integer :: step 
  REAL(kind=rc_kind) ::  const,enkin,enkin2

  !-----------------
  ! N2
  REAL :: dn2

  !-----------------
  ! psi_e
  REAL :: vpbp, wpbp, by, bz
  REAL :: psi_e(1:NJ,1:NK)
  REAL, PARAMETER :: peps=1e-3

  !-----------------
  ! TRACER
  integer :: i,j,k,iy,iz
  integer, PARAMETER :: ntrby=5,ntrbz=3
  REAL :: d1, d2, d3, d4  
  REAL :: ra_trb(ntrby,ntrbz), ra_vol(ntrby,ntrbz), ra_trp(ntrby,ntrbz)
  REAL :: zmin,zmax,ymin,ymax
                                                                 
  const= EPS*EPS*delta*delta 

  enkin = SUM( Jac(1:NI,1:NJ,1:NK) * ( u(1:NI,1:NJ,1:NK,0)**2 + v(1:NI,1:NJ,1:NK,0)**2 + const*w(1:NI,1:NJ,1:NK,0)**2  ))
              
  print*, "#total kinetic energy = ", step, enkin


  ! user2=1e-5

!---------------------------
! psi_e computation
  DO j=1,NI
    DO k=1,NK
 
  vpbp=SUM(v(:,j,k,0)*(rho(:,j,k)-SUM(rho(:,j,k)/NI))/NI)
  wpbp=SUM(w(:,j,k,0)*(rho(:,j,k)-SUM(rho(:,j,k)/NI))/NI)
  by=SUM(freqby(:,j,k))/NI
  bz=SUM(freqbz(:,j,k))/NI

     psi_e(j,k)=peps*( ( ( peps * vpbp * by ) - ( 1./peps * wpbp * by ) ) / ( by**2 + peps**2 * bz**2 ) )
      IF(j==85 .AND. k==16) then
        PRINT*,"P ",psi_e(j,k),vpbp,wpbp,by,bz
      ENDIF
    ENDDO
  ENDDO

IF(.TRUE.) then


!---------------------------
! TRACER diagnosis

IF(step==1 .OR. MOD(step,50)==0) then

  d1=SUM( Jac(1:NI,1:NJ,1:NK) * Tr(1,1:NI,1:NJ,1:NK,0))

  ra_trb=0.; ra_trp=0.; ra_vol=0.;
  DO iy=1,5
    DO iz=1,3
     IF(iy==1) then;ymin=-500.*1e3;ymax=97.5*1e3; ENDIF;
     IF(iy==2) then;ymin=97.5*1e3; ymax=103.5*1e3;ENDIF;
     IF(iy==3) then;ymin=103.5*1e3;ymax=109.5*1e3;ENDIF;
     IF(iy==4) then;ymin=109.5*1e3;ymax=115.5*1e3;ENDIF;
     IF(iy==5) then;ymin=115.5*1e3;ymax=500.*1e3;ENDIF;
     IF(iz==1) then;zmin=-5000.;zmax=-80.;ENDIF;
     IF(iz==2) then;zmin=-80.;zmax=-40.;ENDIF;
     IF(iz==3) then;zmin=-40.;zmax=5000.;ENDIF;

  DO i=1,NI
   DO j=1,NJ
     DO k=1,NK
       IF(zc(i,j,k)*DL<zmax .AND. zc(i,j,k)*DL>zmin .AND. yc(j)*1e3<ymax .AND. yc(j)*1e3>ymin) then
         ra_trb(iy,iz) = ra_trb(iy,iz) + ( Jac(i,j,k) * Tr(1,i,j,k,0) )
         ra_vol(iy,iz) = ra_vol(iy,iz) + ( Jac(i,j,k) * 1. )
       ENDIF 
     ENDDO
   ENDDO
  ENDDO

  ra_trp(iy,iz)=ra_trb(iy,iz)/ra_vol(iy,iz)

  ENDDO
  ENDDO



  PRINT*,"D TR ------------------------------"
  PRINT*,"D TR",step, SUM(ra_trb(:,:))/d1
  PRINT*,"D TR"," PERCENTAGE OF OCCUPATION"
  PRINT'(A94)'," D TR step      < 98              < 104             < 110             < 116             >= 116" 
  PRINT'(A6,I4.4,5(E18.11))'," D TR ",step, ra_trp(1:ntrby,3)
  PRINT'(A6,I4.4,5(E18.11))'," D TR ",step, ra_trp(1:ntrby,2)
  PRINT'(A6,I4.4,5(E18.11))'," D TR ",step, ra_trp(1:ntrby,1)
  PRINT*,"D TR ------------------------------"
  PRINT*,"D TR"," PERCENTAGE OF TRACER"
  PRINT'(A94)'," D TR step      < 98              < 104             < 110             < 116             >= 116" 
  PRINT'(A6,I4.4,5(E18.11))'," D TR ",step, ra_trb(1:ntrby,3)/SUM(ra_trb(:,:))
  PRINT'(A6,I4.4,5(E18.11))'," D TR ",step, ra_trb(1:ntrby,2)/SUM(ra_trb(:,:))
  PRINT'(A6,I4.4,5(E18.11))'," D TR ",step, ra_trb(1:ntrby,1)/SUM(ra_trb(:,:))
  PRINT*,"D TR ------------------------------"


ENDIF


ENDIF


  return 
END subroutine
