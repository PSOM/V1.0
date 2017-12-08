       subroutine rpevalgrad_Sche(nl)

!==========================================
! Adaptation of the baroclinic pressure scheme from ROMS.
! April 2012, JBG.

USE header
     implicit none
      integer istr,istrU,iend,jstr,jstrV,jend, i,j,k, nl, N, imin,imax,jmin,jmax
      real, dimension(0:NI+1,0:NJ+1,NK) :: P_Sche
      real, dimension(0:NI+1,0:NK) :: dR,dZ 
      real, dimension(0:NI+1,0:NJ+1) :: FC,dZx,rx,dRx,dn_u,dm_v
      REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1,  0:NK+1)      :: rhol 
      REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1,  0:NK+1)      :: z_r,Hz 
      REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, -1:NK+1)      :: z_w 
      real GRho, HalfGRho, cff, cfr,g, rho0
      real, parameter :: OneFifth=0.2, OneTwelfth=1./12., epsil=0.
      real :: adim_Sche 


!====================================================================================h

! -----------------------------------------------------------------
! Part I: The variables of PSOM are translated into ROMS variables.
! -----------------------------------------------------------------


! --------
! Indices


imin=1; imax=NI;
jmin=1; jmax=NJ;
N=NK;

istrU=imin; iend=imax;
jstrV=jmin; jend=jmax;
istr=imin;
jstr=jmin;

! --------
! Geometry

dn_u(:,:)=dx; dm_v(:,:)=dy;

z_r=DL*zc; z_w=DL*zf; ! The scheme uses dimensional variables.

do k=1,NK
  Hz(:,:,k)=(z_w(:,:,k)-z_w(:,:,k-1))
enddo


! ----------
! Geophysics

g=gpr*10; rho0=R0;
call evalrho(rho,nl);rhol=rho-rho0;


! PRINT*,zc(NI/2,NJ/2,:)
! PRINT*,zf(NI/2,NJ/2,:)
! PRINT*,Hz(NI/2,NJ/2,:)
! PRINT*,dx,dy
! PRINT*,rho(NI/2,NJ/2,:)+R0

! ------------------------
! Boundary condition in x:

Hz(0,:,:)=Hz(NI,:,:); Hz(NI+1,:,:)=Hz(1,:,:);
rhol(0,:,:)=rhol(NI,:,:); rhol(NI+1,:,:)=rhol(1,:,:);


!====================================================================================h

! -----------------------------------------------------------------
! Part II: The code.
! -----------------------------------------------------------------
! The cpp tests have been removes.
! Except from explicitely designed ligns, nothing has been changed.

!
! Preliminary step (same for XI- and ETA-components):
!------------ ---- ----- --- --- --- ----------------
!
      GRho=g/rho0
      HalfGRho=0.5*GRho
 
      do j=jstrV-1,jend
        do k=1,N-1
          do i=istrU-1,iend
            dZ(i,k)=z_r(i,j,k+1)-z_r(i,j,k)
            dR(i,k)=rhol(i,j,k+1)-rhol(i,j,k)
          enddo
        enddo
        do i=istrU-1,iend
          dR(i,N)=dR(i,N-1)
          dR(i,0)=dR(i,1)
          dZ(i,N)=dZ(i,N-1)
          dZ(i,0)=dZ(i,1)
        enddo
        do k=N,1,-1               !--> irreversible
          do i=istrU-1,iend
            cff=2.*dZ(i,k)*dZ(i,k-1)
            dZ(i,k)=cff/(dZ(i,k)+dZ(i,k-1))
 
            cfr=2.*dR(i,k)*dR(i,k-1)
            if (cfr.gt.epsil) then
              dR(i,k)=cfr/(dR(i,k)+dR(i,k-1))
            else
              dR(i,k)=0.
            endif
          enddo
        enddo
        do i=istrU-1,iend    ! Eq (P.1)
          P_Sche(i,j,N)=0*g*z_w(i,j,N) + GRho*( rhol(i,j,N)                     &
     &       +0.5*(rhol(i,j,N)-rhol(i,j,N-1))*(z_w(i,j,N)-z_r(i,j,N))     &
     &          /(z_r(i,j,N)-z_r(i,j,N-1)) )*(z_w(i,j,N)-z_r(i,j,N))
        enddo
        do k=N-1,1,-1
          do i=istrU-1,iend   ! Eq (P.2)
            P_Sche(i,j,k)=P_Sche(i,j,k+1)+HalfGRho*( (rhol(i,j,k+1)+rhol(i,j,k))    &
     &                                     *(z_r(i,j,k+1)-z_r(i,j,k))   &
 
     &     -OneFifth*( (dR(i,k+1)-dR(i,k))*( z_r(i,j,k+1)-z_r(i,j,k)    &
     &                              -OneTwelfth*(dZ(i,k+1)+dZ(i,k)) )   &
 
     &                -(dZ(i,k+1)-dZ(i,k))*( rhol(i,j,k+1)-rhol(i,j,k)    &
     &                              -OneTwelfth*(dR(i,k+1)+dR(i,k)) )   &
     &                                                             ))   
          enddo
        enddo
      enddo   !<-- j


! Caution:
! The line P_Sche(i,j,N)=... has been modified !
! The term g*z_w(i,j,N) has been set to 0 because PSOM only requires the baroclinic terms and not 
!  the barotropic part linked with surface elevation (g.R0.h).
!


!
! Compute XI-component of pressure gradient term:
!-------- ------------ -- -------- -------- -----
!
      do k=N,1,-1
        do j=jstr,jend
          do i=imin,imax
            FC(i,j)=(z_r(i,j,k)-z_r(i-1,j,k))
            rx(i,j)=(rhol(i,j,k)-rhol(i-1,j,k))

          enddo
        enddo


! This part has been included to take into account the boundary conditions:
!++++ Added by JBG, 20120425.            
       FC(0,:)=FC(NI,:);FC(1,:)=(z_r(1,:,k)-z_r(NI,:,k))
       rx(0,:)=rx(NI,:);rx(1,:)=(rhol(1,:,k)-rhol(NI,:,k))
 
       FC(NI+1,:)=FC(1,:) 
       rx(NI+1,:)=rx(1,:) 

 
!# ifndef EW_PERIODIC
!        if (WESTERN_EDGE) then         ! Extrapolate elementary
!          do j=jstr,jend               ! differences near physical
!            FC(imin-1,j)=FC(imin,j)    ! boundaries to compensate.
!            rx(imin-1,j)=rx(imin,j)    ! for reduced loop ranges.
!          enddo
!        endif
!        if (EASTERN_EDGE) then
!          do j=jstr,jend
!            FC(imax+1,j)=FC(imax,j)
!            rx(imax+1,j)=rx(imax,j)
!          enddo
!        endif
!# endif

!++++
 
        do j=jstr,jend
          do i=istrU-1,iend
            cff=2.*FC(i,j)*FC(i+1,j)
            if (cff.gt.epsil) then
              dZx(i,j)=cff/(FC(i,j)+FC(i+1,j))
            else
              dZx(i,j)=0.
            endif
 
            cfr=2.*rx(i,j)*rx(i+1,j)
            if (cfr.gt.epsil) then
              dRx(i,j)=cfr/(rx(i,j)+rx(i+1,j))
            else
              dRx(i,j)=0.
            endif
          enddo               !--> discard FC, rx

          do i=istrU,iend
            ru_Sche(i,j,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*(         &
     &                              (P_Sche(i-1,j,k)-P_Sche(i,j,k))-HalfGRho*(    &

     &            (rhol(i,j,k)+rhol(i-1,j,k))*(z_r(i,j,k)-z_r(i-1,j,k)) &
 
     &   -OneFifth*( (dRx(i,j)-dRx(i-1,j))*( z_r(i,j,k)-z_r(i-1,j,k)  &
     &                            -OneTwelfth*(dZx(i,j)+dZx(i-1,j)) ) &
 
     &              -(dZx(i,j)-dZx(i-1,j))*( rhol(i,j,k)-rhol(i-1,j,k)  &
     &                            -OneTwelfth*(dRx(i,j)+dRx(i-1,j)) ) &
     &                                                            )))
          enddo
        enddo
!
! ETA-component of pressure gradient term:
!-------------- -- -------- -------- -----
!
        do j=jmin,jmax
          do i=istr,iend
            FC(i,j)=(z_r(i,j,k)-z_r(i,j-1,k))
            rx(i,j)=(rhol(i,j,k)-rhol(i,j-1,k))
          enddo
        enddo

! This part has been altered to take into account the boundary conditions:
!++++ Modified by JBG, 20120425
!# ifndef NS_PERIODIC
!        if (SOUTHERN_EDGE) then
          do i=istr,iend
            FC(i,jmin)=FC(i,jmin+1);FC(i,jmin-1)=FC(i,jmin)
            rx(i,jmin)=rx(i,jmin+1);rx(i,jmin-1)=rx(i,jmin)
          enddo
!        endif
!        if (NORTHERN_EDGE) then
          do i=istr,iend
            FC(i,jmax)=FC(i,jmax-1);FC(i,jmax+1)=FC(i,jmax)
            rx(i,jmax)=rx(i,jmax-1);rx(i,jmax+1)=rx(i,jmax)
          enddo
!        endif
!# endif

!++++
 
        do j=jstrV-1,jend
          do i=istr,iend
            cff=2.*FC(i,j)*FC(i,j+1)
            if (cff.gt.epsil) then
              dZx(i,j)=cff/(FC(i,j)+FC(i,j+1))
            else
              dZx(i,j)=0.
            endif
 
            cfr=2.*rx(i,j)*rx(i,j+1)
            if (cfr.gt.epsil) then
              dRx(i,j)=cfr/(rx(i,j)+rx(i,j+1))
            else
              dRx(i,j)=0.
            endif
          enddo               !--> discard FC, rx
 
          if (j.ge.jstrV) then
            do i=istr,iend     ! Eq. (R.1)
              rv_Sche(i,j,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*(        &
     &                             (P_Sche(i,j-1,k)-P_Sche(i,j,k)) -HalfGRho*(     &

     &            (rhol(i,j,k)+rhol(i,j-1,k))*(z_r(i,j,k)-z_r(i,j-1,k))  &
 
     &   -OneFifth*( (dRx(i,j)-dRx(i,j-1))*( z_r(i,j,k)-z_r(i,j-1,k)   &
     &                            -OneTwelfth*(dZx(i,j)+dZx(i,j-1)) )  &
 
     &              -(dZx(i,j)-dZx(i,j-1))*( rhol(i,j,k)-rhol(i,j-1,k)   &
     &                            -OneTwelfth*(dRx(i,j)+dRx(i,j-1)) )  &
     &                                                            )))
            enddo
          endif
        enddo
      enddo   !<-- k


!====================================================================================h

! -----------------------------------------------------------------
! Part III: The variables of ROMS are now translated into PSOM variables.
! -----------------------------------------------------------------

!=======================================
! Unneeded loop:
!do i=1,NI  
!  do j=1,NJ
!    do k=1,NK
!      ru2_Sche(i,j,k)=-(ru_Sche(i,j,k)*dn_u(i,j)*dm_v(i,j)/Hz(i,j,k))
!      rv2_Sche(i,j,k)=-(rv_Sche(i,j,k)*dn_u(i,j)*dm_v(i,j)/Hz(i,j,k))
!    enddo
!  enddo
!enddo

!--------------------------------
! ru_Sche and rv_Sche are the equivalent to grpifc and grpjfc but 
! - are dimensionalized
! - have a shift of one grid point.
! - the boundary conditions will have to be prescribed.


!-----------------------
!- nondimensionalization

adim_Sche=1/(P1/R0*dx*dy);   ! P1=R0.UL.UL/EPS, where UL=FPAR.LEN.EPS.

! The minus sign comes from the fact that the ROMS code output is: ru = - Hz * dP/dx
ru_Sche=(-1.)*ru_Sche*adim_Sche;
rv_Sche=(-1.)*rv_Sche*adim_Sche;


!-----------------------
!- shift

ru2_Sche=-999.; rv2_Sche=-999.;

do i=1,NI-1
  do j=1,NJ
    do k=1,NK
      ru2_Sche(i,j,k)=ru_Sche(i+1,j,k)
    enddo
  enddo
enddo


do i=1,NI
  do j=1,NJ-1
    do k=1,NK
      rv2_Sche(i,j,k)=rv_Sche(i,j+1,k)
    enddo
  enddo
enddo

!-----------------------
!- boundary conditions

!- in x:
ru2_Sche(0,:,:)=ru_Sche(1,:,:);
ru2_Sche(NI,:,:)=ru_Sche(1,:,:);

!- in y:
rv2_Sche(:,0,:)=0.;
rv2_Sche(:,NJ,:)=0.;


!--------------------------------
! ru2_Sche and rv2_Sche are now totally equivalent to grpifc and grpjfc.
! Now, we have to compute drpx and drpy.
! The relation between grpifc/grpjfc and drpx/drpy is exactly the same as the one in the Song scheme.

do i=1,NI
  do j=1,NJ
    do k=1,NK
      ru3_Sche(i,j,k)=0.5d0*(ru2_Sche(i,j,k)/gi(i,j,k,1)+ru2_Sche(i-1,j,k)/gi(i-1,j,k,1))
    enddo
  enddo
enddo

do i=1,NI
  do j=1,NJ
    do k=1,NK
      rv3_Sche(i,j,k)=0.5d0*(rv2_Sche(i,j,k)/gj(i,j,k,2)+rv2_Sche(i,j-1,k)/gj(i,j-1,k,2))
    enddo
  enddo
enddo

do i=1,NI
  do j=1,NJ
    do k=1,NK
      ru4_Sche(i,j,k)=(ru3_Sche(i,j,k)*ux(i,j)+rv3_Sche(i,j,k)*vx(i,j))
    enddo
  enddo
enddo

do i=1,NI
  do j=1,NJ
    do k=1,NK
      rv4_Sche(i,j,k)=(ru3_Sche(i,j,k)*uy(i,j)+rv3_Sche(i,j,k)*vy(i,j))
    enddo
  enddo
enddo


!---------------------------------
! ru2_Sche and rv2_Sche are the grpifc and grpjfc
! ru4_Sche and rv4_Sche are the drpx and drpy
!---------------------------------


!---------------------------------
! _ End of the code.
!---------------------------------


!---------------------------------
! To illustrate how the code works on a simple case, 
! Let us set:
!  * rho=R0+ constant
!  * h varies in y.
!
! The pressure terms can be decomposed into:
!
! (              Eq. P.2) only pressure terms from NK-1 to 1:                         0           0.737      0.737      0.737    ...
! (part of Eq. P.1 + P.2) pressure terms from NK-1 to 1 + the rho term in NK face:    0.737       1.474      1.474      1.474    ...
! (part of Eq. P.1 + P.2) pressure terms from NK-1 to 1 + the g.z_w term:           -55.75      -55.02     -55.02     -55.02     ...
! (        Eq. P.1 + P.2) all pressure terms:				            -55.01      -54.28     -54.28     -54.28     ...
!                                                                                                        
! Then, a local additional term is added:                                                               
!                                                                                                        
! (part of Eq. R.1)       only this term:				 	      0.737       0.         0.         0.       ...
!                                                                                                        
! (Eq. P.1, P.2 and R.1) In the end, the addition of all these terms gives:         -54.276     -54.285    -54.285   -54.285     ...
!---------------------------------

      return
      end
 
