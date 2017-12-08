subroutine stprofile 
  !     --------------------                                              
  USE header
!     initializes s as pot.density with a single vertical profile       
!     Across-front use the TANH profile with tightness as a measure of f
!     spread                                                            
!     Larger the factor, tighter the front and larger b_xx.             
!     tightness=4, gives the needed b_xx to make b_xx H/ (2f^2) < 1.    
  IMPLICIT NONE 
  integer  i,j,k,iseed,nplm1,iter,iter2 
  integer ic,jc
  real :: jlim,jloc
  real :: coef,coef2,zeff,rv,zv
  real :: intdepth
  integer, parameter :: npl=33
  REAL(kind=rc_kind) :: rhovert(npl),Tvert(npl),svert(npl),dep(npl),depoff(npl) 
  REAL(kind=rc_kind) :: bs1(npl),cs1(npl),ds1(npl),bT1(npl),cT1(npl),dT1(npl),z,seval,zmm, &
 & sbkgrnd,Tbkgrnd,z1,z2,zprev,alpha,dtdz,dsdz
  REAL(kind=rc_kind) :: slfac,dscl,rdscl,yarg,ex2y,thy,ran3,perturb,slfacnew,dz,bfsqbkgrnd,wiggles,amplitude             
                                                                   
  REAL :: dxvor, dyvor, xcenter, ycenter
 
  ! s(:,:,:,0)=1025.5;  s(:,:,:,1)=0.;


  ! BACKGROUND STRATIFICATION
  do j=0,NJ+1 
    do i=0,NI+1 
      do k=0,NK+1 
        z= DL*zc(i,j,k) 
        s(i,j,k,1)=s(i,j,k,1)-1.*z*0.5/100.
      enddo
    enddo
  enddo

 s(:,:,:,0)=s(:,:,:,0)+s(:,:,:,1); s(:,:,:,1)=0.;


  ! BAROCLINIC AREA ON THE SHELF
  do j=1,NJ+1 
    do i=0,NI+1 
      do k=0,NK+1 
        z= DL*zc(i,j,k) 
        jlim=118.
         s(i,j,k,1)=s(i,j-1,k,1)+1.*(tanh((REAL((j-0.5-jlim))*1.5)/NJ*PI)-1.)/2.*2./((D(i,j)+D(i,j-1))*DL)
      enddo
    enddo
  enddo

 s(:,:,:,0)=s(:,:,:,0)+s(:,:,:,1); s(:,:,:,1)=0.;

  ! JET
  do j=1,NJ+1 
    do i=0,NI+1 
      do k=0,NK 
        z= DL*zc(i,j,k) 
           jlim=118.
           rv=20.
           !coef2=exp(zc(i,j,k)/(0.1*DL))*NK/(zf(i,j,k+1)-zf(i,j,k))
           coef2=1.
           s(i,j,k,1)=s(i,j-1,k,1)-1.*(exp(-(REAL(j-0.5-jlim)**2/rv**2)))*3./((D(i,j)+D(i,j-1))*DL)*coef2
      enddo
    enddo
  enddo

 s(:,:,:,0)=s(:,:,:,0)+s(:,:,:,1); s(:,:,:,1)=0.;

!@
! s(:,:,:,0)=0.; s(:,:,:,1)=0.;


!  do j=0,NJ+1 
!    do i=0,NI+1 
!      do k=0,NK+1 

!        xcenter=20.; ycenter=90.; dyvor=10.*dy/1e3; dxvor=10.*dx/1e3;
!        s(i,j,k,1)=s(i,j,k,1)+1.*0.1*dexp(-(((yc(j)-ycenter)/dyvor)**2+((xc(i)-xcenter)/dxvor)**2));

!        xcenter=68.; ycenter=90.; dyvor=10.*dy/1e3; dxvor=10.*dx/1e3;
!        s(i,j,k,1)=s(i,j,k,1)-1.*0.1*dexp(-(((yc(j)-ycenter)/dyvor)**2+((xc(i)-xcenter)/dxvor)**2));

!      enddo
!    enddo
!  enddo



  ! WIGGLE ON THE SHELF
  do j=0,NJ+1 
    do i=0,NI+1 
      do k=0,NK+1 
        z= DL*zc(i,j,k) 
            jlim=92.+(1.+z/100.)*2.+2.*2.5*cos(2*3.14159265358979323*REAL(i)/NI*2.)
            coef=0-1.*(tanh((REAL((j-jlim))*10.)/NJ*PI)-1.)/2.
                          ! *(tanh((REAL((z+100.))*3)/NK*PI)-1.)/2. !* (zc(NI/2,10,1)/zc(i,j,1))**1.
            s(i,j,k,1)=s(i,j,k,1)-0.1*coef
      enddo
    enddo
  enddo



!@
! s(i,j,k,1)=0.
! do j=0,NJ+1; do i=0,NI+1; do k=0,NK+1
!   s(i,j,k,1)=s(i,j,k,1)-0.05*((tanh(REAL((i-45.)*10.)/NI*PI)-1.)/2.-(tanh(REAL((i-85.)*10.)/NI*PI)-1.)/2.)*((tanh(REAL((j-75.)*10.)/NI*PI)-1.)/2.-(tanh(REAL((j-150.)*10.)/NI*PI)-1.)/2.)
! enddo; enddo; enddo;

 s(:,:,:,0)=s(:,:,:,0)+s(:,:,:,1); s(:,:,:,1)=0.;



 s(:,:,:,0)=s(:,:,:,0)+1025.5;  s(:,:,:,1)=0.;

  do k=0,NK+1 
     do i=1,NI 
        s(i,0,k,0)= s(i,1,k,0) 
        s(i,NJ+1,k,0)= s(i,NJ,k,0) 
        T(i,0,k,0)= T(i,1,k,0) 
        T(i,NJ+1,k,0)= T(i,NJ,k,0) 
     end do 
 ! periodicew                                                        
     do j=0,NJ+1 
        s(0,j,k,0)= s(NI,j,k,0) 
        s(NI+1,j,k,0)= s(1,j,k,0) 
        T(0,j,k,0)= T(NI,j,k,0) 
        T(NI+1,j,k,0)= T(1,j,k,0) 
     end do 
  end do 
                                                                    
                                                                    
  return 
  END                                           
