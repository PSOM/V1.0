subroutine ini_st 
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
  real :: coef,coef2,coef3,zeff,rv,zv
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
       s(i,j,k,1)=s(i,j,k,1)-0.4*z*user1/100.
     enddo
   enddo
 enddo

 s(:,:,:,0)=s(:,:,:,0)+s(:,:,:,1); s(:,:,:,1)=0.;


 ! BAROCLINIC AREA ON THE SHELF
 do j=1,NJ+1 
   do i=0,NI+1 
     do k=0,NK+1 
       z= DL*zc(i,j,k) 
       jlim=124.
       coef2=0.5*( tanh(3.*(110.+z)/110.*3.1415926) +1.)
       coef3=((tanh((112.-j)/4.)+1.)*0.5);
       coef2=coef3*1+(1.-coef3)*coef2;
        !s(i,j,k,1)=coef2
        !s(i,j,k,1)=s(i,j-1,k,1)+user7*(tanh((REAL((j-0.5-jlim))*4.5)/NJ*PI)-1.)/2.*2./((D(i,j)+D(i,j-1))*DL)
        s(i,j,k,1)=-user7*1.2*(jlim-j)/jlim*((tanh((jlim-j)/8.)+1.)*0.5)*coef2
     enddo
   enddo
 enddo

 s(:,:,:,0)=s(:,:,:,0)+s(:,:,:,1); s(:,:,:,1)=0.;


 ! JET
 do j=1,NJ+1 
   do i=0,NI+1 
     do k=0,NK 
       z= DL*zc(i,j,k) 
          jlim=120.
          rv=20.
          !coef2=exp(zc(i,j,k)/(0.1*DL))*NK/(zf(i,j,k+1)-zf(i,j,k))

          coef2=0.5*( tanh(3.*(100.+z)/100.*3.1415926) +1.)
          s(i,j,k,1)=0.15*coef2* 0.5*(tanh((j-jlim)/25.*3.1415926)-1.)

          !s(i,j,k,1)=s(i,j-1,k,1)-user2*(exp(-(REAL(j-0.5-jlim)**2/rv**2)))*3./((D(i,j)+D(i,j-1))*DL)*coef2
          ! s(i,j,k,1)=s(i,j-1,k,1)-user2*(exp(-(REAL(j-0.5-jlim)**2/rv**2)))*3.*coef2
          !s(i,j,k,1)=coef2;! s(i,j-1,k,1)-user2*(exp(-(REAL(j-0.5-jlim)**2/rv**2)))*3.*coef2
     enddo
   enddo
 enddo

 s(:,:,:,0)=s(:,:,:,0)+s(:,:,:,1); s(:,:,:,1)=0.;


 PRINT*," W ",user5
 IF (ABS(user5-1000.)<100.) then
   PRINT*," W ",user5
   do j=0,NJ+1 
     do i=0,NI+1 
       do k=0,NK+1 

         xcenter=48.; ycenter=103.+(user5-1000.); dyvor=10.*dy/1e3*user6; dxvor=10.*dx/1e3*user6;
         s(i,j,k,1)=s(i,j,k,1)+REAL(user6)*user8*0.05*dexp(-(((yc(j)-ycenter)/dyvor)**2+((xc(i)-xcenter)/dxvor)**2));

       enddo
     enddo
   enddo

 ENDIF



 IF((user5-2000.)<200.) then

  ! WIGGLE
  do j=0,NJ+1 
    do i=0,NI+1 
      do k=0,NK+1 
        z= DL*zc(i,j,k) 

!         IF(.TRUE.) then
        IF(REAL(i)>REAL(NI)/4. .AND. REAL(i)<(3*REAL(NI)/4.)) then
          jlim=93.+(1.+z/100.)*2.-user6*2.*2.5*cos(2*3.14159265358979323*REAL(i)/96.*2.)+(user5-2000.)
         ELSE
          jlim=93.+(1.+z/100.)*2.-user6*2.*2.5*(-1.)                                    +(user5-2000.)
        ENDIF               

        coef=0-(user8)*(tanh((REAL((j-jlim))*10.)/NJ*PI)-1.)/2.
                    ! *(tanh((REAL((z+100.))*3)/NK*PI)-1.)/2. !* (zc(NI/2,10,1)/zc(i,j,1))**1.
        s(i,j,k,1)=s(i,j,k,1)-0.1*coef
      enddo
    enddo
  enddo

  ELSE





 ENDIF



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
