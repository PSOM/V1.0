subroutine ini_st 
  !     --------------------                                              
  USE header
!     initializes s as pot.density with a single vertical profile       
  IMPLICIT NONE 
  integer  i,j,k,iseed,nplm1,iter,iter2 
  integer, parameter :: npl=33
  REAL(kind=rc_kind) :: rhovert(npl),Tvert(npl),svert(npl),dep(npl),depoff(npl) 
  REAL(kind=rc_kind) :: bs1(npl),cs1(npl),ds1(npl),bT1(npl),cT1(npl),dT1(npl), &
     &     z,seval,zmm,sbkgrnd,Tbkgrnd,z1,z2,zprev,alpha,dtdz,dsdz
  REAL(kind=rc_kind) :: slfac,dscl,rdscl,yarg,ex2y,thy,ran3,         &
     &     perturb,slfacnew,dz,bfsqbkgrnd,wiggles,amplitude             
  data dep / -5500.d0, -5000.d0, -4500.d0, -4000.d0, -3500.d0,      &
     &      -3000.d0, -2500.d0, -2000.d0, -1750.d0, -1500.d0, -1400.d0, &
     &      -1300.d0, -1200.d0, -1100.d0, -1000.d0, -900.d0, -800.d0,   &
     &      -700.d0, -600.d0, -500.d0, -400.d0, -300.d0, -250.d0,       &
     &      -200.d0, -150.d0, -125.d0, -100.d0, -75.d0, -50.d0,         &
     &      -30.d0, -20.d0, -10.d0, -0.d0 /                             
  data svert / 34.76,  34.87, 34.91, 34.94,  35.00,  35.05, 35.08,  &
          35.11,  35.13,  35.13,  35.11,  35.05,  34.99,  34.94,       &
           34.89,  34.91,  34.89,  34.89,  34.88,  34.89,  34.90,       &
           34.88,  34.88,  34.88,  34.88,  34.88,  34.87,  34.81,       &
           34.76,  34.75,  34.75,  34.74,  34.71 /                            
  data tvert / 1.61, 1.72, 1.77, 1.78, 1.87, 2.18, 2.61, 3.05,      &
           3.28, 3.62, 3.80, 4.06, 4.32, 4.57, 4.84, 5.29, 5.71,        &
           6.37, 6.96, 7.83, 8.75, 9.82, 10.42, 11.00, 11.71, 12.12,    &
           12.63, 13.28, 14.46, 16.01, 16.76, 17.04, 17.25/


  !     Specify MLD                                                       
  mldepth= 50.d0 
  dsdz = 0.d0

  do k=1,npl 
    depoff(k)= dep(k)- mldepth 
  end do 
                                                                                                                                                
  do j=0,NJ+1 
    do i=0,NI+1 
      do k=0,NK+1 
        z= DL*zc(i,j,k) 

         T(i,j,k,0) = 29.95d0 + z*0.0454d0        
!         T(i,j,k,0) = ( (5/200)*z    ) + 30.d0         
!          T(i,j,k,0) = T(i,j,k,0) + 0.2d0          
         s(i,j,k,0) = 35.d0            

      end do 
    end do 
  end do 
 
                                                            
  do k=0,NK+1 
    do i=1,NI 
      s(i,0,k,0)= s(i,1,k,0) 
      s(i,NJ+1,k,0)= s(i,NJ,k,0) 
      T(i,0,k,0)= T(i,1,k,0) 
      T(i,NJ+1,k,0)= T(i,NJ,k,0) 
    end do 
    do j=0,NJ+1 
      s(0,j,k,0)= s(NI,j,k,0) 
      s(NI+1,j,k,0)= s(1,j,k,0) 
      T(0,j,k,0)= T(NI,j,k,0) 
      T(NI+1,j,k,0)= T(1,j,k,0) 
    end do 
  end do 
                                                                        
                                                                        
return 
END                                           
