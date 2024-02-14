subroutine ini_st
!    --------------------
!    DESCRIPTION
!    Initializes density for model using an analytic
!    density function, derived from an N^2 profile
!    Across-front use the TANH profile with tightness
!    as a measure of frontal spread
!    Larger the factor, tighter the front and larger b_xx.

   USE header 
   implicit none 
   integer i,j,k,n
   REAL(kind=rc_kind) :: grav,f,N2bkgrnd,dTdz,dz
   REAL(kind=rc_kind) :: thy,z_topog,z_front,z_fac,y_center,y_front,dUdz_max,dT
   parameter (grav=9.81, f=1.3d-4) !@ alpha: thermal expansion coefficient

   !---------------------------------------
   !           Background s
   !---------------------------------------
   s = 35.0     !set salnity constant, front is only in temperature
   !---------------------------------------
   !      Set up constant stratification
   !---------------------------------------
   N2bkgrnd = Nbkgrnd**2        ! in rad/s
   dTdz = -N2bkgrnd*R0/grav ! d(rho)/dz = -N2bkgrnd *R0/grav      
   n=0
   T(:,:,NK+1,n) = R0       ! surface density
   do j=0,NJ+1
      do i=0,NI+1
         do k=NK,0,-1
            dz = (zc(i,j,k+1) - zc(i,j,k))*DL
            T(i,j,k,n) = T(i,j,k+1,n) - dTdz*dz
         end do
      end do
   end do

!    write(6,"(A,4(D15.7))") " # T (top/bot): " ,T(i,j,k,n)-R0,T(i,j,1,n)-R0

! =================================================================================================
!                                          Add front
! =================================================================================================
   z_topog = -depmean_dim/DL        ! non-dimensional z-coordinate of the averaged topography
   z_front = 0.05*z_topog           ! non-dimensional z-coordinate of front's upper bound
   dUdz_max= 2*Umax/(z_topog-z_front)! max. value of dU/dz at the central bottom, m/s per km
   
!    write(6,"(A,1(D15.7))") " # dUdz: " ,dUdz_max/1000

   y_center= 0.5 *(yc(NJ+1) +yc(0)) ! non-dimensional y-coordinate of front's centerline
   y_front = 0.1 *(yc(NJ+1) +yc(0)) ! non-dimensional y-coordinate of front's lateral bound (east)
   tightness= PI/(y_front-y_center) ! per km
   dT=f*dUdz_max*R0/grav/tightness  ! density difference across the front

!    write(6,"(A,4(D15.7))") " # y1, yN, Ly, dy" ,yc(1),yc(NJ),yc(NJ)-yc(0), yc(2)-yc(1)
!    write(6,"(A,2(D15.7))") " # tightness, dT: ", tightness,dT


   do j=0,NJ+1
   do i=0,NI+1
      thy = tanh((yc(j)-y_center)*tightness) ! between +-1
      do k=0,NK+1
         if (zc(i,j,k).lt.z_front) then
            z_fac      = (zc(i,j,k)-z_front)/(z_topog-z_front) ! between 0 and 1
            T(i,j,k,0) = z_fac*dT*thy + T(i,NJ,k,n)
         else
            T(i,j,k,0) = T(i,NJ,k,n)
         end if
      end do ! end k
   end do ! end i
   end do ! end j

! =================================================================================================
!                                        Save initial s, T
! ================================================================================================= 
   do j=0,NJ+1
   do k=0,NK+1
      T_ini(j,k)= 0.d0
      do i=0,NI+1
         T_ini(j,k)= T_ini(j,k) + T(i,j,k,0)
      enddo
      T_ini(j,k)= T_ini(j,k)/(dble(NI)+2)
      s_ini(j,k)= 35.0
   enddo
   enddo
   
!    write(6,"(A,4(D15.7))") " # T (east/west): " ,T(NI/2,1,1,0)-R0,T(NI/2,NJ,1,0)-R0

   
   return
end