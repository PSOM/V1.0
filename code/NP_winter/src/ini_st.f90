subroutine ini_st
    ! --------------------
    ! initializes s and t with a single vertical profile of Salinity and Temperature
    ! Across-front use the TANH profile with tightness as a measure of frontal spread
    ! Larger the factor, tighter the front and larger the gradient.

    USE header
    implicit none
    ! Number of datapoints in the vertical profile
    integer npl
    parameter (npl=1001)
    double precision Svert(npl),Tvert(npl),dep(npl),bS1(npl),cS1(npl),dS1(npl), &
        & bT1(npl),cT1(npl),dT1(npl),z,seval,z1,z2,yfront1,yfront2, &
        & yfront3
    double precision sbkgrnd, tbkgrnd
    integer  i,j,k
    double precision slfac,Tlfac,thy,slfacnew,Tlfacnew,wiggles,amplitude
    integer njf1, njf2, njf3, njf4

    ! tightness = 10 represents a very tight front, =1 loose front
    ! mldepth is the desired ML depth

    ! Read in Temperature profile
    open(unit=31,file='./input/init_T_profile.in')
    do k=1,npl
        read(31,*) dep(k),Tvert(k)
    end do
    close(31)
    ! Read in salinity profile
    open(unit=31,file='./input/init_S_profile.in')
    do k=1,npl
        read(31,*) dep(k),Svert(k)
    end do
    close(31)
        
    ! Controls depth of frontal zone
    mldepth= 90.d0
    z1= mldepth
    z2= 300.d0  ! increased to get larger eddies for C export

    ! Tightness of the front. tightness should be tuned to get the appropriate 
    ! density gradient throughout the domain.
    tightness= 0.3d0

    call spline (npl,dep,Tvert,bT1,cT1,dT1)
    call spline (npl,dep,Svert,bS1,cS1,dS1)

    ! sclwidth is the scale width of the channel (i.e. orig width = 48km)
    ! yfront is the distance of the front from the coast
    sclwidth = 48.0

    yfront2 = 0.5*(yc(NJ-1) +yc(0))
    yfront1 = 0.5*yfront2
    yfront3 = 1.5*yfront2

!     z1 is the depth of the ml (diff rho on both sides above this depth)
!     z2 is the vertical extent of the density anamoly (it is gradually 
!        linearly anihillated with depth).

    ! half of the range covered by 1 front (I.e. 1/6 of the background S,T,
    ! or density range over the domain). Temperature coef negative as T decreases with latitude
    slfac = 0.009418
    Tlfac = -0.234327

    ! The density fronts are created starting from northern boundary, 
    ! but the profiles provided are from Station Papa, 
    ! i.e., the middle of the domain. The profiles are then 
    ! corrected by a factor to obtained the profiles at the northern boundary.
    Tvert = Tvert+3*Tlfac
    Svert = Svert+3*slfac
    
      do j=0,NJ+1
          do i=0,NI+1
              do k=0,NK+1
                  z= DL*zc(i,j,k)
                  z=-z
                  sbkgrnd =  seval(npl,z,dep,Svert,bS1,cS1,dS1) 
                  s(i,j,k,0)=  sbkgrnd
                  tbkgrnd =  seval(npl,z,dep,Tvert,bT1,cT1,dT1) 
                  T(i,j,k,0)=  tbkgrnd
              end do
          end do
      end do

      njf1 = 0.375*NJ
      njf2 = 0.625*NJ
      njf3 = njf1
      njf4 = njf2

      wiggles=2.d0
      amplitude= 7.d0

      !----------
      ! FRONT 3
      ! ---------
      do j = njf2+1,NJ+1
          do i=0,NI+1
              ! WIGGLE in FRONT
              yfront= yfront3 + amplitude*dsin(2.d0*PI*xc(i)/(0.5*(xc(NI+1)+xc(NI))/wiggles))
              thy = tanh(tightness*(yc(j)-yfront)*PI)
              do k=1,NK
                  z= DL*zc(i,j,k)
                  if (z.gt.-z2) then
                      if (z.ge.-z1) then 
                          slfacnew = slfac
                          Tlfacnew = Tlfac
                      else if (z.ge.-z2) then
                          slfacnew = slfac*(z+z2)/(z2-z1)
                          Tlfacnew = Tlfac*(z+z2)/(z2-z1)
                      else
                          slfacnew = 0.d0
                          Tlfacnew = 0.d0
                      end if

                      s(i,j,k,0)= slfacnew*(thy-1.d0) + s(i,NJ,k,0)
                      T(i,j,k,0)= Tlfacnew*(thy-1.d0) + T(i,NJ,k,0)
                  end if
              end do
          end do
      end do

      !----------
      ! FRONT 2
      ! ---------
      do j = njf1+1,njf2
          do i=0,NI+1
              ! WIGGLE in FRONT
              yfront = yfront2 + amplitude*dsin(2.d0*PI*xc(i)/(0.5*(xc(NI+1)+xc(NI))/wiggles))
              thy = tanh(tightness*(yc(j)-yfront)*PI)
              do k=1,NK
                  z= DL*zc(i,j,k)
                  if (z.gt.-z2) then
                      if (z.ge.-z1) then 
                          slfacnew = slfac
                          Tlfacnew = Tlfac
                      else if (z.ge.-z2) then
                          slfacnew = slfac*(z+z2)/(z2-z1)
                          Tlfacnew = Tlfac*(z+z2)/(z2-z1)
                      else
                          slfacnew = 0.d0
                          Tlfacnew = 0.d0
                      end if

                      s(i,j,k,0)= slfacnew*(thy-1.d0) + s(i,njf2+1,k,0)
                      T(i,j,k,0)= Tlfacnew*(thy-1.d0) + T(i,njf2+1,k,0)
                  end if
              end do
          end do
      end do

      ! ----------
      ! FRONT 1
      ! ---------
      do j=0,njf1
          do i=0,NI+1
              ! WIGGLE in FRONT
              yfront= yfront1 + amplitude*dsin(2.d0*PI*xc(i)/(0.5*(xc(NI+1)+xc(NI))/wiggles))
              thy = tanh(tightness*(yc(j)-yfront)*PI)
              do k=1,NK
                  z= DL*zc(i,j,k)
                  if (z.gt.-z2) then
                      if (z.ge.-z1) then 
                          slfacnew = slfac
                          Tlfacnew = Tlfac
                      else if (z.ge.-z2) then
                          slfacnew = slfac*(z+z2)/(z2-z1)
                          Tlfacnew = Tlfac*(z+z2)/(z2-z1)
                      else
                          slfacnew = 0.d0
                          Tlfacnew = 0.d0
                      end if

                      s(i,j,k,0)= slfacnew*(thy-1.d0) + s(i,njf1+1,k,0)
                      T(i,j,k,0)= Tlfacnew*(thy-1.d0) + T(i,njf1+1,k,0)
                  end if
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
          ! periodicew
          do j=0,NJ+1
              s(0,j,k,0)= s(NI,j,k,0)
              s(NI+1,j,k,0)= s(1,j,k,0)
              T(0,j,k,0)= T(NI,j,k,0)
              T(NI+1,j,k,0)= T(1,j,k,0)              
          end do
      end do

      return
end

