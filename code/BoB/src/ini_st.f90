!subroutine stprofile 
subroutine ini_st
  !     --------------------                                              
  USE header
!     initializes s as pot.density with a single vertical profile       
!     Across-front use the TANH profile with tightness as a measure of f
!     spread                                                            
!     Larger the factor, tighter the front and larger b_xx.             
!     tightness=4, gives the needed b_xx to make b_xx H/ (2f^2) < 1.    
  implicit none
  INTEGER, PARAMETER                 :: nps=148, nwork =50000
  real(kind=rc_kind) depth(nps),S_south(nps),S_north(nps),T_south(nps),T_north(nps)
!  REAL*8 :: bs(np),cs(np),ds(np), bt(np),ct(np),dt(np),                 &
!       xs(np),ys(np),xt(np),yt(np)
  !integer  i,j,k,it,iseed,npm1,n
  integer  i,j,k,iseed,npm1,n

  real(kind=rc_kind) :: dep(nps)
  real(kind=rc_kind) ::  bs1(nps),cs1(nps),ds1(nps),bT1(nps),cT1(nps),dT1(nps),  &
       sal(nps),temp(nps),z,seval,         &
       zmm,sbkgrnd,z1,z2,potdens
  !real(kind=rc_kind) :: slfac,dum,dscl,rdscl,yarg,ex2y,thy,ran3,              &
  real(kind=rc_kind) :: slfac,dscl,rdscl,yarg,ex2y,thy,ran3,              &
       perturb,slfacnew,zoffset,dz,bfsqbkgrnd,efold,dsalt(0:NJ+1),       &
       sd,am,da(NI)
!     tightness = 10 represents a very tight front, =1 loose front      
!=      parameter (zoffset= 200.)  also used zoffset=50 ! if zoffset=0 m
!=       parameter (zoffset= 200.d0, tightness=0.03)                    
!=       parameter (zoffset= -10.d0)                                    
!       parameter (zoffset= -20.d0) 
  real(kind=rc_kind) :: dwork(nwork)
  parameter (zoffset= 0.d0)
  tightness=0.01d0
  open(unit=40,file='./BoB/Bengal_sT_south_north.in')
  do i=1,nps
     read(40,*) dep(i),s_south(i),T_south(i),s_north(i),T_north(i)
     write(41,11) dble(i),dep(i),s_south(i),T_south(i),s_north(i),T_north(i)
11 format (6(F10.4))
  end do
  close(40)

  !     The z values must be in STRICTLY INCREASING ORDER and hence should
  !     be written from ocean bottom up.                                  
                                                                        
  n=0
  it=1
  call spline (nps,dep,s_south,bs1,cs1,ds1) 
  call spline (nps,dep,T_south,bt1,ct1,dt1) 

  i=NI/2
  j=NJ/2
  do k=1,NK
     z= DL*zc(i,j,k) 
     ssouth(k)= seval(nps,z,dep,S_south,bs1,cs1,ds1)   
     Tsouth(k)= seval(nps,z,dep,T_south,bt1,ct1,dt1) 
     write(60,*) k,ssouth(k),Tsouth(k)
  end do

  call spline (nps,dep,s_north,bs1,cs1,ds1) 
  call spline (nps,dep,T_north,bt1,ct1,dt1) 

  i=NI/2
  j=NJ/2
  do k=1,NK
     z= DL*zc(i,j,k) 
     snorth(k)= seval(nps,z,dep,S_north,bs1,cs1,ds1)   
     Tnorth(k)= seval(nps,z,dep,T_north,bt1,ct1,dt1) 
     write(61,*) k,snorth(k),Tnorth(k)
  end do

!  do i=1,NI
!     do k=1,NK
!        do j=0,NJ/3
!           s(i,j,k,n)= ssouth(k)
!          T(i,j,k,n)= Tsouth(k)
!        end do
!        do j=2*NJ/3, NJ
!           s(i,j,k,n)= snorth(k)
!           T(i,j,k,n)= Tnorth(k)
!        end do
!        do j=NJ/3+1,2*NJ/3-1
!           s(i,j,k,n)= snorth(k)
!           T(i,j,k,n)= Tnorth(k)
!        end do
!     end do
!  end do
 
  ! Perturb the FRONT                                                   
  ! Add random perturbation to front (white noise)
  iseed= 44294
  am=  0.d0    ! mean
  sd= 0.5d0   !standard deviaton
  call rannw (am, sd, iseed, da, NI, dwork, nwork)
  !     This generates a ran number array NI long, with mean 0 and sd=.1        

  do j=0,NJ+1 
     !     do i=0,NI+1 
     yfront= 0.5d0*(yc(NJ+1) +yc(0))                           !  &
!     &           + da(i)
     thy = 0.5 + 0.5*tanh(tightness*(yc(j)-yfront)*2.*PI) 
     !now thy goes from 0 to 1 ,  thy=0 at j=0
!     write(6,*) j,thy
     do i=0,NI+1
        do k=1,NK
            s(i,j,k,0)= ssouth(k) *(1.d0-thy) + snorth(k)*thy
            T(i,j,k,0)= Tsouth(k) *(1.d0-thy) + Tnorth(k)*thy
         end do
	 s(i,j,0,0)=s(i,j,1,0)
	 T(i,j,0,0)=T(i,j,1,0)
	 s(i,j,NK+1,0)=s(i,j,NK,0)
	 T(i,j,NK+1,0)=T(i,j,NK,0)
      end do
  end do
  !call evalrho(rho,0) 

  return 
END subroutine ini_st
