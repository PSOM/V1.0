subroutine tracersource(n,dtimel) 
  !     -----------------------                                           
  USE header

  !     tracer sources                                                    
  !                                                                       
  integer  i,j,k,n,m 
  REAL(kind=rc_kind) :: dtimel,fac,tau,tauinv,vol
!  REAL(kind=rc_kind) :: trinit,facmortal,phyconvert 
  !                                                                       
  !     consump is the uptake (or if negative, then addition) of tracer   

  !tau is the time scale for consump = 1day = 86400 seconds
  tauinv = 1.d0/(1.d0*86400.d0)    ! tauinv is in (per second)
  fac= dtimel/(1.d0*86400.d0*UL/LEN) ! dtime is non-dim time, TL=LEN/UL
  ! fac is non-dim

  it=1
  do k=0,NK+1 
     do j=0,NJ+1 
        do i=0,NI+1 
           vol= Jac(i,j,k)*LEN*LEN*DL   !vol in m^3
           consump(i,j,k,it)= vol*tauinv*                           &
                (Tr(it,i,j,k,n)- trinit(j,k))  !Tr in milimoles per m^3
           ! consump in milimoles/ sec
           Tr(it,i,j,k,n)= Tr(it,i,j,k,n) - fac*(Tr(it,i,j,k,n)- trinit(j,k))
        end do
     end do
  end do
  !
  ! prod is consump averaged in the x direction
  ! Multiply by 1.d-3 to get this in MOLES per SECOND
  do k=1,NK
     do j=1,NJ
        prod(j,k)= 0.d0
        do i=1,NI
           prod(j,k)= prod(j,k)+ consump(i,j,k,it)
        end do
        prod(j,k)= prod(j,k)*1.d-3/dble(NI)
     end do
  end do
  !                                                                       
  do m=0,1 
     do k=0,NK+1 
        do j=0,NJ+1 
           do it=1,ntr 
              tr(it,0,j,k,m)= tr(it,NI,j,k,m) 
              tr(it,NI+1,j,k,m)= tr(it,1,j,k,m) 
           end do
        end do
     end do
  end do
  !                                                                       
  return 
end subroutine tracersource
