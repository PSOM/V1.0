subroutine conadjust(stepl,n) 
  !----------------------------------------------------                   
  USE header
  !     T(i,j,k) is set to 0, so Treal is T0                              
  !     Performs convective adjustment                                    
  integer i,j,k,n,stepl 
  REAL(kind=rc_kind) :: rh,rhtop,zt,zb,dz,zinv 
  REAL(kind=rc_kind) :: dzupper 
  REAL(kind=rc_kind) :: sreal,Treal,rhreal,pbar 
  REAL(kind=rc_kind) ::  potdens
  ! assume that evalrho has been called
                                                                          
  do j=1,NJ 
    do i=1,NI 
      do k=NK-1,1,-1 
        conv(i,j,k)= 0 
     
        ! top level                                                                   
        rhtop= rho(i,j,k+1); 
        dzupper = zf(i,j,k+1) -zf(i,j,k) 

        ! this level                                                        
        rh =  rho(i,j,k) 
        dz = zf(i,j,k) -zf(i,j,k-1) 
                                                                        
        if (rh.lt.rhtop) then !     mix                                                               
          conv(i,j,k)= 1 
          zinv= 1.d0/(dzupper +dz) 
          s(i,j,k,n)=(dzupper*s(i,j,k+1,n) +dz*s(i,j,k,n))*zinv 
          s(i,j,k+1,n)= s(i,j,k,n) 
          T(i,j,k,n)=(dzupper*T(i,j,k+1,n) +dz*T(i,j,k,n))*zinv 
          T(i,j,k+1,n)= T(i,j,k,n) 
          rho(i,j,k)= potdens(s(i,j,k,n),T(i,j,k,n)) 
          rho(i,j,k+1) = rho(i,j,k)
          do it=1,ntr 
             Tr(it,i,j,k,n)=(dzupper*Tr(it,i,j,k+1,n)+dz*Tr(it,i,j,k,n))*zinv                       
             Tr(it,i,j,k+1,n)= Tr(it,i,j,k,n) 
          end do
        end if
                                                                       
       end do
    end do
  end do
                                                                        
  if (mod((stepl-1),100).eq. 0) then 
     con100= 0
  else 
     con100= con100 +conv
  end if

  return 
END subroutine conadjust
