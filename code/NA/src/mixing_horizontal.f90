subroutine mixing_horizontal(udif,vdif,wdif,sdif,Tdif,Trdif,n) 
  !     ---------------------------------------------                     
  USE header
  !     use level n                                                       
  !     computes d2u/dx2 + d2u/dy2 at the cell centers.                   
  !                                                                       
  IMPLICIT NONE 
  integer i,j,k,n,nj2
  REAL(kind=rc_kind) :: udif(NI,NJ,NK),vdif(NI,NJ,NK),wdif(NI,NJ,NK),    &
       sdif(NI,NJ,NK),Tdif(NI,NJ,NK),Trdif(ntr,NI,NJ,NK)                            
  REAL(kind=rc_kind) ::                                                  &
       dudxfc(0:NI),dvdxfc(0:NI),dwdxfc(0:NI),dsdxfc(0:NI),dTdxfc(0:NI),dTrdxfc(ntr,0:NI), &
   &   dudyfc(0:NJ),dvdyfc(0:NJ),dwdyfc(0:NJ),dsdyfc(0:NJ),dTdyfc(0:NJ),dTrdyfc(ntr,0:NJ)
  REAL(kind=rc_kind) :: vyj,uxi,fac,ubylsq,wbylsq 
  REAL(kind=rc_kind) :: Kx_l,Ky_l(NJ),e1inv,yfrac,y0,y0inv,yend,yy,KxPassTr                                                     
                                                                        
!     constant eddy viscosity (in m2/s)                                 
!                                                                       
!      Kx= 1.d-1                                                        
!      KxPassTr= 1.d-1                                                  
  Kx_l= Kx 
  KxPassTr= Kx 
                                                                        
  do j=1,NJ 
     Ky_l(j) = Ky 
  end do
                                                                        
  do k=1,NK 
    do i=1,NI 
      do j=1,NJ-1 
        vyj= 0.5*(vy(i,j)+vy(i,j+1)) 
        dudyfc(j)= vyj*(u(i,j+1,k,n)-u(i,j,k,n)) 
        dvdyfc(j)= vyj*(v(i,j+1,k,n)-v(i,j,k,n)) 
        dwdyfc(j)= vyj*(w(i,j+1,k,n)-w(i,j,k,n)) 
        dsdyfc(j)= vyj*(s(i,j+1,k,n)-s(i,j,k,n)) 
        dTdyfc(j)= vyj*(T(i,j+1,k,n)-T(i,j,k,n)) 
        do it=1,ntr 
           dTrdyfc(it,j)= vyj*(Tr(it,i,j+1,k,n)-Tr(it,i,j,k,n)) 
        end do
      end do
      do j=0,NJ,NJ 
         dudyfc(j)= 0.d0 
         dvdyfc(j)= 0.d0 
         dwdyfc(j)= 0.d0 
         dsdyfc(j)= 0.d0 
         dTdyfc(j)= 0.d0 
         do it=1,ntr 
            dTrdyfc(it,j)= 0.d0 
         end do
      end do
      do j=1,NJ 
        udif(i,j,k)= Ky_l(j)*vy(i,j)*(dudyfc(j)-dudyfc(j-1)) 
        vdif(i,j,k)= Ky_l(j)*vy(i,j)*(dvdyfc(j)-dvdyfc(j-1)) 
        wdif(i,j,k)= Ky_l(j)*vy(i,j)*(dwdyfc(j)-dwdyfc(j-1)) 
        sdif(i,j,k)= Ky_l(j)*vy(i,j)*(dsdyfc(j)-dsdyfc(j-1)) 
        Tdif(i,j,k)= Ky_l(j)*vy(i,j)*(dTdyfc(j)-dTdyfc(j-1)) 
        do it=1,ntr 
           Trdif(it,i,j,k)= KxPassTr*vy(i,j)*(dTrdyfc(it,j)-dTrdyfc(it,j-1))   
        end do
      end do
    end do
  end do
                                                                        
  do k=1,NK 
    do j=1,NJ 
      do i=1,NI-1 
        uxi= 0.5*(ux(i,j)+ux(i+1,j)) 
        dudxfc(i)= uxi*(u(i+1,j,k,n)-u(i,j,k,n)) 
        dvdxfc(i)= uxi*(v(i+1,j,k,n)-v(i,j,k,n)) 
        dwdxfc(i)= uxi*(w(i+1,j,k,n)-w(i,j,k,n)) 
        dsdxfc(i)= uxi*(s(i+1,j,k,n)-s(i,j,k,n)) 
        dTdxfc(i)= uxi*(T(i+1,j,k,n)-T(i,j,k,n)) 
        do it=1,ntr 
           dTrdxfc(it,i)= uxi*(Tr(it,i+1,j,k,n)-Tr(it,i,j,k,n)) 
        end do
      end do
      !     periodic-ew boundaries                                            
      do i=0,NI,NI 
        uxi= 0.5*(ux(NI,j)+ux(1,j)) 
        dudxfc(i)= uxi*(u(1,j,k,n)-u(NI,j,k,n)) 
        dvdxfc(i)= uxi*(v(1,j,k,n)-v(NI,j,k,n)) 
        dwdxfc(i)= uxi*(w(1,j,k,n)-w(NI,j,k,n)) 
        dsdxfc(i)= uxi*(s(1,j,k,n)-s(NI,j,k,n)) 
        dTdxfc(i)= uxi*(T(1,j,k,n)-T(NI,j,k,n)) 
        do it=1,ntr 
           dTrdxfc(it,i)= uxi*(Tr(it,1,j,k,n)-Tr(it,NI,j,k,n)) 
        end do
      end do
      do i=1,NI 
        udif(i,j,k)= udif(i,j,k)+ Kx_l*ux(i,j)*(dudxfc(i)-dudxfc(i-1))                
        vdif(i,j,k)= vdif(i,j,k)+ Kx_l*ux(i,j)*(dvdxfc(i)-dvdxfc(i-1))                
        wdif(i,j,k)= wdif(i,j,k)+ Kx_l*ux(i,j)*(dwdxfc(i)-dwdxfc(i-1))                
        sdif(i,j,k)= sdif(i,j,k)+ Kx_l*ux(i,j)*(dsdxfc(i)-dsdxfc(i-1))                
        Tdif(i,j,k)= Tdif(i,j,k)+ Kx_l*ux(i,j)*(dTdxfc(i)-dTdxfc(i-1))                
        do it=1,ntr 
           Trdif(it,i,j,k)= Trdif(it,i,j,k)+ KxPassTr*ux(i,j)*(dTrdxfc(it,i)-dTrdxfc(it,i-1))
        end do
      end do
    end do
  end do
                                                                        
  !     Save to fricu,fricv,fricw before multiplying by Jac*fac           
  !     multiply by UL/(LEN*LEN) or WL/(LEN*LEN) to put fricu,fricv,fricw 
  !     in dimensional form                                               
                                                                        
  fac= 1.0/(UL*LEN) 
  ubylsq= UL/(LEN*LEN) 
  wbylsq= WL/(LEN*LEN) 
  do k=1,NK 
    do j=1,NJ 
      do i=1,NI 
        fricu(i,j,k)= fricu(i,j,k)+ udif(i,j,k)*ubylsq 
        fricv(i,j,k)= fricw(i,j,k)+ vdif(i,j,k)*ubylsq 
        fricw(i,j,k)= fricw(i,j,k)+ wdif(i,j,k)*wbylsq 
                                                                     
        udif(i,j,k)= fac*Jac(i,j,k)*udif(i,j,k) 
        vdif(i,j,k)= fac*Jac(i,j,k)*vdif(i,j,k) 
        wdif(i,j,k)= fac*Jac(i,j,k)*wdif(i,j,k) 
        sdif(i,j,k)= fac*Jac(i,j,k)*sdif(i,j,k) 
        Tdif(i,j,k)= fac*Jac(i,j,k)*Tdif(i,j,k) 
        do it=1,ntr 
           TRdif(it,i,j,k)= fac*Jac(i,j,k)*TRdif(it,i,j,k) 
        end do
      end do
    end do
  end do
                                                                        
  return 
END subroutine 
