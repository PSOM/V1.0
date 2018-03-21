subroutine calcskfc 
  !----------------------------------------------------                   
  USE header
  !     after h is computed                                               
                                                                        
!     We interpolate the source terms onto the cell faces               
!     It is important that u,v,rp have been filled outside the boundarie
!      implicit REAL(kind=rc_kind) :: (a-h,o-z) 
  integer i,j,k 
  REAL(kind=rc_kind) :: be2,Jack,kaph1 
  REAL(kind=rc_kind) :: hxi, heta, gradh, hy, py,px 
  REAL(kind=rc_kind) ::  hx, wxk, wyk, wxsk, wysk, wzsk                                                                        
  !     We are using the values at the ghost points.                      
                                                                    
  be2= beta*EPS*EPS 
  kaph1= 1.d0 -kappah 
                                                                    
  do i=1,NI 
    do  j=1,NJ 

      hxi= 0.5d0*( h(i+1,j)-h(i-1,j) ) 
      heta= 0.5d0*( h(i,j+1)-h(i,j-1) ) 

      hx= ux(i,j)*hxi +vx(i,j)*heta 
      hy= uy(i,j)*hxi +vy(i,j)*heta 

      hx = gpr*(kappah*hx + kaph1*gradhn(i,j,1)) 
      hy = gpr*(kappah*hy + kaph1*gradhn(i,j,2)) 

      do  k=1,NK-1 

        wxk = 0.5d0*(wx(i,j,k+1) +wx(i,j,k)) 
        wxsk= wxk*(0.5d0*(si(i,j,k+1)+si(i,j,k)) + hx) 

        wyk = 0.5d0*(wy(i,j,k+1) +wy(i,j,k)) 
        wysk= wyk*(0.5d0*(sj(i,j,k+1)+sj(i,j,k)) + hy) 

        wzsk= wzk(i,j,k)*(sk(i,j,k+1)+ sk(i,j,k) )*0.5 

        Jack= 0.5d0*(Jac(i,j,k+1) + Jac(i,j,k)) 
        skfc(i,j,k)= ( be2*wzsk +wxsk +wysk )*Jack 

      end do 

      k=0 
      ! linear extrapolation for wzsk                                     
      wzsk= wzk(i,j,k)*0.5*(3.0*sk(i,j,k+1)-sk(i,j,k+2)) 
      wxsk= wx(i,j,k+1)*(si(i,j,k+1) +hx) 
      wysk= wy(i,j,k+1)*(sj(i,j,k+1) +hy) 

      Jack= Jac(i,j,k+1) 
      skfc(i,j,k)= ( be2*wzsk +wxsk +wysk )*Jack 
                                                                  
      k= NK 
      ! linear extrap                                                     
      wzsk= wzk(i,j,k)*0.5*(3.0*sk(i,j,k)-sk(i,j,k-1)) 
      wxsk= 0.5d0*(wx(i,j,k+1) +wx(i,j,k))*(si(i,j,k)+hx) 
      wysk= 0.5d0*(wy(i,j,k+1) +wy(i,j,k))*(sj(i,j,k)+hy) 

      Jack= Jac(i,j,k) 
      skfc(i,j,k)= ( be2*wzsk +wxsk +wysk )*Jack 

    end do 
  end do 
  return 
 END                                           
