 subroutine cdiv(dtimel,maxdiv,n) 
!     -------------------------------                                   
 USE header
!     checks divergence using u,v,w at cell centers                     
 integer i,j,k,n,imax,jmax,kmax 
 REAL(kind=rc_kind) :: Udxi,Vdeta,Wdsig,div,maxdiv,dtimel 
!      edt= EPS/dtime                                                   
!     We have now changed the defn of tol so that it is just (u_x +v_y +
!     We do not need edd any more.                                      

 maxdiv= 0.d0 
 do k=1,NK 
   do j=1,NJ 
     do i=1,NI 
       Udxi= 0.5*( (u(i+1,j,k,n)*ux(i+1,j)   + v(i+1,j,k,n)*uy(i+1,j))*Jac(i+1,j,k) & 
               & - (u(i-1,j,k,n)*ux(i-1,j)   + v(i-1,j,k,n)*uy(i-1,j))*Jac(i-1,j,k) )             
       Vdeta= 0.5*( (u(i,j+1,k,n)*vx(i,j+1)   + v(i,j+1,k,n)*vy(i,j+1))*Jac(i,j+1,k) & 
                & - (u(i,j-1,k,n)*vx(i,j-1)   + v(i,j-1,k,n)*vy(i,j-1))*Jac(i,j-1,k) )             
       Wdsig= 0.5*( (u(i,j,k+1,n)*wx(i,j,k+1) + v(i,j,k+1,n)*wy(i,j,k+1) + EPS*w(i,j,k+1,n)*wz(i,j,k+1))*Jac(i,j,k+1) &
                & - (u(i,j,k-1,n)*wx(i,j,k-1) + v(i,j,k-1,n)*wy(i,j,k-1) + EPS*w(i,j,k-1,n)*wz(i,j,k-1))*Jac(i,j,k-1) )  
       div= abs(Udxi+ Vdeta +Wdsig) 
       if (div.gt.maxdiv) then 
          maxdiv=div; imax=i; jmax=j; kmax=k; 
          end if
       enddo
     enddo
  enddo 
  return 
  END                                           
