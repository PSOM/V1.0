subroutine uvchy(dtimel) 
!----------------------------------------------------                   
  USE header
!  !USE rpgrads
!     compute the tilde velocities at cell centers and store them in n=1
  integer i,j,k 
  REAL(kind=rc_kind) :: dte,dtimel,hxi,heta,hx,hy,kaph1,                  &
  &     Ub,Vb,Wb,d2,c1,c2,epsinv,temp1,temp2,wfk                     
!     operation:  u(n)= u(0) + dtime* f(u(m))                           
!                                                                       
  REAL(kind=rc_kind) :: pxi,peta,psig,px,py 
                                                                       
  dte= dtimel/EPS 
  epsinv= 1.d0/EPS 
  kaph1= 1.d0 - kappah 
                                                                   
  do j=1,NJ 
    do i=1,NI 
      hxi= 0.5d0*( h(i+1,j)-h(i-1,j) ) 
      heta= 0.5d0*( h(i,j+1)-h(i,j-1) ) 
      hx= ux(i,j)*hxi +vx(i,j)*heta 
      hy= uy(i,j)*hxi +vy(i,j)*heta 
      do k=1,NK 
        pxi= 0.5d0*(p(i+1,j,k)-p(i-1,j,k)) 
        !     IF STATEMENTS ADDED jUN 8,2005 (doesn't make much diff)           
        if (j.eq.1) then 
           peta= p(i,j+1,k)-p(i,j,k) 
        else if (j.eq.NJ) then 
           peta= p(i,j,k)-p(i,j-1,k) 
        else 
           peta= 0.5d0*(p(i,j+1,k)-p(i,j-1,k)) 
        endif 
        psig= 0.5d0*(p(i,j,k+1)-p(i,j,k-1)) 
        px= ux(i,j)*pxi +vx(i,j)*peta +wx(i,j,k)*psig 
        py= uy(i,j)*pxi +vy(i,j)*peta +wy(i,j,k)*psig 
                                                                 
        !     cx and cy contain the convective terms.           
        cx(i,j,k)=  cx(i,j,k) -dte*( gpr*(kappah*hx+ kaph1*gradhn(i,j,1)) + si(i,j,k) + qpr *px ) 
        cy(i,j,k)=  cy(i,j,k) -dte*( gpr*(kappah*hy+ kaph1*gradhn(i,j,2)) + sj(i,j,k) + qpr *py ) 
       !cz(i,j,k)=  cz(i,j,k)                                   
      enddo
     gradhn(i,j,1)= hx 
     gradhn(i,j,2)= hy 
   enddo
 enddo

  !     czf is already computed                                           
  do j=1,NJ 
    do i=1,NI 
      do k=1,NK 
        wfk= 0.5d0*(czf(i,j,k) +czf(i,j,k-1)) 
        cz(i,j,k)= (wfk/Jac(i,j,k) -cx(i,j,k)*wx(i,j,k) -cy(i,j,k)*wy(i,j,k) )/(EPS*wz(i,j,k))              
      enddo
    enddo
  enddo




 return 



 END                                           
