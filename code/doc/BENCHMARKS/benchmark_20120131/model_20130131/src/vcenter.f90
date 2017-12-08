  subroutine vcenter(pf,dtimel,n) 
!----------------------------------------------------                   
  USE header
!     compute the final vel (n+1 th step) at the cell centers           
!     operation:  u(n)= u(0) + dtime* f(u(m))                           
  integer i,j,k,n 
  REAL(kind=rc_kind) :: dtbeta,dtimel,pxi,peta,psig,c1,c2,dte,  &
   &     Ub,Vb,Wb,d2,Wr,Wtemp(0:NI+1,0:NJ+1),                   &
   &     px,py,pz,fac,pf(0:NI+1,0:NJ+1,0:NK+1),wfk                    
  REAL(kind=rc_kind) :: z,seval,bu0(NK),cu0(NK),du0(NK),        &
   &     bv0(NK),cv0(NK),dv0(NK),bw0(NK),cw0(NK),dw0(NK),       &
   &     dep(NK),uvert(NK),vvert(NK),wvert(NK)                        

  dte= dtimel/EPS 
  dtbeta= dtimel*beta 
                                                                     
  do j=1,NJ 
    do i=1,NI 
      do k=1,NK 
        pxi= 0.5d0*(pf(i+1,j,k)-pf(i-1,j,k)) 
        peta= 0.5d0*(pf(i,j+1,k)-pf(i,j-1,k)) 
        psig= 0.5d0*(pf(i,j,k+1)-pf(i,j,k-1)) 
        px= ux(i,j)*pxi +vx(i,j)*peta +wx(i,j,k)*psig 
        py= uy(i,j)*pxi +vy(i,j)*peta +wy(i,j,k)*psig 
        pz= wz(i,j,k)*psig 

        u(i,j,k,n)=  cx(i,j,k) -dte*(qpr*px +si(i,j,k) ) 
        v(i,j,k,n)=  cy(i,j,k) -dte*(qpr*py +sj(i,j,k) ) 
       !w(i,j,k,n)=  cz(i,j,k) -dtbeta*( pz +sk(i,j,k) )       

      enddo
    enddo
  enddo 

  !  vface is already called - compute w from wf                       
  do j=1,NJ 
    do i=1,NI 
      do k=1,NK 
         wfk= 0.5d0*(wf(i,j,k) +wf(i,j,k-1)) 
         w(i,j,k,n)= (wfk/Jac(i,j,k) -u(i,j,k,n)*wx(i,j,k)-v(i,j,k,n)*wy(i,j,k) )/(EPS*wz(i,j,k))             
      enddo
    enddo
  enddo
  
  return 

  END                                           
