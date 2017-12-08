subroutine coriolis(n) 
!----------------------------------------------------                   
 USE header
  !     Evaluate the Coriolis terms si,sj                                 
 integer i,j,k,n 
 REAL(kind=rc_kind) :: pxi,peta,psig,px,py,fac,fac2,ainv,betainv 
                                                                   
 fac= EPS*delta 
 fac2= EPS*lambda 
 ainv= 1.d0/apr 
 betainv= 1.0/beta 
                                                                   
 !  We are using the values at the ghost points.                      
 do k=1,NK 
   do j=1,NJ 
     do i=1,NI 
       si(i,j,k)= -ffc(i,j)*v(i,j,k,n) +fac*bbc(i,j)*w(i,j,k,n) + drpx(i,j,k) - EPS* uvis(i,j,k)                                   
       sj(i,j,k)=  ffc(i,j)*u(i,j,k,n)                          + drpy(i,j,k) - EPS* vvis(i,j,k)                                  
       sk(i,j,k)=  fnhhy*(-bbc(i,j)*u(i,j,k,n) - fac2*(u(i,j,k,n)*u(i,j,k,n)+v(i,j,k,n)*v(i,j,k,n))*ainv) - betainv*wvis(i,j,k) 
     enddo
   enddo
  enddo
 return 
END                                           
