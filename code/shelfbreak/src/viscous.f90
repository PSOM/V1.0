subroutine viscous(var,vardif) 
  !     ---------------------------------------------                     
  USE header
  !     use level n                                                       
  !     computes d2u/dx2 + d2u/dy2 at the cell centers.                   
  !                                                                       
  IMPLICIT NONE 
  integer i,j,k,n,nj2
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1) :: var
  REAL(kind=rc_kind) :: vardif(NI,NJ,NK)
  REAL(kind=rc_kind) :: dvardxfc(0:NI), dvardyfc(0:NJ)
  REAL(kind=rc_kind) :: vyj,uxi,fac,ubylsq,wbylsq 
  REAL(kind=rc_kind) :: Kx_l,Ky_l(NJ),e1inv,yfrac,y0,y0inv,yend,yy,KxPassTr                                                     

  integer :: selvar
                                                                        
  Kx_l= Kx; KxPassTr= Kx; 
                                                                        
  Ky_l(1:NJ) = Ky 
                                                                       
  do k=1,NK 
    do i=1,NI 
      do j=1,NJ-1 
        dvardyfc(j)= 0.5*(vy(i,j)+vy(i,j+1))*(var(i,j+1,k)-var(i,j,k)) 
      end do
      dvardyfc(0)= 0.d0 
      dvardyfc(NJ)= 0.d0 
      do j=1,NJ 
        vardif(i,j,k)= Ky_l(j)*vy(i,j)*(dvardyfc(j)-dvardyfc(j-1)) 
      end do
    end do
  end do
                                                                        
  do k=1,NK 
    do j=1,NJ 
      do i=1,NI-1 
        dvardxfc(i)= 0.5*(ux(i,j)+ux(i+1,j))*(var(i+1,j,k)-var(i,j,k)) 
      end do
      !     periodic-ew boundaries                                            
      dvardxfc(0)= 0.5*(ux(NI,j)+ux(1,j))*(var(1,j,k)-var(NI,j,k)) 
      dvardxfc(NI)= 0.5*(ux(NI,j)+ux(1,j))*(var(1,j,k)-var(NI,j,k)) 
      do i=1,NI 
        vardif(i,j,k)= vardif(i,j,k)+ Kx_l*ux(i,j)*(dvardxfc(i)-dvardxfc(i-1))                
      end do
    end do
  end do
                                                                        
                                                                        
  fac= 1.0/(UL*LEN) 

  vardif(1:NI,1:NJ,1:NK)=fac*Jac(1:NI,1:NJ,1:NK)*vardif(1:NI,1:NJ,1:NK)

                                                                        
  return 
END subroutine viscous
