subroutine vort(n) 
  !     ---------------------------------------------------               
  USE header,only : NI,NJ,NK,u,v,w,ux,vx,wx,uy,vy,wy,wz,vor,shear,strain,UL,LEN,DL,DLinv,rc_kind
  !     computes the vorticity (for now just the vertical component)      

      integer i,j,k,n 
      REAL(kind=rc_kind) :: vdx,udy,udx,vdy,LENinv 
      REAL(kind=rc_kind) ::  vz,uz,udz,vdz
                                                                       

! ----------------------------------------------------------------------------

  !  j=0,NJ                                                            
  do i=1,NI 
     do k=1,NK 
        u(i,0,k,n)= 2.0*u(i,1,k,n) -u(i,2,k,n) 
        v(i,0,k,n)= 2.0*v(i,1,k,n) -v(i,2,k,n) 
        w(i,0,k,n)= 2.0*w(i,1,k,n) -w(i,2,k,n) 
                                                                    
        u(i,NJ+1,k,n)= 2.0*u(i,NJ,k,n) -u(i,NJ-1,k,n) 
        v(i,NJ+1,k,n)= 2.0*v(i,NJ,k,n) -v(i,NJ-1,k,n) 
        w(i,NJ+1,k,n)= 2.0*w(i,NJ,k,n) -w(i,NJ-1,k,n) 
     end do 
  end do 


  !  k=0,NK+1                                                            
  do j=0,NJ+1 
     do i=1,NI 
        u(i,j,0,n)= 2.0*u(i,j,1,n) -u(i,j,2,n) 
        v(i,j,0,n)= 2.0*v(i,j,1,n) -v(i,j,2,n) 
        w(i,j,0,n)= 2.0*w(i,j,1,n) -w(i,j,2,n) 
                                                                    
        u(i,j,NK+1,n)= 2.0*u(i,j,NK,n) -u(i,j,NK-1,n) 
        v(i,j,NK+1,n)= 2.0*v(i,j,NK,n) -v(i,j,NK-1,n) 
        w(i,j,NK+1,n)= 2.0*w(i,j,NK,n) -w(i,j,NK-1,n) 
     end do 
  end do 

  !     periodic-ew boundaries                                            
  do k=0,NK+1 
     do j=0,NJ+1 
        u(0,j,k,n)= u(NI,j,k,n) 
        v(0,j,k,n)= v(NI,j,k,n) 
        w(0,j,k,n)= w(NI,j,k,n) 
                                                                    
        u(NI+1,j,k,n)= u(1,j,k,n) 
        v(NI+1,j,k,n)= v(1,j,k,n) 
        w(NI+1,j,k,n)= w(1,j,k,n) 
     end do 
  end do 
             
! ----------------------------------------------------------------------------

                                                       
  DLinv= 1.0/DL 
  LENinv=1.0/LEN 

  ! Computation of vorticity and strain
  do j=1,NJ 
    do i=1,NI 
      do k=1,NK 
        vdx= 0.5*( (v(i+1,j,k,0) - v(i-1,j,k,0)) * ux(i,j) +     &
             &     (v(i,j+1,k,0) - v(i,j-1,k,0)) * vx(i,j) +     &
             &     (v(i,j,k+1,0) - v(i,j,k-1,0)) * wx(i,j,k) )            
        udy= 0.5*( (u(i+1,j,k,0) - u(i-1,j,k,0)) * uy(i,j) +     &
             &     (u(i,j+1,k,0) - u(i,j-1,k,0)) * vy(i,j) +     &
             &     (u(i,j,k+1,0) - u(i,j,k-1,0)) * wy(i,j,k) )            
        vor(i,j,k)= (vdx -udy)*UL*LENinv 
                                                                 
        udx= 0.5*( (u(i+1,j,k,0) - u(i-1,j,k,0)) * ux(i,j) +     &
             &     (u(i,j+1,k,0) - u(i,j-1,k,0)) * vx(i,j) +     &
             &     (u(i,j,k+1,0) - u(i,j,k-1,0)) * wx(i,j,k) )          
        vdy= 0.5*( (v(i+1,j,k,0) - v(i-1,j,k,0)) * uy(i,j) +     &
             &     (v(i,j+1,k,0) - v(i,j-1,k,0)) * vy(i,j) +     &
             &     (v(i,j,k+1,0) - v(i,j,k-1,0)) * wy(i,j,k) )          
        strain(i,j,k)= sqrt((udx - vdy)**2 +(vdx+ udy)**2) * UL * LENinv                                          
      end do 
    end do 
  end do


! ----------------------------------------------------------------------------


  ! Computation of the vertical shear                                                    
  do j=1,NJ 
    do i=1,NI 
      k=1 
      vdz= (v(i,j,k+1,0)-v(i,j,k,0))*wz(i,j,k) 
      udz= (u(i,j,k+1,0)-u(i,j,k,0))*wz(i,j,k) 
      shear(i,j,k)= sqrt(udz*udz + vdz*vdz)*UL/DL 
      do k=2,NK-1 
        vdz = 0.5*( (v(i,j,k+1,0) -v(i,j,k-1,0))*wz(i,j,k) ) 
        udz = 0.5*( (u(i,j,k+1,0) -u(i,j,k-1,0))*wz(i,j,k) ) 
        shear(i,j,k)= sqrt(udz*udz + vdz*vdz)*UL/DL 
      end do 
      k=NK 
      vdz= (v(i,j,k,0)-v(i,j,k-1,0))*wz(i,j,k) 
      udz= (u(i,j,k,0)-u(i,j,k-1,0))*wz(i,j,k) 
      shear(i,j,k)= sqrt(udz*udz + vdz*vdz)*UL/DL 
    end do 
  end do 

                                                                    
! ----------------------------------------------------------------------------

                                                                    
  !     boundary points                                                   
  !     j=0,NJ                                                            
  do i=1,NI 
    do k=1,NK 
      vor(i,0,k)= 2.0*vor(i,1,k) -vor(i,2,k) 
      vor(i,NJ+1,k)= 2.0*vor(i,NJ,k) -vor(i,NJ-1,k) 
      shear(i,0,k)= 2.0*shear(i,1,k) -shear(i,2,k) 
      shear(i,NJ+1,k)= 2.0*shear(i,NJ,k) -vor(i,NJ-1,k) 
    end do 
  end do 

  do j=0,NJ+1 
    do i=1,NI 
      vor(i,j,0)= 2.0*vor(i,j,1) -vor(i,j,2) 
      vor(i,j,NK+1)= 2.0*vor(i,j,NK) -vor(i,j,NK-1) 
      shear(i,j,0)= 2.0*shear(i,j,1) -shear(i,j,2) 
      shear(i,j,NK+1)= 2.0*shear(i,j,NK) -shear(i,j,NK-1) 
    end do 
  end do 

  !     periodic-ew boundaries                                            
  do k=0,NK+1 
    do j=0,NJ+1 
      vor(0,j,k)= vor(NI,j,k) 
      vor(NI+1,j,k)= vor(1,j,k) 
      shear(0,j,k)= shear(NI,j,k) 
      shear(NI+1,j,k)= shear(1,j,k) 
    end do 
  end do 
                                                                    
  return 
  END                                           
