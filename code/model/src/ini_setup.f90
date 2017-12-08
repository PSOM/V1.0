SUBROUTINE ini_setup(pcorr)  

#include "cppdefs.h"
  !     subroutine init(p,vfent,ufex)                                    
  !     ------------------------------------------------                  
  !     For the curvilinear grid.                                         
  !     sigma levels are evenly spaced.  Hence wz is a function of        
  !     time but not of z.   Here u,v,w refer to the xi,eta,sigma direcito
  !     Those metric quantities that do not change with time are evaluated
  !     The rest are evaluated in "sigma.f" which will be called at every 
  !     step.                                                             
  !     Also it is absolutely essential that the integral of the flux     
  !     over the entrance section is equal to that over the exit section. 
  !     -No longer necessary with the free surface.                       
  !     This may require some special attention when the free surface is  
  !     moving at the exit. At the entrance vf is held fixed.             
  !                                                                       
 USE header 

 REAL(kind=rc_kind) ::   xdu(0:NI+1,0:NJ+1),ydu(0:NI+1,0:NJ+1),xdv(0:NI+1,0:NJ+1),ydv(0:NI+1,0:NJ+1)

 REAL(kind=rc_kind) :: pcorr(maxout)

 REAL(kind=rc_kind) :: fconst,phi0,cosphi0,sinphi0,dthet,dtheta,dphi,dep,y
 
 INTEGER i,j,k,jmid,iseed 

 DLinv= 1.d0/DL 
 phi0= phi0deg*PI/180.d0  !     phi0 is the central lat, dphi,dtheta grid spacing(angle)          
 dphi = dy/(apr*AL) 



 !--------------------------------------------------
 ! INITIALIZATION OF THE GRID
 !


 ! HORIZONTAL GRID

 xc(0) = -0.5*dx*1.d-3 
 yc(0) = -0.5*dy*1.d-3; 

 do i=1,NI+1; xc(i)= xc(i-1) + dx*1.d-3; end do;
 do j=1,NJ+1; yc(j)= yc(j-1) + dy*1.d-3; end do; 

 DO  j=0,NJ+1 
   DO  i=0,NI+1 
     xdu(i,j)= dx/LEN 
     ydu(i,j)= 0.d0 
     xdv(i,j)= 0.d0 
     ydv(i,j)= dy/LEN !- used for constant dy
     ! ydv(i,j)= dyM(j)/LEN 
   ENDDO
 ENDDO

 DO j=0,NJ+1 
   DO i=0,NI+1 
     J2d(i,j)=  xdu(i,j)*ydv(i,j) - xdv(i,j)*ydu(i,j) 
     ux(i,j) =  ydv(i,j)/J2d(i,j) 
     vx(i,j) = -ydu(i,j)/J2d(i,j) 
     uy(i,j) = -xdv(i,j)/J2d(i,j) 
     vy(i,j) =  xdu(i,j)/J2d(i,j) 
     g11(i,j)=  ux(i,j)*ux(i,j) + uy(i,j)*uy(i,j) 
     g12(i,j)=  ux(i,j)*vx(i,j) + uy(i,j)*vy(i,j) 
     g22(i,j)=  vx(i,j)*vx(i,j) + vy(i,j)*vy(i,j) 
   ENDDO
 ENDDO
    
    
 ! DEPTH OF EVERY MODEL COLUMN

 if(.NOT.(lv_flat_bottom)) then
   call ini_topog
  else
   dep = total_depth
   D(:,:)= -dep*DLinv   ! D(i,j) is -ve and  non-dim by DL
   CALL smooth          ! Compute partial derivatives dD/dx,dD/dy 
                        ! Does not smooth anything.   
 endif

 !----------------------------------------------
 ! COMPUTATION OF THE PLANETARY VORTICITY.
 !

 ! fplane is read from namelist.
 ! 1: f-plane
 ! 0: f has a latitude-dependance.

 fconst= 2.d0*OMEGA/FPAR 
 jmid= NJ/2 
 DO  j=0,NJ+1 
   latrad(j)= phi0 + REAL(1-fplane) * ( DBLE(j-jmid)*dphi ) 
   !latrad(j)= phi0                                               
 ENDDO

 DO j=1,NJ
   ffi(:,j)= fconst*sin(latrad(j)) 
   bbi(:,j)= fnhhy*fconst*cos(latrad(j)) 
 ENDDO
 DO j=0,NJ
   ffj(:,j)= fconst*sin(0.5d0*(latrad(j+1)+latrad(j))) 
   bbj(:,j)= fnhhy*fconst*cos(0.5d0*(latrad(j+1)+latrad(j))) 
 ENDDO
 DO j=0,NJ+1
   ffc(:,j)= fconst*sin(latrad(j)) 
   bbc(:,j)= fnhhy*fconst*cos(latrad(j)) 
 ENDDO



 !--------------------------------------------------
 ! INITIALIZATION OF THE VARIABLES
 !

 u(:,:,:,0) = 0d0
 v(:,:,:,0) = 0d0
 w(:,:,:,0) = 0d0
 s(:,:,:,0) = 0d0
 T(:,:,:,0) = 0d0
 Tr(:,:,:,:,0) = 0d0
 conv(:,:,:)    = 0
 con100(:,:,:)  = 0

 pcorr(:)= 0.d0 

 uvis(:,:,:) = 0.d0
 vvis(:,:,:) = 0.d0
 wvis(:,:,:) = 0.d0

 h = 0d0
 uf = 0d0
 vf = 0d0
 wf = 0d0

 ufbce(:,:)= uf(NI,:,:) 
 ufbcw(:,:)= uf(0,:,:) 


 vfbcn(:,:)= vf(:,NJ,:) 
 vfbcs(:,:)= vf(:,0,:) 

 wfbcb(:,:) = 0.d0    ! at the sea bed, wf=0                                              

 !---------------------

 CALL findzall        ! finds the vertical grid
                         
 CALL ini_st          ! initialize s,T                                                    
 CALL evalrho(rho,0)  ! deduce rho from s,T 

 CALL ini_h           ! initialize free surface h

 CALL findzall        ! find z again


#ifdef file_output                     
 if(lv_flat_bottom) then                 
   OPEN (unit=60,file=TRIM(dirout)//'zgrid.out')  !     write out z-grid                                                  
     WRITE(60,*) '# vertical grid' 
     do k=0,NK+1; WRITE(60,*) k,zc(10,10,k)*1000.; end do;
     WRITE(60,*) '# face values' 
     do k=-1,NK+1;WRITE(60,*) k,zf(10,10,k)*1000.; end do;
   CLOSE(60) 
 endif
#endif

 advecpv(:)= 0.d0 
 friction(:)= 0.d0 
 diabatic(:)= 0.d0 

 RETURN 
END SUBROUTINE ini_setup
