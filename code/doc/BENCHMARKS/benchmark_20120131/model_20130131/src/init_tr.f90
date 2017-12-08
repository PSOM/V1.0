SUBROUTINE init(pcorr)  

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
!  use grids

  REAL(kind=rc_kind) ::   xdu(0:NI+1,0:NJ+1),ydu(0:NI+1,0:NJ+1),xdv(0:NI+1,0:NJ+1),ydv(0:NI+1,0:NJ+1)

  REAL(kind=rc_kind) :: pcorr(maxout)

  REAL(kind=rc_kind) :: C,fconst,phi0,cosphi0,sinphi0,sumuf,sumvf,dthet,dtheta,dphi,          &
       &     gin,c1,c2,temp,temp2,dnkinv,tmplon,tmplat,cjent,dep,y
 

  ! added 20120302
  REAL(kind=rc_kind) ::  dscl, rdscl, yarg,ex2y,thy,sechy2,ran3, hlfac

  !                                                                       
  !     FPAR is the value used for scaling Coriolis parameters f and b    
  INTEGER i,j,k,jmid,iseed 

  DLinv= 1.d0/DL 
  C= PI/180.d0 
  gin= 1.d0/gpr 
  dnkinv= 1.d0/DBLE(NK) 
                                                                       
  phi0= phi0deg*C  !     phi0 is the central lat, dphi,dtheta grid spacing(angle)          
  !print*, "# phi0deg = ",phi0deg

  dphi = dy/(apr*AL) 
  !print*, '# dphi = ',dphi, 'in deg ',dphi/C 

  jmid= NJ/2 
  !print*,'# jmid=',jmid 

  ! fplane is read from namelist.
  ! 1: fplane
  ! 0: f has a latitude-dependance.
  DO  j=0,NJ+1 
     latrad(j)= phi0 + REAL(1-fplane) * ( DBLE(j-jmid)*dphi ) 
     !latrad(j)= phi0                                               
  enddo

  !read dy information from dyM.data   
  !     call load_dy(dyM)
  ! print*, dyM

  DO  j=0,NJ+1 
    DO  i=0,NI+1 
      xdu(i,j)= dx/LEN 
      ydu(i,j)= 0.d0 
      xdv(i,j)= 0.d0 
      ydv(i,j)= dy/LEN !- used for constant dy
      ! ydv(i,j)= dyM(j)/LEN 
    ENDDO
  ENDDO

  !---- the following block is used for varying dy -moved to mod_grids.f90 -J.Wang 21-Dec-2011
  !     yc(0) = -0.5*dyM(0)*1.d-3 
  !     DO j=1,NJ+1 
  !        yc(j)= yc(j-1) + dyM(j)*1.d-3 
  !     END DO
  !---- the following block is used for constant dy
  !     yc(0) = -0.5*dy*1.d-3 
  !     DO j=1,NJ+1 
  !        yc(j)= yc(j-1) + dy*1.d-3 
  !     END DO

  xc(0) = -0.5*dx*1.d-3 
  DO i=1,NI+1 
    xc(i)= xc(i-1) + dx*1.d-3 
  END DO

  yc(0) = -0.5*dy*1.d-3 
  do j=1,NJ+1 
     yc(j)= yc(j-1) + dy*1.d-3 
  end do 



  DO  j=0,NJ+1 
    DO  i=0,NI+1 
      J2d(i,j)= xdu(i,j)*ydv(i,j) -xdv(i,j)*ydu(i,j) 
      ux(i,j)= ydv(i,j)/J2d(i,j) 
      vx(i,j)= -ydu(i,j)/J2d(i,j) 
      uy(i,j)= -xdv(i,j)/J2d(i,j) 
      vy(i,j)= xdu(i,j)/J2d(i,j) 
      g11(i,j)= ux(i,j)*ux(i,j) +uy(i,j)*uy(i,j) 
      g12(i,j)= ux(i,j)*vx(i,j) +uy(i,j)*vy(i,j) 
      g22(i,j)= vx(i,j)*vx(i,j) +vy(i,j)*vy(i,j) 
    ENDDO
  ENDDO
        

   if(.NOT.(lv_flat_bottom)) then
       call topog
      else
       !     D(i,j) is -ve and  non-dim by DL    
       dep = total_depth
       D(:,:)= -dep*DLinv
       !     Compute partial  dD/dx,dD/dy    
       CALL smooth 
     endif


     ! 101  call findzall                                                    
     !101  CONTINUE
     !      close(55)                                                        
     !c      open (unit= 31, file='xyD.dat')                                 
     !c      do j=0,NJ+1                                                     
     !c         tmplat= phi0deg + dble(j-jmid)*dthet                         
     !c         yc(j)= tmplat                                                
     !c         do i=0,NI+1                                                  
     !c            tmplon= 300.d0 +dthet*dble(i)                             
     !c            write(31,*) tmplon,tmplat,D(i,j)                          
     !c            if (j.eq.0) xc(i)= tmplon                                 
     !c         end do                                                       
     !c      end do                                                          
     !c      close(31)                                                       
     !     In smooth we smooth D and then use central differencing to evaluat
     !     Ddx and Ddy. The values obtained are much lower than those obtaine
     !     by smoothing Ddx, Ddy evaluated from the bspline - probably bec   
     !     we had a bug with a factor of 3 which is now corrected. (values   
     !     differ even after the correction and we use the numerial values.) 
     !c      call smooth                                                     
     !*    flat bottom                                                       
     ! 999  do j=0,NJ+1                                                      
     !         do i=0,NI+1                                                   
     !            Ddx(i,j)= 0.d0                                             
     !            Ddy(i,j)= 0.d0                                             
     !            D(i,j)= 1.d0                                               
     !            write(150,*) Ddx(i,j)                                      
     !            write(250,*) Ddy(i,j)                                      
     !            write(350,*) D(i,j)                                        
     !         end do                                                        
     !      end do                                                           
     !      stop                                                             
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
     !                                                                       
     !     specify the initial fields                                        
     !     read in the cell-centered pressure and u,v velocities             
     h = 0d0
     uf = 0d0
     vf = 0d0
     wf = 0d0
     !     fill in the vfbc arrays                                           

     ufbce(:,:)= uf(NI,:,:) 
     ufbcw(:,:)= uf(0,:,:) 


     vfbcn(:,:)= vf(:,NJ,:) 
     vfbcs(:,:)= vf(:,0,:) 

     !     at the sea bed, wf=0                                              

     wfbcb(:,:) = 0.d0 
     !     Evaluate the Coriolis parameter f' and fill the arrays ffi,ffj,ffc
     !     f= f0+By= 2Omega(Sin(phi0) +(phi-phi0)Cos(phi0))                  
     !     b= b0+By= 2Omega(Cos(phi0) -(phi-phi0)Sin(phi0))                  
     fconst= 2.d0*OMEGA/FPAR 

    
     DO j =1, NJ
        ffi(:,j)= fconst*sin(latrad(j)) 
        bbi(:,j)= fnhhy*fconst*cos(latrad(j)) 
     ENDDO
     DO j =0, NJ
        ffj(:,j)= fconst*sin(0.5d0*(latrad(j+1)+latrad(j))) 
        bbj(:,j)= fnhhy*fconst*cos(0.5d0*(latrad(j+1)+latrad(j))) 
     ENDDO
     DO j=0, NJ+1
        ffc(:,j)= fconst*sin(latrad(j)) 
        bbc(:,j)= fnhhy*fconst*cos(latrad(j)) 
     ENDDO
     !                                                                       
     !     Initialize s,T                                                    
     CALL findzall   ! finds the vertical grid
     !print*,zc(10,10,:)
     !print*, 'findzall'
!     CALL ini_rho()
     call stprofile

     call evalrho(rho,0) 

     CALL inith           ! initialize free surface h
     !print*, 'hmean',hmean
     CALL findzall        ! find z again
     !     write out z-grid                                                  
     !     -----------------           
#ifdef file_output                     
     if(lv_flat_bottom) then                 
       OPEN (unit=60,file=TRIM(dirout)//'zgrid.out')
       WRITE(60,*) '# vertical grid' 
       DO k=0,NK+1 
          WRITE(60,*) k,zc(10,10,k)*1000. 
       END DO
       WRITE(60,*) '# face values' 
       DO k=-1,NK+1 
          WRITE(60,*) k,zf(10,10,k)*1000. 
       END DO
       CLOSE(60) 
    endif
#endif
     !     -----------                                                       


     advecpv(:)= 0.d0 
     friction(:)= 0.d0 
     diabatic(:)= 0.d0 

!if (pickup_step == 0) then
     !     Perturb the front boundary if the model does not start from pickup                                     
     !     for a quicker onset of instability                                
     !     ----------------------------------------------------              
!     iseed= 44294 
!     dum = ran3(iseed) 
!     !                                                                       
!     DO i=0,NI+1 
!        perturb =  1.d-3*ran3(iseed)
!        DO j=1,NJ 
!           h(i,j)= h(i,j) +perturb 
!        END DO
!     END DO
!endif

     RETURN 

   END SUBROUTINE init
