  subroutine ini_uv(n) 
!     ---------------------------------------------                     
  USE header
!  USE rpgrads
!     modified for periodicew bc                                        
!     ---------------------------                                       
!     Sets up the initial velocities so that they are in geostrophic bal
!     with the initial density field.                                   
!      IMPLICIT NONE                                           
      integer i,j,k,n,iter,imax,jmax,kmax,m 
      REAL(kind=rc_kind) :: uxi,vyj,hxi,heta,hx,hy,px,py,ujfc,pyjfc 
      REAL(kind=rc_kind) :: res,resmax,fiter,fzero 
      REAL(kind=rc_kind) :: ainv,be2,fac2,wzsk,wxsk,wysk,Jack,pxi,peta,      &
     &     psig,pgrad,con,pz                                            
      REAL(kind=rc_kind) ::  kaph1
      ! pyjfc added on 2012/02/14                                                                      


imax=0;jmax=0;kmax=0;

 if(.NOT.(use_Shchepetkin)) then ! calculates baroclinic press gradients
    CALL rpevalgrad_Song(n) 
   else
    CALL rpevalgrad_Song(n) 
    CALL rpevalgrad_Sche(n);
! #include "test_Sche.f90"
    drpx(1:NI,1:NJ,1:NK)=ru4_Sche(1:NI,1:NJ,1:NK);  drpy(1:NI,1:NJ,1:NK)=rv4_Sche(1:NI,1:NJ,1:NK); 
    grpifc(0:NI,1:NJ,1:NK)=ru2_Sche(0:NI,1:NJ,1:NK); grpjfc(1:NI,0:NJ,1:NK)=rv2_Sche(1:NI,0:NJ,1:NK);
 endif

      kaph1= 1.d0 -kappah 
      fzero= 0.d0 
      fiter= 0.0 
      iter = 1 
      con = 1.0 -qpr 
  101 resmax = 1.d-20 
      do j=1,NJ 
         do i=1,NI 
            hxi= 0.5d0*( h(i+1,j)-h(i-1,j) ) 
            heta= 0.5d0*( h(i,j+1)-h(i,j-1) ) 
            hx= ux(i,j)*hxi +vx(i,j)*heta 
            hy= uy(i,j)*hxi +vy(i,j)*heta 
            do k=1,NK 
               pxi= 0.5d0*(p(i+1,j,k)-p(i-1,j,k)) 
               peta= 0.5d0*(p(i,j+1,k)-p(i,j-1,k)) 
               psig= 0.5d0*(p(i,j,k+1)-p(i,j,k-1)) 
               px= (ux(i,j)*pxi +vx(i,j)*peta +wx(i,j,k)*psig) 
               py= (uy(i,j)*pxi +vy(i,j)*peta +wy(i,j,k)*psig) 
                                                                        
!               res = u(i,j,k,n) + (qpr*py +drpy(i,j,k)+gpr*hy)/ffc(i,j)
!c               res = ffc(i,j)*u(i,j,k,n)*con +                        
!c     &              (qpr*py +drpy(i,j,k)*con+gpr*hy)                  
               res = ffc(i,j)*u(i,j,k,n) +                              &
     &              (qpr*py +drpy(i,j,k)+gpr*hy)                        
               if (abs(res).gt.abs(resmax)) then 
                  resmax = res 
                  imax = i 
                  jmax = j 
                  kmax = k 
               end if 
  
!               u(i,j,k,n) = - (qpr*py +con*drpy(i,j,k)+gpr*hy)/        
!     &              (ffc(i,j)*con)                                     
!               v(i,j,k,n) = (qpr*px +con*drpx(i,j,k) +gpr*hx)/         
!     &              (ffc(i,j)*con)                                     
               u(i,j,k,n) = - (qpr*py +drpy(i,j,k)+gpr*hy)/             &
     &              (ffc(i,j))                                          

                  !IF(i==20 .AND. j==120 .AND. k==1) then
                  !  PRINT*,"WQ2",u(i,j,k,n),drpy(i,j,k),gpr*hy,ffc(i,j)
                  !ENDIF


               v(i,j,k,n) = (qpr*px +drpx(i,j,k) +gpr*hx)/              &
     &              (ffc(i,j))                                          
!               if ((i.eq.48).and.(k.eq.30)) then
!                  write(120,*) drpx(i,j,k),drpy(i,j,k)
!               end if
            end do 
         end do 
         do k=1,NK 
            u(0,j,k,n)= u(NI,j,k,n) 
            v(0,j,k,n)= v(NI,j,k,n) 
            u(NI+1,j,k,n) = u(1,j,k,n) 
            v(NI+1,j,k,n) = v(1,j,k,n) 
         end do 
      end do 

!      write(6,*) iter,'geostroph',imax,jmax,kmax,resmax 

!      write(6,*) 'stop in geostroph'
!      stop
                                                                        
!     COMPUTE NH PRESSURE                                               
      if (fnhhy.gt.0.9) then 
         call coriolis(n) 
         call srcface_nopy(n) 
         call calcskfc 
                                                                        
         fac2= EPS*lambda 
         be2= beta*EPS*EPS 
         ainv= 1.d0/apr 
         do j=1,NJ 
            do i=1,NI 
               do k= NK,0,-1 
!     pz= (p(i,j,k+1) -p(i,j,k))*gqk(i,j,k,3)                           
!     pz = -skfc                                                        
                  pgrad = fiter* (                                      &
     &              0.25*(p(i+1,j,k+1)                                  &
     &              +p(i+1,j,k)-p(i-1,j,k+1)-p(i-1,j,k))*gqk(i,j,k,1)   &
     &              +0.25*(p(i,j+1,k+1)+p(i,j+1,k)-p(i,j-1,k+1)         &
     &              -p(i,j-1,k))*gqk(i,j,k,2)  )                        
                                                                        
                  if (k.eq.NK) then 
!     Use the BC that p(NK)+p(NK+1) = 0.0  (p=0 at the free-surf)       
                     p(i,j,k) = 0.5*(pgrad + skfc(i,j,k))/gqk(i,j,k,3) 
                  else 
                     p(i,j,k) = (p(i,j,k+1)*gqk(i,j,k,3) + pgrad        &
     &                    + skfc(i,j,k))/gqk(i,j,k,3)                   
                  end if 
               end do 
            end do 
         end do 
         call mgpfill(fzero,p) 
!     filter                                                            
         do m=1,8 
            do k=0,NK+1 
               do j=1,NJ 
                  do i=0,NI+1 
                     p(i,j,k)= 0.5*p(i,j,k)+0.25*(p(i,j-1,k)            &
     &                    +p(i,j+1,k))                                  
                  end do 
               end do 
            end do 
         end do 
         goto 301 
!     -------------- p at boundaries (mgpfill not suitable)-----        
!     faces j=0,j=NJ                                                    
!     --------------                                                    
         j=0 
         do k=1,NK 
            do i=1,NI 
               vyj= 0.5*(vy(i,j)+vy(i,j+1)) 
!               ujfc=(15.0*u(i,1,k,n)-10.*u(i,2,k,n)                    
!     &              +3.0*u(i,3,k,n))/8.0                               
!               ujfc= u(i,1,k,n)                                        
               u(i,0,k,n)= 2.*u(i,1,k,n) - u(i,2,k,n) 
               ujfc= 0.5*(u(i,1,k,n)+u(i,0,k,n)) 
!               if ((i.eq.8).and.(k.eq.4)) write(6,*)                   
!     &              ujfc,u(i,1,k,n),u(i,2,k,n),u(i,3,k,n)              
               pyjfc= -(gpr*hyn(i,j,k) +grpjfc(i,j,k) +ujfc*            &
     &              Jjfc(i,j,k)*vyj*ffc(i,j) )                          
               p(i,j,k) = p(i,j+1,k) -pyjfc/gqj(i,j,k,2) 
               p(i,0,k)= 2.*p(i,1,k) -p(i,2,k) 
            end do 
         end do 
         j=NJ 
         do k=1,NK 
            do i=1,NI 
               vyj= 0.5*(vy(i,j)+vy(i,j+1)) 
!               ujfc=(15.0*u(i,NJ,k,n)-10.*u(i,NJ-1,k,n)                
!     &              +3.0*u(i,NJ-2,k,n))/8.0                            
!               ujfc=u(i,NJ,k,n)                                        
               u(i,NJ+1,k,n)= 2.*u(i,NJ,k,n) - u(i,NJ-1,k,n) 
               ujfc= 0.5*(u(i,NJ,k,n)+u(i,NJ+1,k,n)) 
!               if ((i.eq.8).and.(k.eq.2)) write(6,*)                   
!     &              ujfc,u(i,NJ,k,n),u(i,NJ-1,k,n),u(i,NJ-2,k,n)       
               pyjfc= -(gpr*hyn(i,j,k) +grpjfc(i,j,k) +ujfc*            &
     &              Jjfc(i,j,k)*vyj*ffc(i,j) )                          
               p(i,j+1,k) = p(i,j,k) + pyjfc/gqj(i,j,k,2) 
               p(i,NJ+1,k)= 2.*p(i,NJ,k) -p(i,NJ-1,k) 
            end do 
         end do 
         k=NK 
         do j=0,NJ+1 
            do i=1,NI 
               p(i,j,k+1)= - p(i,j,k) 
            end do 
         end do 
!     k= 0 (already done)                                               
!     periodic e-w boundaries                                           
         do k=0,NK+1 
            do j=1,NJ 
               p(NI+1,j,k)= p(1,j,k) 
               p(0,j,k)= p(NI,j,k) 
            end do 
         end do 
  301    continue 
!     --------------------------------                                  
         iter = iter +1 
         fiter = 1.0 
      end if 
!     end computation of NH pressure                                    
!====                                                                   
!      if (iter.eq.24) goto 202                                         
!==                                                                     
      if ((iter.lt.50).and.(fnhhy.ge.0.001).and.                        &
     &     (abs(resmax).ge.1.d-15))                                    &
     &     goto 101                                                     
!                                                                       
!     -------                                                           
!     filter                                                            
!      do m=1,3                                                         
!      do k=0,NK+1                                                      
!         do j=1,NJ                                                     
!            do i=0,NI+1                                                
!               p(i,j,k)= 0.5*p(i,j,k)+0.25*(p(i,j-1,k)+p(i,j+1,k))     
!            end do                                                     
!         end do                                                        
!      end do                                                           
!      end do                                                           
! 202  write(6,*) 'not stopping in geostroph'                           
!      stop                                                             
      do j=1,NJ 
         do i=0,NI 
!            vyj= 0.5*(vy(i,j+1)+vy(i,j))                               
            do k=1,NK 
            uf(i,j,k)= (ux(i+1,j)*u(i+1,j,k,n) +ux(i,j)*u(i,j,k,n)      &
     &              +uy(i+1,j)*v(i+1,j,k,n)+uy(i,j)*v(i,j,k,n))*0.5     &
     &             *Jifc(i,j,k)                                         
!               uf(i,j,k)= grpjfc(i,j,k)/vyj*ffj(i,j)*Jjfc(i,j,k)       
            end do 
         end do 
      end do 
      do j=0,NJ 
         do i=1,NI 
!            hxi= 0.25*(h(i+1,j+1) +h(i+1,j) -h(i-1,j+1) -h(i-1,j))     
!            heta= h(i,j+1) -h(i,j)                                     
!            uxi= 0.5*(ux(1,j)+ux(NI,j))                                
            do k=1,NK 
!               hy= gj(i,j,k,1)*hxi +gj(i,j,k,2)*heta                   
!               gradh= gpr*(kappah*hy + kaph1*hyn(i,j,k))               
!               vf(i,j,k)= -gradh -grpjfc(i,j,k)/vxj*ffi(i,j)*Jifc(i,j,k
               vf(i,j,k)= 0.5*Jjfc(i,j,k)*                              &
     &              (vx(i,j+1)*u(i,j+1,k,n)+                            &
     &              vx(i,j)*u(i,j,k,n) +vy(i,j+1)*v(i,j+1,k,n)+         &
     &              vy(i,j)*v(i,j,k,n))                                 
            end do 
         end do 
      end do 
                                                                        
      if (rect.eqv.(.false.)) then 
         write(6,*) 'need to modify compuation of face velocities in    &
     &        subroutine geostroph'                                     
         stop 
      endif 
                                                                        
      return 
      END                                           
