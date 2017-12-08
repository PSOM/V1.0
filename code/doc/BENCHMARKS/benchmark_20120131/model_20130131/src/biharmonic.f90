subroutine biharmonic(udif,vdif,m,step) 
!     ---------------------------------------------                     
  USE header
!     use level m                                                       
!     computes d4(u,v)/dx4, d4(u,v)/dy4  at the cell centers.           
!     The biharmonic operator d4/dx4 is applied along sigma levels      
!     We assume that vx(i,j)=0, uy(i,j)=0, i.e. the grid is rectilinear.
!     Biharmonic diffusion  for u,v                                     
!     The biharmonic viscosity bi_nu, used by Chapman and Lentz is 5.d9m
!     and is constant in space and time.                                
!     We extend the u,v, arrays by 2 points outside the boundary        
                                                                        
!                                                                       
      integer i,j,k,m,step 
      REAL(kind=rc_kind) :: binu, fac,z,seval 
      REAL(kind=rc_kind) :: dep0(NK),dep1(NK),uvert0(NK),uvert1(NK),         &
     &     bu0(NK),cu0(NK),du0(NK),bu1(NK),cu1(NK),du1(NK),             &
     &     vvert0(NK),vvert1(NK),                                       &
     &     bv0(NK),cv0(NK),dv0(NK),bv1(NK),cv1(NK),dv1(NK)              
      REAL(kind=rc_kind) :: udif(NI,NJ,NK), vdif(NI,NJ,NK) 
      REAL(kind=rc_kind) :: ui(-1:NI+2),vi(-1:NI+2),                         &
     &     uxi(-1:NI+1),uxi2(0:NI+1),uxi3(0:NI),uxi4(1:NI),             &
     &     vxi(-1:NI+1),vxi2(0:NI+1),vxi3(0:NI),vxi4(1:NI),             &
     &     u4x(NI)                                                      
      REAL(kind=rc_kind) :: uj(-1:NJ+2),vj(-1:NJ+2),                         &
     &     ueta(-1:NJ+1),ueta2(0:NJ+1),ueta3(0:NJ),ueta4(1:NJ),         &
     &     veta(-1:NJ+1),veta2(0:NJ+1),veta3(0:NJ),veta4(1:NJ),         &
     &     v4y(NJ)                                                      
      parameter (binu=5.d7) 
      fac= binu/(UL*LEN*LEN*LEN) 
                                                                        
!     i-th direction                                                    
!     ----------------                                                  
      do j=1,NJ 
         do k=1,NK 
            dep0(k)= zc(1,j,k) 
            uvert0(k)= u(1,j,k,m) 
            vvert0(k)= v(1,j,k,m) 
            dep1(k)= zc(NI,j,k) 
            uvert1(k)= u(NI,j,k,m) 
            vvert1(k)= v(NI,j,k,m) 
         end do 
         call spline(NK,dep0,uvert0,bu0,cu0,du0) 
         call spline(NK,dep1,uvert1,bu1,cu1,du1) 
!                                                                       
         call spline(NK,dep0,vvert0,bv0,cv0,dv0) 
         call spline(NK,dep1,vvert1,bv1,cv1,dv1) 
!                                                                       
         do i=1,NI 
            u4x(i)= ux(i,j)*ux(i,j)*ux(i,j)*ux(i,j) 
         end do 
                                                                        
         do k=1,NK 
            if (abs(zc(0,j,k)-zc(1,j,k)).lt.1.d-16) then 
               ui(0)= u(1,j,k,m) 
               vi(0)= v(1,j,k,m) 
               ui(-1)= u(1,j,k,m) 
               vi(-1)= v(1,j,k,m) 
            else 
               ui(0)= seval(NK,zc(0,j,k),dep0,uvert0,bu0,cu0,du0) 
               vi(0)= seval(NK,zc(0,j,k),dep0,vvert0,bv0,cv0,dv0) 
               z= 2.d0*zc(0,j,k) -zc(1,j,k) 
               ui(-1)= seval(NK,z,dep0,uvert0,bu0,cu0,du0) 
               vi(-1)= seval(NK,z,dep0,vvert0,bv0,cv0,dv0) 
            end if 
!                                                                       
            if (abs(zc(NI+1,j,k)-zc(NI,j,k)).lt.1.d-16) then 
               ui(NI+1)= u(NI,j,k,m) 
               vi(NI+1)= v(NI,j,k,m) 
               ui(NI+2)= u(NI,j,k,m) 
               vi(NI+2)= v(NI,j,k,m) 
            else 
               ui(NI+1)= seval(NK,zc(NI+1,j,k),dep1,uvert1,bu1,cu1,du1) 
               vi(NI+1)= seval(NK,zc(NI+1,j,k),dep1,vvert1,bv1,cv1,dv1) 
               z= 2.d0*zc(NI+1,j,k) -zc(NI,j,k) 
               ui(NI+2)= seval(NK,z,dep1,uvert1,bu1,cu1,du1) 
               vi(NI+2)= seval(NK,z,dep1,vvert1,bv1,cv1,dv1) 
            end if 
!                                                                       
            do i=1,NI 
               ui(i)= u(i,j,k,m) 
               vi(i)= v(i,j,k,m) 
            end do 
!                                                                       
            goto 101 
!     Now compute du/dxi,d2u/dxi2,d3u/dxi3,d4u/dxi4                     
!     at cell faces...                                                  
            do i=-1,NI+1 
               uxi(i)= ui(i+1) -ui(i) 
               vxi(i)= vi(i+1) -vi(i) 
            end do 
!     at cell centers                                                   
            do i=0,NI+1 
               uxi2(i)= uxi(i) -uxi(i-1) 
               vxi2(i)= vxi(i) -vxi(i-1) 
            end do 
!     at cell faces                                                     
            do i=0,NI 
               uxi3(i)= uxi2(i+1) -uxi2(i) 
               vxi3(i)= vxi2(i+1) -vxi2(i) 
            end do 
!     at cell centers                                                   
            do i=1,NI 
               uxi4(i)= uxi3(i) -uxi3(i-1) 
               vxi4(i)= vxi3(i) -vxi3(i-1) 
            end do 
  101       do i=1,NI 
               udif(i,j,k)= udif(i,j,k) +                               &
     &              (ui(i+2) -4.0*ui(i+1) +6.0*ui(i)                    &
     &              -4.0*ui(i-1) + ui(i-2))*u4x(i)*fac*Jac(i,j,k)       
               vdif(i,j,k)= vdif(i,j,k) +                               &
     &              (vi(i+2) -4.0*vi(i+1) +6.0*vi(i)                    &
     &              -4.0*vi(i-1) + vi(i-2))*u4x(i)*fac*Jac(i,j,k)       
!=               udif(i,j,k)= udif(i,j,k) + u4x(i)*fac*Jac(i,j,k)*uxi4(i
!=               vdif(i,j,k)= vdif(i,j,k) + u4x(i)*fac*Jac(i,j,k)*vxi4(i
!              udif(i,j,k)=  u4x(i)*fac*uxi4(i)                         
!              vdif(i,j,k)=  u4x(i)*fac*vxi4(i)                         
            end do 
         end do 
      end do 
                                                                        
!      if (step.eq.2) then                                              
!      i=25                                                             
!      do k=1,NK                                                        
!         write(250,*) 'k=',k                                           
!         do j=1,NJ                                                     
!            do i=1,NI                                                  
!            write(250,*) j,udif(i,j,k),vdif(i,j,k)                     
!            end do                                                     
!         end do                                                        
!      end do                                                           
!      end if                                                           
!                                                                       
!     j-th direction                                                    
!     ----------------                                                  
      do i=1,NI 
         do k=1,NK 
            dep0(k)= zc(i,1,k) 
            uvert0(k)= u(i,1,k,m) 
            vvert0(k)= v(i,1,k,m) 
            dep1(k)= zc(i,NJ,k) 
            uvert1(k)= u(i,NJ,k,m) 
            vvert1(k)= v(i,NJ,k,m) 
         end do 
         call spline(NK,dep0,uvert0,bu0,cu0,du0) 
         call spline(NK,dep1,uvert1,bu1,cu1,du1) 
!                                                                       
         call spline(NK,dep0,vvert0,bv0,cv0,dv0) 
         call spline(NK,dep1,vvert1,bv1,cv1,dv1) 
!                                                                       
         do j=1,NJ 
            v4y(j)= vy(i,j)*vy(i,j)*vy(i,j)*vy(i,j) 
         end do 
!                                                                       
         do k=1,NK 
            if (abs(zc(i,1,k)-zc(i,0,k)).lt.1.d-16) then 
               uj(0)= u(i,1,k,m) 
               vj(0)= v(i,1,k,m) 
               uj(-1)= u(i,1,k,m) 
               vj(-1)= v(i,1,k,m) 
            else 
               uj(0)= seval(NK,zc(i,0,k),dep0,uvert0,bu0,cu0,du0) 
               vj(0)= seval(NK,zc(i,0,k),dep0,vvert0,bv0,cv0,dv0) 
               z= 2.d0*zc(i,0,k) -zc(i,1,k) 
               uj(-1)= seval(NK,z,dep0,uvert0,bu0,cu0,du0) 
               vj(-1)= seval(NK,z,dep0,vvert0,bv0,cv0,dv0) 
            end if 
!                                                                       
            if (abs(zc(i,NJ,k)-zc(i,NJ+1,k)).lt.1.d-16) then 
               uj(NJ+1)= u(i,NJ,k,m) 
               vj(NJ+1)= v(i,NJ,k,m) 
               uj(NJ+2)= u(i,NJ,k,m) 
               vj(NJ+2)= v(i,NJ,k,m) 
            else 
               uj(NJ+1)= seval(NK,zc(i,NJ+1,k),dep1,uvert1,bu1,cu1,du1) 
               vj(NJ+1)= seval(NK,zc(i,NJ+1,k),dep1,vvert1,bv1,cv1,dv1) 
               z= 2.d0*zc(i,NJ+1,k) -zc(i,NJ,k) 
               uj(NJ+2)= seval(NK,z,dep1,uvert1,bu1,cu1,du1) 
               vj(NJ+2)= seval(NK,z,dep1,vvert1,bv1,cv1,dv1) 
            end if 
!                                                                       
            do j=1,NJ 
               uj(j)= u(i,j,k,m) 
               vj(j)= v(i,j,k,m) 
            end do 
                                                                        
                                                                        
            goto 201 
!     at cell faces...                                                  
            do j=-1,NJ+1 
               ueta(j)= uj(j+1) -uj(j) 
               veta(j)= vj(j+1) -vj(j) 
            end do 
!     at cell centers                                                   
            do j=0,NJ+1 
               ueta2(j)= ueta(j) -ueta(j-1) 
               veta2(j)= veta(j) -veta(j-1) 
            end do 
!     at cell faces                                                     
            do j=0,NJ 
               ueta3(j)= ueta2(j+1) -ueta2(j) 
               veta3(j)= veta2(j+1) -veta2(j) 
            end do 
!     at cell centers                                                   
            do j=1,NJ 
               ueta4(j)= ueta3(j) -ueta3(j-1) 
               veta4(j)= veta3(j) -veta3(j-1) 
            end do 
  201       do j=1,NJ 
               udif(i,j,k)= udif(i,j,k) +                               &
     &              (uj(j+2) -4.0*uj(j+1) +6.0*uj(j)                    &
     &              -4.0*uj(j-1) + uj(j-2))*v4y(j)*fac*Jac(i,j,k)       
               vdif(i,j,k)= vdif(i,j,k) +                               &
     &              (vj(j+2) -4.0*vj(j+1) +6.0*vj(j)                    &
     &              -4.0*vj(j-1) + vj(j-2))*v4y(j)*fac*Jac(i,j,k)       
!=               udif(i,j,k)= udif(i,j,k) +v4y(j)*fac*Jac(i,j,k)*ueta4(j
!=               vdif(i,j,k)= vdif(i,j,k) +v4y(j)*fac*Jac(i,j,k)*veta4(j
!              udif(i,j,k)=  v4y(j)*fac*ueta4(j)                        
!              vdif(i,j,k)=  v4y(j)*fac*veta4(j)                        
!-               if ((step.eq.2).and.(i.eq.25)) then                    
!-                  write(60,*) v4y(j)*fac,ueta4(j),veta4(j)            
!-                  write(70,*) vj(j),veta(j),veta2(j),veta3(j)         
!-               end if                                                 
            end do 
         end do 
      end do 
!      if (step.eq.2) then                                              
!      i=25                                                             
!      do k=1,NK                                                        
!         write(350,*) 'k=',k                                           
!         do j=1,NJ                                                     
!c            do i=1,NI                                                 
!            write(350,*) j,udif(i,j,k),vdif(i,j,k)                     
!c            end do                                                    
!         end do                                                        
!      end do                                                           
!      stop                                                             
!      endif                                                            
                                                                        
      return 
      END                                           
