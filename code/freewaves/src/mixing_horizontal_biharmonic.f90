subroutine mixing_horizontal_biharmonic(var,vardif) 
!     ---------------------------------------------                     
  USE header
!     use level m                                                       
!     computes d4(var)/dx4, d4(var)/dy4  at the cell centers.           
!     The biharmonic operator d4/dx4 is applied along sigma levels      
!     We assume that vx(i,j)=0, uy(i,j)=0, i.e. the grid is rectilinear.
!     Biharmonic diffusion  for var                                     
!     The biharmonic viscosity bi_nu, used by Chapman and Lentz is 5.d9m
!     and is constant in space and time.                                

      integer i,j,k,m,step
!       REAL(kind=rc_kind) :: facx,facy
      REAL(kind=rc_kind), dimension(0:NI+1,0:NJ+1,0:NK+1) :: var
      REAL(kind=rc_kind) :: vardif(NI,NJ,NK),facx(NI,NJ,NK),facy(NI,NJ,NK)
      REAL(kind=rc_kind) :: vari(-1:NI+2),u4x(NI)                                                      
      REAL(kind=rc_kind) :: varj(-1:NJ+2),v4y(NJ)                                                      
      INTEGER OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS,CHUNK,NPROC

   if(iv_compute_kh==0) then
      Kx(:,:,:)= B0x
      Ky(:,:,:)= B0y
   else
      Kx(:,:,:)= B0x*Kz(:,:,:)/KzMom
      Ky(:,:,:)= B0y*Kz(:,:,:)/KzMom
   endif

   facx(:,:,:)= -Kx(1:NI,1:NJ,1:NK)/(UL*LEN*LEN*LEN) 
   facy(:,:,:)= -Ky(1:NI,1:NJ,1:NK)/(UL*LEN*LEN*LEN) 
      
      vardif(:,:,:) = 0.d0
!     j-th direction                                                    
!     ----------------                                                  
!$OMP PARALLEL DO PRIVATE(k,i,j,v4y,varj) SHARED(vardif) COLLAPSE(2) 
      do i=1,NI                                                                                                                 
         do j=1,NJ 
            v4y(j)= vy(i,j)*vy(i,j)*vy(i,j)*vy(i,j) 
         end do 

         do k=1,NK 
            do j=0,NJ+1 
               varj(j)= var(i,j,k) 
            end do 

!             do j=2,5
!                vardif(i,j,k)= (varj(j+2) -4.0*varj(j+1) +6.0*varj(j) &
!      &              -4.0*varj(j-1) + varj(j-2))*v4y(j)*facy*Jac(i,j,k)*100
!             end do
!             do j=6,NJ-5
!                vardif(i,j,k)= (varj(j+2) -4.0*varj(j+1) +6.0*varj(j) &
!      &              -4.0*varj(j-1) + varj(j-2))*v4y(j)*facy*Jac(i,j,k)
!             end do
!             do j=NJ-4,NJ
!                vardif(i,j,k)= (varj(j+2) -4.0*varj(j+1) +6.0*varj(j) &
!      &              -4.0*varj(j-1) + varj(j-2))*v4y(j)*facy*Jac(i,j,k)*100
!             end do

            do j=2,NJ-1
               vardif(i,j,k)= (varj(j+2) -4.0*varj(j+1) +6.0*varj(j) &
     &              -4.0*varj(j-1) + varj(j-2))*v4y(j)*facy(i,j,k)*Jac(i,j,k)                                                   
            end do 
            
            vardif(i,1, k) = 0.d0
            vardif(i,NJ,k) = 0.d0
         end do 
      end do 
!$OMP END PARALLEL DO

!     i-th direction                                                    
!     ----------------                                                  
!$OMP PARALLEL DO PRIVATE(k,j,i,u4x,vari) SHARED(vardif) COLLAPSE(2)                                                                       	                                           
      do j=1,NJ                                                                      
         do i=1,NI 
            u4x(i)= ux(i,j)*ux(i,j)*ux(i,j)*ux(i,j) 
         end do 

         do k=1,NK 
            do i=0,NI+1 
               vari(i)= var(i,j,k) 
            end do 
            vari(0)   = var(NI  ,j,k)
            vari(-1)  = var(NI-1,j,k) 
            vari(NI+2)= var(2,   j,k)
            vari(NI+1)= var(1,   j,k) 

            do i=1,NI 
               vardif(i,j,k)= vardif(i,j,k) + (vari(i+2) -4.0*vari(i+1) +6.0*vari(i) &
     &              -4.0*vari(i-1) + vari(i-2))*u4x(i)*facx(i,j,k)*Jac(i,j,k)                       
            end do 
         end do 
      end do 
!$OMP END PARALLEL DO

      return 
      END