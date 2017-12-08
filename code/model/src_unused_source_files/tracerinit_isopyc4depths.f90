      subroutine tracerinit(step)
!     --------------------                                              
!     initializes tracer fields                                         
!                                                                       
      USE header
      integer  i,j,k,l,it,step,ktr2,ktr3,ktr4,jrelease 
      REAL(kind=rc_kind) :: zt(ntr),zb(ntr),zmm(ntr),dsal 
!                                                                       
      if (step.gt.1) then 
         write(6,*) 'in tracerinit, trinit needs to be set correctly' 
         stop 
      end if 
                                                                        
!     if this is a continuation of a run, only trinit is to be initializ
!     Tr 1 is based on density.                                         
                                                                        
                                                                        
!     TRACER 1                                                          
!     =========                                                         
      sigrelease(1)=0
      deprelease(1)=0
      it=1 
      do k=0,NK+1 
         do j=0,NJ+1 
            do i=0,NI+1 
               if (rho(i,j,k).lt.1025.2) then   !USUAL - GOOD
!               if (rho(i,j,k).lt.1025.3) then   ! for byby5
                  T(it,i,j,k,0)= 1.d0 
               else 
                  T(it,i,j,k,0)= 0.d0 
               end if 
            end do 
         end do 
      end do 

                                                                        
!     Release tracers 2,3,4,5 streaks at 4 depths                         
      ktr2=11    ! at 216m
      ktr3 =9    ! at 257 m
      ktr4=7     ! at 303 m
      ktr5=5     ! at 355 m
!     shallow ML                                                        
!      ktr2= 23                                                         
!      ktr3=17                                                          
!      ktr4=10                                                          
!     TRACER 2,3,4                                                      
!     =============================                                     
      do k=0,NK+1 
         do j=0,NJ+1 
            do i=0,NI+1 
               do it=2,ntr
                  T(it,i,j,k,0) = 0.d0 
               end do 
            end do 
         end do 
      end do 


!     Release at center in a streak
!      jrelease=(NJ/2)- 0.08333333*NJ 
      jrelease=(NJ+1)/2

      if (step.ne.ntrstart) goto 101

!     Release with gradient in y, tr=1 at j=1, tr=0 at j=NJ
      do i=0,NI+1 
         do j=1,NJ
            T(2,i,j,ktr2,0)= 1.0 - dble(j-1)/dble(NJ)
            T(3,i,j,ktr3,0)= 1.d0 - dble(j-1)/dble(NJ)
            T(4,i,j,ktr4,0)= 1.d0 - dble(j-1)/dble(NJ)
            T(5,i,j,ktr5,0)= 1.d0 - dble(j-1)/dble(NJ)
         end do
         T(6,i,jrelease,ktr2,0)= 1.0
         T(7,i,jrelease,ktr3,0)= 1.d0
         T(8,i,jrelease,ktr4,0)= 1.d0
         T(9,i,jrelease,ktr5,0)= 1.d0
      end do

101   i=NI/2 
      j=NJ/2
      sigrelease(2)=rho(i,j,ktr2) 
      sigrelease(3)=rho(i,j,ktr3) 
      sigrelease(4)=rho(i,j,ktr4) 
      sigrelease(5)=rho(i,j,ktr5) 
      sigrelease(6)= sigrelease(2) 
      sigrelease(7)= sigrelease(3) 
      sigrelease(8)= sigrelease(4)
      sigrelease(9)= sigrelease(5)

      deprelease(2)= zc(NI/2,NJ/2,ktr2)*1000.
      deprelease(3)= zc(NI/2,NJ/2,ktr3)*1000.
      deprelease(4)= zc(NI/2,NJ/2,ktr4)*1000.
      deprelease(5)= zc(NI/2,NJ/2,ktr5)*1000.
      deprelease(6)= deprelease(2) 
      deprelease(7)= deprelease(3) 
      deprelease(8)= deprelease(4) 
      deprelease(9)= deprelease(5) 
                                                                        
      write(6,*) 'sigrelease ', sigrelease 
                                                                        
!     No gradient across solid boundaries
      do n=0,1
         do k=0,NK+1
            do i=1,NI
               do it=1,ntr
                  t(it,i,0,k,n)=t(it,i,1,k,n)
                  t(it,i,NJ+1,k,n)=t(it,i,NJ,k,n)
               end do
            end do
         end do
      end do
 
!     Periodic Boundaries (periodic e-w)                                
      do n=0,1 
         do k=0,NK+1 
            do j=0,NJ+1 
               do it=1,ntr 
                  t(it,0,j,k,n)= t(it,NI,j,k,n) 
                  t(it,NI+1,j,k,n)= t(it,1,j,k,n) 
               end do 
            end do 
         end do 
      end do 

                                                                        
                                                                        
      return 
      END                                           
