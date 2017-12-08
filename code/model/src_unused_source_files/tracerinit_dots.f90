subroutine tracerinit(step) 
!     --------------------                                              
  USE header
!     initializes tracer fields                                         
!                                                                       
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
      it=1 
      do k=0,NK+1 
         do j=0,NJ+1 
            do i=0,NI+1 
               if (rho(i,j,k).lt.1025.2) then 
                  T(it,i,j,k,0)= 1.d0 
               else 
                  T(it,i,j,k,0)= 0.d0 
               end if 
            end do 
         end do 
      end do 
                                                                        
!     Release tracers 2,3,4 streaks at 3 depths                         
!      ktr2=11     !is 98.6 m    use k=23, z=59 m                       
!      ktr3 =9    !is 150 m     use k=17  z=123 m                       
!      ktr4=7     !is 236.1 m  isopycn is flat                          
      ktr2=9 
      ktr3=9 
!     shallow ML                                                        
!      ktr2= 23                                                         
!      ktr3=17                                                          
!      ktr4=10                                                          
!     TRACER 2,3,4                                                      
!     =============================                                     
      do k=0,NK+1 
         do j=0,NJ+1 
            do i=0,NI+1 
               do it=2,6 
!               do it=2,4                                               
                  T(it,i,j,k,0) = 0.d0 
               end do 
            end do 
         end do 
      end do 
!     Release at center                                                 
      jrelease=(NJ/2)- 0.08333333*NJ 
      write(6,*) 'jrelease =', jrelease 
      j=jrelease 
      do i=NI/2-1,NI/2+1 
         T(2,i,j-10,ktr2,0)= 1.d0 
         T(3,i,j+10,ktr3,0)= 1.d0 
!         T(4,i,j,ktr4,0)= 1.d0                                         
      end do 
      i=NI/2 
      sigrelease(2)=rho(i,j,ktr2) 
      sigrelease(3)=rho(i,j,ktr3) 
      sigrelease(4)= sigrelease(3) 
      sigrelease(5)= sigrelease(3) 
      sigrelease(6)= sigrelease(3) 
                                                                        
      write(6,*) 'sigrelease ', sigrelease 
                                                                        
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
