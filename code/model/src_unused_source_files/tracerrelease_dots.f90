subroutine tracerrelease(it,jrelease) 
  !     --------------------                                              
  USE header
  !     initializes tracer fields                                         
  !                                                                       
      implicit REAL(kind=rc_kind) :: (a-h,o-z) 
      integer  i,j,k,l,it,step,ktr2,ktr3,ktr4,jrelease 
      REAL(kind=rc_kind) :: zt(ntr),zb(ntr),zmm(ntr),dsal 
!                                                                       
                                                                        
!     Release tracers 4,5,6 very short streaks at 256 m depth           
!     Release at center                                                 
!     -------------------                                               
      write(6,*) 'sigrelease',sigrelease(it),'at j=', jrelease 
!      jrelease=(NJ/2)- 0.08333333*NJ                                   
      j=jrelease 
!     FIND sigma surface                                                
                                                                        
!      do i=1,NI                                                        
      do i=NI/2-1,NI/2+1 
         do k=2,NK 
            if (rho(i,j,k).le.sigrelease(it)) then 
               dsig= rho(i,j,k-1)-rho(i,j,k) 
               sigup= (rho(i,j,k-1) - sigrelease(it))/dsig 
               sigdn= (sigrelease(it)-rho(i,j,k))/dsig 
               T(it,i,j,k,0)= sigup 
               T(it,i,j,k-1,0)= sigdn 
               goto 101 
            end if 
         end do 
!     If NK reached and rho above not less than sig, then sig is at surf
         T(it,i,j,NK,n)=1.d0 
  101 continue 
      end do 
                                                                        
                                                                        
!     Periodic Boundaries (periodic e-w)                                
!      n=0                                                              
!      do k=0,NK+1                                                      
!         do j=0,NJ+1                                                   
!            do it=1,ntr                                                
!               t(it,0,j,k,n)= t(it,NI,j,k,n)                           
!               t(it,NI+1,j,k,n)= t(it,1,j,k,n)                         
!            end do                                                     
!         end do                                                        
!      end do                                                           
                                                                        
                                                                        
      return 
      END                                           
