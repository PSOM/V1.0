subroutine  correctbc 
  !     ------------------------------------------------------------      
  USE header
!                                                                       
      integer i,j,k,nexit 
      REAL(kind=rc_kind) :: flowin,flowout,const 
                                                                        
!     Outflow :  correct and update                                     
!     ---------                                                         
      flowin= 0.d0 
      flowout= 0.d0 
      nexit= 0 
      do k=1,NK 
         do i=1,NI 
!            if (vfbcs(i,k).ge.0.d0) then                               
               flowin= flowin +vfbcs(i,k) 
!            else                                                       
!               flowout= flowout - vfbcs(i,k)                           
!               nexit= nexit +1                                         
!            endif                                                      
         end do 
      end do 
      do k=NK-1,NK 
         do i=1,NI 
!            if (vfbcn(i,k).le.0.d0) then                               
!               flowin= flowin -vfbcn(i,k)                              
!            else                                                       
               flowout= flowout +vfbcn(i,k) 
               nexit= nexit +1 
!            endif                                                      
         end do 
      end do 
!                                                                       
!     MODIFY OUTFLOW VELOCITIES - assumes rectilinear grid in computing 
      const= (flowin - flowout)/nexit 
      !write(6,*) 'flowin',flowin,'flowout',flowout 
      !write(6,*) 'nexit',nexit,' const ',const 
!     modify only vfbcn in top 2 layers                                 
!                                                                       
      do k=NK-1,NK 
         do i=1,NI 
!            if (vfbcn(i,k).gt.0.d0) then                               
               vfbcn(i,k)= vfbcn(i,k) +const 
               vnorth(i,k)= vfbcs(i,k)/(vy(i,NJ)*Jjfc(i,NJ,k)) 
!            end if                                                     
         end do 
      end do 
!                                                                       
      return 
      END                                           
