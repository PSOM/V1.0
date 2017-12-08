subroutine calc_dvdy(dvdy,n) 
  !     ---------------------------------------------------               
  USE header
  !     computes the vorticity (for now just the vertical component)      
      integer i,j,k,n 
      REAL(kind=rc_kind) :: vdy,LENinv, dvdy(0:NI+1,0:NJ+1,0:NK+1) 
                                                                        
!     j=0,NJ                                                            
      do i=1,NI 
         do k=1,NK 
            v(i,0,k,n)= 2.0*v(i,1,k,n) -v(i,2,k,n) 
!                                                                       
            v(i,NJ+1,k,n)= 2.0*v(i,NJ,k,n) -v(i,NJ-1,k,n) 
         end do 
      end do 
      do j=0,NJ+1 
         do i=1,NI 
            v(i,j,0,n)= 2.0*v(i,j,1,n) -v(i,j,2,n) 
!                                                                       
            v(i,j,NK+1,n)= 2.0*v(i,j,NK,n) -v(i,j,NK-1,n) 
         end do 
      end do 
!     periodic-ew boundaries                                            
      do k=0,NK+1 
         do j=0,NJ+1 
            v(0,j,k,n)= v(NI,j,k,n) 
!                                                                       
            v(NI+1,j,k,n)= v(1,j,k,n) 
         end do 
      end do 
                                                                        
      DLinv= 1.0/DL 
      LENinv=1.0/LEN 
      do j=1,NJ 
         do i=1,NI 
            do k=1,NK 
               vdy= 0.5*((v(i+1,j,k,0)-v(i-1,j,k,0))*uy(i,j) +          &
     &              (v(i,j+1,k,0) -v(i,j-1,k,0))*vy(i,j) +              &
     &              (v(i,j,k+1,0) -v(i,j,k-1,0))*wy(i,j,k) )            
               dvdy(i,j,k)= (vdy)*UL*LENinv 
            end do 
         end do 
      end do 
                                                                        
!     boundary points                                                   
!     j=0,NJ                                                            
      do i=1,NI 
         do k=1,NK 
            dvdy(i,0,k)= 2.0*dvdy(i,1,k) -dvdy(i,2,k) 
            dvdy(i,NJ+1,k)= 2.0*dvdy(i,NJ,k) -dvdy(i,NJ-1,k) 
         end do 
      end do 
      do j=0,NJ+1 
         do i=1,NI 
            dvdy(i,j,0)= 2.0*dvdy(i,j,1) -dvdy(i,j,2) 
            dvdy(i,j,NK+1)= 2.0*dvdy(i,j,NK) -dvdy(i,j,NK-1) 
         end do 
      end do 
!     periodic-ew boundaries                                            
      do k=0,NK+1 
         do j=0,NJ+1 
            dvdy(0,j,k)= dvdy(NI,j,k) 
            dvdy(NI+1,j,k)= dvdy(1,j,k) 
         end do 
      end do 
                                                                        
      return 
      END                                           
