subroutine surfaceflux(dtimel,n) 
  !     ----------------------------                                      
  USE header
  integer i,j,k,n 
  REAL(kind=rc_kind) :: dtimel,sflux,temp 
                                                                        
      sflux= 0.d0 
      do k=1,NK 
         do i=1,NI 
            sflux= sflux +vfbcs(i,k)*(sbackgrnd - ssouth(i,k)) 
         end do 
      end do 
      sflux= sflux/dble(NI*NJ)*dtime 
                                                                        
      do j=1,NJ 
         do i=1,NI 
            s(i,j,k,n)= s(i,j,k,n) + sflux/Jac(i,j,NK) 
!            temp= sflux/Jac(i,j,NK)                                    
!            write(100,*) temp                                          
         end do 
      end do 
!      stop                                                             
                                                                        
      return 
      END                                           
