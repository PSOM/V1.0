subroutine pfilter(p) 
  !----------------------------------------------------                   
  USE header, only : NI,NJ,NK,rc_kind
  !     filter the correction to the NH presure                           
  !     --------------------------                                        
  !                                                                       
      integer i,j,k,m 
      REAL(kind=rc_kind) :: p(0:NI+1,0:NJ+1,0:NK+1) 
                                                                        
      do m=1,3 
      do k=0,NK+1 
         do j=1,NJ 
            do i=0,NI+1 
               p(i,j,k)= 0.5*p(i,j,k)+0.25*(p(i,j-1,k)                  &
     &              +p(i,j+1,k))                                        
            end do 
         end do 
      end do 
      end do 
      return 
      END                                           
