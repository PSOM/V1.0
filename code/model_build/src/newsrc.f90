subroutine newsrc 
!----------------------------------------------------                   
  USE header
!     We interpolate the source terms onto the cell faces               
  integer i,j,k 
     
 
  do k=1,NK 
    do j=1,NJ 
      do i=0,NI 
        sifc(i,j,k) = 0.0 
      end do 
    end do 
  end do 
                                                                    
  do k=1,NK 
    do j=0,NJ 
      do i=1,NI 
        sjfc(i,j,k)= 0.0 
      end do 
    end do 
  end do 
                                                                    
  do k=0,NK 
    do j=1,NJ 
      do i=1,NI 
        skfc(i,j,k)= 0.0 
      end do 
    end do 
  end do 

  return 

  END                                           
