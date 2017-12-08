 subroutine meanh(NI,NJ,h,hmean) 
!     -------------------------------                                   
!     computes the mean free-surface elevation within the domain        
 use header, only : rc_kind
 integer NI,NJ,i,j 
 REAL(kind=rc_kind) :: h(0:NI+1,0:NJ+1),hsum,hmean 
                                                                       
 hsum=0.d0 
 do i=1,NI 
   do j=1,NJ 
     hsum= hsum +h(i,j) 
   enddo
 enddo

 hmean= hsum/dble(NI*NJ) 

 if (abs(hmean).gt.100.d0) then 
   write(6,*) 'error, hmean=',hmean 
   stop 
 end if 

 return 
END                                           
