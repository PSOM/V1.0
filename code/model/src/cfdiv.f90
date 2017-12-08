subroutine cfdiv(maxdiv) 
  !     ------------------                                                
  USE header
  !     checks divergence of cxf,cyf,czf                                  
  !     czf does not contain  J*wt and is not exactly the contravariant ve
  integer i,j,k,imax,jmax,kmax 
  REAL(kind=rc_kind) :: cxfdx,cyfdy,czfdz,div,maxdiv 
!     edd is a scaling factor that makes this divergence                
!     the same as the residual from the pressure-Poisson equation.      
!     This is a check on the code.                                      
!                                                                       
   maxdiv= 0.d0 
   do k=1,NK 
     do j=1,NJ 
       do i=1,NI 
         cxfdx= (cxf(i,j,k)-cxf(i-1,j,k)) 
         cyfdy= (cyf(i,j,k)-cyf(i,j-1,k)) 
         czfdz= (czf(i,j,k)-czf(i,j,k-1)) 
                                                                  
         div= abs(cxfdx+ cyfdy + czfdz) 
         if (div.gt.maxdiv) then 
           maxdiv=div; imax=i; jmax=j; kmax=k; 
         end if 
      enddo
    enddo
  enddo
!      write(6,*) 'in facediv, i,j,k',imax,jmax,kmax                    
      return 
      END                                           
