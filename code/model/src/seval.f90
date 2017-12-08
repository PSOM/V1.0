function seval(n, u, x, y, b, c, d) 
!-----------------------------------------------------------            
     use header, only : rc_kind
  integer n 
  REAL(kind=rc_kind) :: seval
  REAL(kind=rc_kind) ::  u, x(n), y(n), b(n), c(n), d(n) 
!                                                                       
!  this subroutine evaluates the cubic spline function                  
!                                                                       
!    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3 
!                                                                       
!    where  x(i) .lt. u .lt. x(i+1), using horner's rule                
!                                                                       
!  if  u .lt. x(1) then  i = 1  is used.                                
!  if  u .ge. x(n) then  i = n  is used.                                
!                                                                       
!  input..                                                              
!                                                                       
!    n = the number of data points                                      
!    u = the abscissa at which the spline is to be evaluated            
!    x,y = the arrays of data abscissas and ordinates                   
!    b,c,d = arrays of spline coefficients computed by spline           
!                                                                       
!  if  u  is not in the same interval as the previous call, then a      
!  binary search is performed to determine the proper interval.         
!                                                                       
  integer i, j, k 
  REAL(kind=rc_kind) :: dx 
  data i/1/ 
  if ( i .ge. n ) i = 1 
  if ( u .lt. x(i) ) go to 10 
  if ( u .le. x(i+1) ) go to 30 
!                                                                       
!  binary search                                                        
!                                                                       
10 i = 1 
  j = n+1 
20 k = (i+j)/2 
  if ( u .lt. x(k) ) j = k 
  if ( u .ge. x(k) ) i = k 
  if ( j .gt. i+1 ) go to 20 
  !                                                                       
!  evaluate spline                                                      
!                                                                       
30 dx = u - x(i) 
  seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i))) 
END function seval
