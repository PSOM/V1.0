 subroutine facediv(dtimel,maxdiv) 
!     ------------------                                                
!     checks divergence                                                 
!     wf does not contain  J*wt and is not exactly the contravariant vel
 USE header
 integer i,j,k,imax,jmax,kmax 
 REAL(kind=rc_kind) :: ufdx,vfdy,wfdz,div,maxdiv,dtimel
!     edd is a scaling factor that makes this divergence                
!     the same as the residual from the pressure-Poisson equation.      
!     This is a check on the code.                                      
!                                                                       
  maxdiv= 0.d0 
  do k=1,NK 
    do j=1,NJ 
      do i=1,NI 
        ufdx= (uf(i,j,k)-uf(i-1,j,k)) 
        vfdy= (vf(i,j,k)-vf(i,j-1,k)) 
        wfdz= (wf(i,j,k)-wf(i,j,k-1)) 
        div= abs(ufdx+ vfdy + wfdz) 
        if (div.gt.maxdiv) then 
              maxdiv=div; imax=i; jmax=j; kmax=k; 
        end if
     enddo
   enddo
 enddo 
 return 
 END                                           
