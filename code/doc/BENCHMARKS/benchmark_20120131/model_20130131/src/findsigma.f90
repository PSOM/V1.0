subroutine findsigma(i,j,z,sigma,sigup)

USE header

 INTEGER :: i,j,sigma
 REAL :: z


 REAL :: zs, dsig,sigup, sigdn
 INTEGER :: k,l 

zs=z
    

if (zc(i,j,1)*DL .le. z) then 

  do k=2,NK 
    if (zc(i,j,k)*DL .ge. z) then 
      dsig= (zc(i,j,k)-zc(i,j,k-1) )*DL
      sigup= -(zc(i,j,k-1)*DL - z)/dsig 
      sigdn=  (zc(i,j,k)*DL   - z)/dsig 
      sigma=k
      ! zs(i,j)      = sigup*zc(i,j,k)     + sigdn*zc(i,j,k-1) 
      goto 102 
    end if 
  end do 

 else
  sigma=1;sigup=0.

end if 


102 l=0


end subroutine
