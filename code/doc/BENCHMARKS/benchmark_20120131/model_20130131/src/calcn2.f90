subroutine calcn2 
 !     ---------------------------                                       
 USE header
 !     Calculate N2  (units: per second)                                 
                                                                       
 implicit none 
 integer i,j,k 
 REAL(kind=rc_kind) :: rdy,rdz,rpdz,depth,dz 
 REAL(kind=rc_kind) :: avezf(NK),aveN2(NK) 
                                                                       
 !  evalrho was called from rpevalgrad                                
  DLinv= 1.d0/DL 
  do k=1,NK 
    do j=1,NJ 
      do i=1,NI 

        rdy= 0.5*(  ( rho(i,j+1,k) - SUM(rho(:,j+1,k))/NI ) - ( rho(i,j-1,k) - SUM(rho(:,j-1,k))/NI ) ) * vy(i,j)/LEN
       !  PRINT*,"P2 ",LEN,vy(i,j),rdy                             
         if (k.eq.NK) then 
           rdz= (rho(i,j,k) -rho(i,j,k-1))*wz(i,j,k)*DLinv 
           rpdz= (rho(i,j,k) - SUM(rho(:,j,k))/NI -( rho(i,j,k-1)- SUM(rho(:,j,k-1))/NI))*wz(i,j,k)*DLinv 
          else if (k.eq.1) then 
           rdz= (rho(i,j,k+1) -rho(i,j,k))*wz(i,j,k)*DLinv 
           rpdz= (rho(i,j,k+1) - SUM(rho(:,j,k+1))/NI -( rho(i,j,k)- SUM(rho(:,j,k))/NI))*wz(i,j,k)*DLinv 
          else 
           rdz= 0.5*(rho(i,j,k+1) -rho(i,j,k-1))*wz(i,j,k)*DLinv 
           rpdz= 0.5* (rho(i,j,k+1) - SUM(rho(:,j,k-1))/NI -( rho(i,j,k-1)- SUM(rho(:,j,k-1))/NI))*wz(i,j,k)*DLinv 
         end if 
                                                   
         freqby(i,j,k)=(-gpr*10.d0/R0)*rdy                 
         freqbz(i,j,k)=(-gpr*10.d0/R0)*rpdz                 
         freqN2(i,j,k)=(-gpr*10.d0/R0)*rdz 
                                                                   
      end do 
    end do 
  end do 
                                                                    
  return 
END                                           
