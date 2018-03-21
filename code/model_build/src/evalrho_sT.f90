subroutine evalrho_sT(rhonew,nl2) 
!     ---------------------------------------------                     
  USE header,only:s,T,NI,NJ,NK,rc_kind
!  s is salinity, T is temp, rhonew is pot density
  implicit none
  REAL(kind=rc_kind) :: rhonew(0:NI+1,0:NJ+1,0:NK+1),potdens
  integer i,j,k,nl2 
  REAL(kind=rc_kind) :: se, Te
 
  do k=0,NK+1 
     do j=0,NJ+1 
        do i=0,NI+1 
           se=s(i,j,k,nl2)
           Te=T(i,j,k,nl2)
  
           rhonew(i,j,k)= potdens(s(i,j,k,nl2),T(i,j,k,nl2))
           !rhonew(i,j,k)= rho(i,j,k)
!            rhonew(i,j,k)=potdens(se,Te)
!           rhonew(i,j,k)=potdens(REAL(NK+1-k),20.)
        end do
     end do
  end do



  return 
END subroutine evalrho_sT
