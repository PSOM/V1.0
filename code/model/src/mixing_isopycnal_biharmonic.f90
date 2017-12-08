subroutine mixing_isopycnal_biharmonic(var,vardif,biharmonic_coef) 
  !     ---------------------------------------------                     
  USE header
  !     use level n                                                       
  !     computes d2u/dx2 + d2u/dy2 at the cell centers.                   
  !                                                                       
  IMPLICIT NONE 

REAL :: biharmonic_coef
REAL :: kredi
REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1) :: var      ! var=(s,T,u,v,w,Tr(it,:,:,:)                    
REAL(kind=rc_kind), dimension(    1:NI  ,1:NJ  , 1:NK  ) :: vardif,vardif2   ! vardif is the source term from diabatic processes 
REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1) :: vardp
             

kredi=sqrt(biharmonic_coef)
PRINT*,kredi
vardif=0.;vardif2=0.;vardp=0.;
CALL mixing_isopycnal(var,vardif2,kredi)

vardp(1:NI,1:NJ,1:NK)=-vardif2(1:NI,1:NJ,1:NK)*Jacinv(1:NI,1:NJ,1:NK)
CALL mixing_isopycnal(vardp,vardif,kredi);

                                                          
 
  return 
END subroutine
