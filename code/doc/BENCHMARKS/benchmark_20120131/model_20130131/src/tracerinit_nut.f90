subroutine tracerinit(step) 
  !     --------------------                                              
  !     initializes tracer fields                                         
  !                                                                       
  USE header
  integer  i,j,k,l,it,step 
  REAL(kind=rc_kind) :: zt(ntr),zb(ntr),zmm(ntr),dsal,tauinv,vol
!       trinit(ntr,0:NJ+1,0:NK+1),   ! trinit is needed if tracersource used
  !      common/trsource/trinit 
!                                                                       
  if (step.gt.1) then 
     write(6,*) 'in tracerinit, trinit needs to be set correctly' 
     stop 
  end if
                                                                        
!     if this is a continuation of a run, only trinit is to be initializ
!     Tr 1 is like nutrient                                             
                                                                        
!    TRACER   beneath ML - correlated with density
!=============================================
  it=1
  rhomax= rho(NI/2,NJ/2,1) 
  rhomin= 1025.3
  i=NI/2
  do k=1,NK
     do j=1,NJ
        if (rho(i,j,k).gt.rhomin) then
           Trinit(j,k)= (1.-(rhomax- rho(i,j,k) )/(rhomax-rhomin))*40.
        else
           Trinit(j,k) = 0.d0
        end if
     end do
  end do
  do k=1,NK
     do j=1,NJ
        do i=1,NI
           Tr(it,i,j,k,0) = Trinit(j,k)
        end do
     end do
  end do

  !============================================
  !     Periodic Boundaries (periodic e-w)                                
  n=0
  do k=0,NK+1 
     do j=0,NJ+1 
        do it=1,ntr 
           Tr(it,0,j,k,n)= Tr(it,NI,j,k,n) 
           Tr(it,NI+1,j,k,n)= Tr(it,1,j,k,n) 
        end do
     end do
  end do
!

  return 
end subroutine tracerinit

  
