subroutine eddydivergence(dtime)
  !     --------------                                                    
  USE header
  !     Calculate divergence of mean and eddy fluxes
                                                                        
  implicit NONE 
  integer it,i,j,k
  REAL(kind=rc_kind) :: vfbar(0:NJ,NK),wfbar(NJ,0:NK),vfcpbar(0:NJ,NK),  &
       wfcpbar(NJ,0:NK),cbarold(0:NJ,0:NK),vfp,wfp,cp,dtime,fac
! vfcpbar is bar(vf' c') and wfcpbar is bar(wf' c')

! First save the old value of cbar to cbarold (cbar init in tracerinit)
! -------------------------------------------
  cbarold=cbar

  it=1
  !vfbar      
  do k=1,NK
     do j=0,NJ
        vfbar(j,k)= 0.
        do i=1,NI
           vfbar(j,k) = vf(i,j,k)+ vfbar(j,k)
        end do
        vfbar(j,k)= vfbar(j,k)/dble(NI)
     end do
  end do
  !wfbar
  do k=0,NK
     do j=1,NJ
        wfbar(j,k)= 0.
        do i=1,NI
           wfbar(j,k) = wf(i,j,k)+ wfbar(j,k)
        end do
        wfbar(j,k)= wfbar(j,k)/dble(NI)
     end do
  end do
  !cbar
  do k=0,NK
     do j=0,NJ
        !cbar at k=0, j=0 not used, because wfbar(k=0)=0 and vfbar(j=0)=0
        cbar(j,k)= 0.
        do i=1,NI
           cbar(j,k)= Tr(it,i,j,k,0) + cbar(j,k)
        end do
        cbar(j,k) = cbar(j,k)/dble(NI)
     end do
  end do
  
  !Divergence of Mean flux      
  !---------------------
  do k=1,NK
     do j=1,NJ
        divmean(j,k)= vfbar(j,k)*cbar(j,k)- vfbar(j-1,k)*cbar(j-1,k)  &
             + wfbar(j,k)*cbar(j,k) - wfbar(j,k-1)*cbar(j,k-1)
     end do
  end do

  !Divergence of Reynolds averaged flux      
  !--------------------------------------
  !calculate vfcpbar
  do k=1,NK
     do j=0,NJ
        vfcpbar(j,k)=0.
        do i=1,NI
           vfp= vf(i,j,k)-vfbar(j,k)
           cp= Tr(it,i,j,k,0)-cbar(j,k)
           vfcpbar(j,k)= vfcpbar(j,k)+ vfp*cp
        end do
        vfcpbar(j,k)= vfcpbar(j,k)/dble(NI)
     end do
  end do

  !calcualte wfcpbar
  do k=0,NK
     do j=1,NJ
        wfcpbar(j,k)=0.
        do i=1,NI
           wfp= wf(i,j,k)-wfbar(j,k)
           cp= Tr(it,i,j,k,0)-cbar(j,k)
           wfcpbar(j,k)= wfcpbar(j,k)+ wfp*cp
        end do
        wfcpbar(j,k)= wfcpbar(j,k)/dble(NI)
     end do
  end do

  !calculate divergence
  do k=1,NK
     do j=1,NJ
        divreyn(j,k)= vfcpbar(j,k)- vfcpbar(j-1,k)                       &
             + wfcpbar(j,k) - wfcpbar(j,k-1)
     end do
  end do

  !To convert to Milimoles per second  (from non-dim concen per non-dim time)
  !Multiply by (LEN*LEN*DL / TL)  = LEN*DL*UL 
  !TO outut in MOLES / Second, divide by 10^3
  fac= LEN*DL*UL*1.d-3
  do k=1,NK
     do j=1,NK
        divreyn(j,k)= divreyn(j,k)*fac
        divmean(j,k)= divmean(j,k)*fac

     end do
  end do

! CALCULATE dcdt =  vol *(dc/dt)/ TL,  where vol= Jac*LEN*LEN*DL
! dcdt = (dc/dt ) *jac*LEN*DL*UL
! fac=LEN*DL*UL*1.d-3   1.d-3 converts to mols.  So dcdt is in MOLES / SECOND

  i=NI/2
  do k=1,NK
     do j=1,NJ
        dcdt(j,k)= Jac(i,j,k)*fac*(cbar(j,k)- cbarold(j,k))/dtime
     end do
  end do

  return
end subroutine eddydivergence
