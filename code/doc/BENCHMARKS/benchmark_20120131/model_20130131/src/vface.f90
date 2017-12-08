  subroutine vface(pf,dtimel) 
!----------------------------------------------------                   
  USE header
!     modified for periodicew bc                                        
!     --------------------------                                        
!     compute the final vel (n+1 th step) at the cell faces             
!     applies the boundary conditions                                   
  integer i,j,k,nexit,iflag 
  REAL(kind=rc_kind) :: dte,px,py,pz,dtimel,kaph1,flowin,flowout,const,pf(0:NI+1,0:NJ+1,0:NK+1)                                     
!                                                                       
  dte= dtimel/EPS 
  kaph1= 1.d0 -kappah 
! 

! Be careful: pf is pcorr. 
                                                                      
!     uf compute here is used as the flux in convecn. In intpol it is   
!     used just at the outflow as the value of Uf at the prev time step.
!     In hbc and pbc, hnbc and pnbc, ufbce,ufbcw,vfbcn,vfbcs are used fo


!     x-direction                                                       
!     ----------- 
  do j=1,NJ 
    do i=1,NI 
      do k=1,NK 
        px=        (pf(i+1,j,k) -pf(i,j,k))                              * gqi(i,j,k,1) &
          & + 0.25*(pf(i+1,j+1,k)+pf(i,j+1,k)-pf(i+1,j-1,k)-pf(i,j-1,k)) * gqi(i,j,k,2) &
          & + 0.25*(pf(i+1,j,k+1)+pf(i,j,k+1)-pf(i+1,j,k-1)-pf(i,j,k-1)) * gqi3(i,j,k)                
        uf(i,j,k)= cxf(i,j,k) -dte*( px +sifc(i,j,k)) 
      enddo
    enddo
   enddo

  do k=1,NK 
     do j=1,NJ 
        uf(0,j,k)=uf(NI,j,k) 
     end do 
  end do 


!     y-direction                                                       
!     -----------                                                       
  do i=1,NI 
    do j=0,NJ 
      do k=1,NK 
        if (j.eq.0) then 
           vf(i,j,k)= vfbcs(i,k) 
        else if (j.eq.NJ) then 
           vf(i,j,k)= vfbcn(i,k) 
        else 
        py=        (pf(i,j+1,k)-pf(i,j,k))                               * gqj(i,j,k,2) &
          & + 0.25*(pf(i+1,j+1,k)+pf(i+1,j,k)-pf(i-1,j+1,k)-pf(i-1,j,k)) * gqj(i,j,k,1) &
          & + 0.25*(pf(i,j+1,k+1)+pf(i,j,k+1)-pf(i,j+1,k-1)-pf(i,j,k-1)) * gqj3(i,j,k)                               
        vf(i,j,k)= cyf(i,j,k) -dte*( py +sjfc(i,j,k)) 
        endif 
     enddo
   enddo
  enddo


!     z-direction                                                       
!     -----------                                                       
  do j=1,NJ 
    do i=1,NI 
      wf(i,j,0)= wfbcb(i,j) 
      do k=1,NK 
        pz=        (pf(i,j,k+1) -pf(i,j,k))                              * gqk(i,j,k,3)  &
          & + 0.25*(pf(i+1,j,k+1)+pf(i+1,j,k)-pf(i-1,j,k+1)-pf(i-1,j,k)) * gqk(i,j,k,1)  &
          & + 0.25*(pf(i,j+1,k+1)+pf(i,j+1,k)-pf(i,j-1,k+1)-pf(i,j-1,k)) * gqk(i,j,k,2)   
        wf(i,j,k)= czf(i,j,k) -dte*(pz +skfc(i,j,k)) 
     enddo
   enddo
  enddo

 return 
 END                                           
