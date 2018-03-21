subroutine vhydro(dtimel) 
!     -----------------------------------------------                   
  USE header
!     modified for periodicew                                           
!     computes, u,v,w -tilde and writes these over cxf,cyf,czf.         
!     czf is unchanged.                                                 
  integer i,j,k 
  REAL(kind=rc_kind) :: hxi,heta,dte,dtimel,gradh,kaph1,temp,pz,          &
                  &     sumuf(0:NI,NJ),sumvf(NI,0:NJ),hx,hy,Jack,wfbct,ufdu,vfdv     
  dte= dtimel/EPS 
  kaph1= 1.d0 -kappah 
                                                                   
  do j=1,NJ 
    do i=1,NI 
      hxi= h(i+1,j)- h(i,j) 
      heta= 0.25*( h(i+1,j+1) + h(i,j+1) - h(i+1,j-1) -h (i,j-1))                                            
      do k=1,NK 
        hx= gi(i,j,k,1)*hxi +gi(i,j,k,2)*heta 
        gradh= gpr*(kappah*hx + kaph1*hxn(i,j,k)) 
        cxf(i,j,k)= cxf(i,j,k) -dte*(gradh +sifc(i,j,k)) 
        hxn(i,j,k)= hx 
      enddo
    enddo
                                                                   
    do k=1,NK 
      cxf(0,j,k)= cxf(NI,j,k) 
      hxn(0,j,k)= hxn(NI,j,k) 
    end do 
  end do 
                                                                   
  do i=1,NI 
    do j=0,NJ 
      hxi= 0.25*(h(i+1,j+1) +h(i+1,j) -h(i-1,j+1) -h(i-1,j)) 
      heta= h(i,j+1) -h(i,j) 
      do k=1,NK 
        hy= gj(i,j,k,1)*hxi +gj(i,j,k,2)*heta 
        gradh= gpr*(kappah*hy + kaph1*hyn(i,j,k)) 
        cyf(i,j,k)= cyf(i,j,k) -dte*(gradh +sjfc(i,j,k)) 
        hyn(i,j,k)= hy 
      enddo
    enddo
  enddo

  do i=1,NI 
    do k=1,NK 
      cyf(i,NJ,k)= vfbcn(i,k) 
    enddo   
    do k=1,NK 
      cyf(i,0,k)= vfbcs(i,k) 
    enddo   
  enddo   

  do j=1,NJ 
    do i=1,NI 
      do k=0,NK 
        pz=         (p(i,j,k+1) -p(i,j,k))                            * gqk(i,j,k,3) &
           & + 0.25*(p(i+1,j,k+1)+p(i+1,j,k)-p(i-1,j,k+1)-p(i-1,j,k)) * gqk(i,j,k,1) &
           & + 0.25*(p(i,j+1,k+1)+p(i,j+1,k)-p(i,j-1,k+1)-p(i,j-1,k)) * gqk(i,j,k,2) 
        czf(i,j,k)= czf(i,j,k) -dte*(pz +skfc(i,j,k)) 
     end do 
   end do 
  end do 
 
  call uvchy(dtimel) 
                                                                        
  return 
  END                                           
