      subroutine prolong(nxm,nym,nzm,cor,fin) 
!     ------------------------- --------------                          
!     call prolong(m,nx(m),ny(m),nz(m),p(loco(m)),p(loco(m-1)) )        
!     performs coarse to fine grid interpolation or "prolongation"      
!     m is the coarse grid level                                        
     use header, only : rc_kind

      integer i,j,k,iw,js,kb,nxm,nym,nzm 
      REAL(kind=rc_kind) :: fin(0:2*nxm+1,0:2*nym+1,0:2*nzm+1),              &
     &     cor(0:nxm+1,0:nym+1,0:nzm+1)                                 
      REAL(kind=rc_kind) :: an(8),bn(8),cn(8),dn(8),en(8),fn(8),             &
     &     gn(8),hn(8),x,y,z                                            
!                                                                       
!                                                                       
      x=0.25 
      y=0.25 
      z=0.25 
      call findn(x,y,z,an) 
!      write(6,*) an                                                    
      x= x+.5 
      call findn(x,y,z,bn) 
!      write(6,*) bn                                                    
      y= y+0.5 
      call findn(x,y,z,cn) 
!      write(6,*) cn                                                    
      x= x-0.5 
      call findn(x,y,z,dn) 
!      write(6,*) dn                                                    
      y=y-0.5 
      z= z+0.5 
      call findn(x,y,z,en) 
!      write(6,*) en                                                    
      x= x+.5 
      call findn(x,y,z,fn) 
!      write(6,*) fn                                                    
      y= y+0.5 
      call findn(x,y,z,gn) 
!      write(6,*) gn                                                    
      x= x-0.5 
      call findn(x,y,z,hn) 
!      write(6,*) hn                                                    
!                                                                       
!                                                                       
!      write(44,*) cor                                                  
      do k=0,nzm 
         kb= 2*k 
         do j=0,nym 
            js= 2*j 
            do i=0,nxm 
               iw = 2*i 
               fin(iw,js,kb)= fin(iw,js,kb) +                           &
     &              an(1)*cor(i,j,k) +an(2)*cor(i+1,j,k) +              &
     &              an(3)*cor(i+1,j+1,k) +an(4)*cor(i,j+1,k) +          &
     &              an(5)*cor(i,j,k+1) +an(6)*cor(i+1,j,k+1) +          &
     &              an(7)*cor(i+1,j+1,k+1) + an(8)*cor(i,j+1,k+1)       
               fin(iw+1,js,kb)= fin(iw+1,js,kb) +                       &
     &              bn(1)*cor(i,j,k) +bn(2)*cor(i+1,j,k) +              &
     &              bn(3)*cor(i+1,j+1,k) +bn(4)*cor(i,j+1,k) +          &
     &              bn(5)*cor(i,j,k+1) +bn(6)*cor(i+1,j,k+1) +          &
     &              bn(7)*cor(i+1,j+1,k+1) + bn(8)*cor(i,j+1,k+1)       
               fin(iw+1,js+1,kb)= fin(iw+1,js+1,kb) +                   &
     &              cn(1)*cor(i,j,k) +cn(2)*cor(i+1,j,k) +              &
     &              cn(3)*cor(i+1,j+1,k) +cn(4)*cor(i,j+1,k) +          &
     &              cn(5)*cor(i,j,k+1) +cn(6)*cor(i+1,j,k+1) +          &
     &              cn(7)*cor(i+1,j+1,k+1) + cn(8)*cor(i,j+1,k+1)       
               fin(iw,js+1,kb)= fin(iw,js+1,kb) +                       &
     &              dn(1)*cor(i,j,k) +dn(2)*cor(i+1,j,k) +              &
     &              dn(3)*cor(i+1,j+1,k) +dn(4)*cor(i,j+1,k) +          &
     &              dn(5)*cor(i,j,k+1) +dn(6)*cor(i+1,j,k+1) +          &
     &              dn(7)*cor(i+1,j+1,k+1) + dn(8)*cor(i,j+1,k+1)       
               fin(iw,js,kb+1)= fin(iw,js,kb+1) +                       &
     &              en(1)*cor(i,j,k) +en(2)*cor(i+1,j,k) +              &
     &              en(3)*cor(i+1,j+1,k) +en(4)*cor(i,j+1,k) +          &
     &              en(5)*cor(i,j,k+1) +en(6)*cor(i+1,j,k+1) +          &
     &              en(7)*cor(i+1,j+1,k+1) + en(8)*cor(i,j+1,k+1)       
               fin(iw+1,js,kb+1)= fin(iw+1,js,kb+1) +                   &
     &              fn(1)*cor(i,j,k) +fn(2)*cor(i+1,j,k) +              &
     &              fn(3)*cor(i+1,j+1,k) +fn(4)*cor(i,j+1,k) +          &
     &              fn(5)*cor(i,j,k+1) +fn(6)*cor(i+1,j,k+1) +          &
     &              fn(7)*cor(i+1,j+1,k+1) + fn(8)*cor(i,j+1,k+1)       
               fin(iw+1,js+1,kb+1)= fin(iw+1,js+1,kb+1) +               &
     &              gn(1)*cor(i,j,k) +gn(2)*cor(i+1,j,k) +              &
     &              gn(3)*cor(i+1,j+1,k) +gn(4)*cor(i,j+1,k) +          &
     &              gn(5)*cor(i,j,k+1) +gn(6)*cor(i+1,j,k+1) +          &
     &              gn(7)*cor(i+1,j+1,k+1) + gn(8)*cor(i,j+1,k+1)       
               fin(iw,js+1,kb+1)= fin(iw,js+1,kb+1) +                   &
     &              hn(1)*cor(i,j,k) +hn(2)*cor(i+1,j,k) +              &
     &              hn(3)*cor(i+1,j+1,k) +hn(4)*cor(i,j+1,k) +          &
     &              hn(5)*cor(i,j,k+1) +hn(6)*cor(i+1,j,k+1) +          &
     &              hn(7)*cor(i+1,j+1,k+1) + hn(8)*cor(i,j+1,k+1)       
            end do 
         end do 
      end do 
!      write(34,*) fin                                                  
      return 
      END                                           
                                                                        
                                                                        
      subroutine findn(x,y,z,pn) 
     use header, only : rc_kind

      REAL(kind=rc_kind) :: pn(8),x,y,z,sum 
      pn(1)= (1.d0-x)*(1.d0-y)*(1.d0-z) 
      pn(2)= x*(1.d0-y)*(1.d0-z) 
      pn(3)= x*y*(1.d0-z) 
      pn(4)= (1.d0-x)*y*(1.d0-z) 
      pn(5)= (1.d0-x)*(1.d0-y)*z 
      pn(6)= x*(1.d0-y)*z 
      pn(7)= x*y*z 
      pn(8)= (1.d0-x)*y*z 
      sum= pn(1)+pn(2)+pn(3)+pn(4)+pn(5)+pn(6)+pn(7)+pn(8) 
      if (dabs(sum-1.d0).gt.0.0001) write(6,*) 'coefficients wrong' 
      return 
      END                                           
