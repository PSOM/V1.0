      subroutine resid(m,nxm,nym,nzm,cp,p,fn,res,maxres) 
!     --------------------------------------------------                
!     call resid(m,nx(m),ny(m),nz(m),cp(loccp(m)),p(loco(m)),           
!    &       rhs(loci(m)),res )                                         
!                                                                       
!     check the residual                                                
     use header, only : rc_kind
      integer i,j,k,nxm,nym,nzm,m 
      REAL(kind=rc_kind) :: p(0:nxm+1,0:nym+1,0:nzm+1),                      &
     &     fn(nxm,nym,nzm),cp(19,nxm,nym,nzm),                          &
     &     res(nxm,nym,nzm),maxres                                      
!                                                                       
!                                                                       
      maxres= 0.d0 
!$DOACROSSSHARE(cp,p,fn,maxres,res),                                    
!$&   LOCAL(i,j,k)                                                      
      do 110 k=1,nzm 
         do 120 j=1,nym 
            do 130 i=1,nxm 
!     cpim1*p(i-1) +cpip1*p(i+1) +cpjm1p(j-1) ...+ fn = cp0*p(i,j,k)    
!     res = b - A x'                                                    
               res(i,j,k)= fn(i,j,k)-                                   &
     &              (cp(1,i,j,k)*p(i,j,k)                               &
     &              +cp(2,i,j,k)*p(i+1,j,k)                             &
     &              +cp(3,i,j,k)*p(i,j+1,k)                             &
     &              +cp(4,i,j,k)*p(i-1,j,k)                             &
     &              +cp(5,i,j,k)*p(i,j-1,k)                             &
     &              +cp(6,i,j,k)*p(i,j,k+1)                             &
     &              +cp(7,i,j,k)*p(i,j,k-1)                             &
     &              +cp(8,i,j,k)*p(i-1,j+1,k)                           &
     &              +cp(9,i,j,k)*p(i-1,j-1,k)                           &
     &              +cp(10,i,j,k)*p(i+1,j-1,k)                          &
     &              +cp(11,i,j,k)*p(i+1,j+1,k)                          &
     &              +cp(12,i,j,k)*p(i-1,j,k-1)                          &
     &              +cp(13,i,j,k)*p(i+1,j,k-1)                          &
     &              +cp(14,i,j,k)*p(i+1,j,k+1)                          &
     &              +cp(15,i,j,k)*p(i-1,j,k+1)                          &
     &              +cp(16,i,j,k)*p(i,j-1,k-1)                          &
     &              +cp(17,i,j,k)*p(i,j+1,k-1)                          &
     &              +cp(18,i,j,k)*p(i,j+1,k+1)                          &
     &              +cp(19,i,j,k)*p(i,j-1,k+1))                         
               if (abs(res(i,j,k)).gt.maxres) then 
                  maxres= abs(res(i,j,k)) 
               end if 
  130       continue 
  120    continue 
  110 continue 
!                                                                       
      if (maxres.gt.3000.d0) then 
         write(6,*) 'STOP. res too large, i,j,k,maxres=',               &
     &        maxres                                                    
         stop 
      end if 
!      write(6,*) 'level',m,'     maxres',maxres                        
  202 return 
      END                                           
