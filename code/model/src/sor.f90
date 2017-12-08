      subroutine sor(nxm,nym,nzm,cp,p,fn) 
!     -----------------------------------------                         
!     call relax(nx(m),ny(m),nz(m),cp(loccp(m)),p(loco(m)),rhs(loci(m)))
!     solve  Grad p = fn                                                
     use header, only : rc_kind
      integer i,j,k,nxm,nym,nzm,iter,niter 
      parameter(niter=1) 
      REAL(kind=rc_kind) :: p(0:nxm+1,0:nym+1,0:nzm+1),                      &
     &     fn(nxm,nym,nzm),cp(19,nxm,nym,nzm),                          &
     &     rlx,pstar                                                    
                                                                        
!      rlx= 1.75                                                        
      rlx= 1.8 
      do 1000 iter=1,8 
!$DOACROSSSHARE(cp,p,rlx,fn),                                           
!$&LOCAL(i,j,k,pstar)                                                   
      do 110 k=1,nzm 
         do 120 j=1,nym 
            do 130 i=1,nxm 
               pstar= ( cp(2,i,j,k)*p(i+1,j,k)                          &
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
     &              +cp(19,i,j,k)*p(i,j-1,k+1)                          &
     &              - fn(i,j,k) )/(-cp(1,i,j,k))                        
               p(i,j,k)= (1.d0 -rlx)*p(i,j,k) +rlx*pstar 
  130       continue 
  120    continue 
  110 continue 
 1000 continue 
      return 
      END                                           
