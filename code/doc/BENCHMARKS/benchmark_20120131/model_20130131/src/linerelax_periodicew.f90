subroutine linerelax(nxm,nym,nzm,cp,p,fn) 
!     -----------------------------------------                         
  USE header, only : NI,NJ,NK,rc_kind
!     This routine is called in place of relax from the multigrid solver
!     "mgrid".  It does not work (don't  know why).                     
!     call relax(nx(m),ny(m),nz(m),cp(loccp(m)),p(loco(m)),rhs(loci(m)))
!     solve  Grad p = fn                                                
      integer i,j,k,nxm,nym,nzm,info,im1,ip1,ii,istart,iseed 
      REAL(kind=rc_kind) :: p(0:nxm+1,0:nym+1,0:nzm+1),                      &
     &     fn(nxm,nym,nzm),cp(19,nxm,nym,nzm),subd(NK),supd(NK),        &
     &     dia(NK),rh(NK),dum,ran3                                      
!                                                                       
!==                                                                     
      iseed= 91294 
      dum = ran3(iseed) 
                                                                        
!      do k=0,nzm+1                                                     
!         do j=0,nym+1                                                  
!            do i=0,nxm +1                                              
!               p(i,j,k)= 0.0                                           
!            end do                                                     
!         end do                                                        
!      end do                                                           
!==                                                                     
!$DOACROSSSHARE(cp,p,fn,nxm,nym,nzm),                                   
!$&LOCAL(i,j,k,rh,subd,supd,dia,info)                                   
      do 110 j=1,nym 
         istart = int(ran3(iseed)*(nxm+1)) 
         do 120 ii=1,nxm 
!            write(6,*) 'istart = ',istart                              
            i = ii + istart 
            if (i.gt.nxm) i = i -nxm 
            im1= i-1 
            ip1= i+1 
!c periodic-ew boundaries                                               
            if (i.eq.1) im1= nxm 
            if (i.eq.nxm) ip1= 1 
            do 130 k=1,nzm 
               rh(k)= -( cp(2,i,j,k)*p(ip1,j,k)                         &
     &              +cp(3,i,j,k)*p(i,j+1,k)                             &
     &              +cp(4,i,j,k)*p(im1,j,k)                             &
     &              +cp(5,i,j,k)*p(i,j-1,k)                             &
     &              +cp(8,i,j,k)*p(im1,j+1,k)                           &
     &              +cp(9,i,j,k)*p(im1,j-1,k)                           &
     &              +cp(10,i,j,k)*p(ip1,j-1,k)                          &
     &              +cp(11,i,j,k)*p(ip1,j+1,k)                          &
     &              +cp(12,i,j,k)*p(im1,j,k-1)                          &
     &              +cp(13,i,j,k)*p(ip1,j,k-1)                          &
     &              +cp(14,i,j,k)*p(ip1,j,k+1)                          &
     &              +cp(15,i,j,k)*p(im1,j,k+1)                          &
     &              +cp(16,i,j,k)*p(i,j-1,k-1)                          &
     &              +cp(17,i,j,k)*p(i,j+1,k-1)                          &
     &              +cp(18,i,j,k)*p(i,j+1,k+1)                          &
     &              +cp(19,i,j,k)*p(i,j-1,k+1)                          &
     &              - fn(i,j,k) )                                       
!     &              +cp(6,i,j,k)*p(i,j,k+1)                            
!     &              +cp(7,i,j,k)*p(i,j,k-1)                            
!     &              - fn(i,j,k) )/(-cp(1,i,j,k))                       
               if (k.ne.1) subd(k)= cp(7,i,j,k) 
               if (k.ne.nzm) supd(k)= cp(6,i,j,k) 
               dia(k)= cp(1,i,j,k) 
  130       continue 
            rh(1)= rh(1) -cp(7,i,j,1)*p(i,j,0) 
            rh(nzm)= rh(nzm) -cp(6,i,j,nzm)*p(i,j,nzm+1) 
!                                                                       
!==                                                                     
!            if (j.eq.40) write(6,*) i, rh(12),fn(i,j,12),              
!     &           cp(1,i,j,12)                                          
!==                                                                     
            call dgtsl(nzm,subd,dia,supd,rh,info) 
            if (info.ne.0) then 
               write(6,*) 'error in linerelax' 
               stop 
            endif 
            do 140 k=1,nzm 
               p(i,j,k)= rh(k) 
  140       continue 
  120    continue 
  110 continue 
!                                                                       
      return 
!==                                                                     
      stop 
      k=12 
      j=40 
      write(6,*) 'linerelax' 
      do i=1,NI 
         write(6,*) p(i,j,k),fn(i,j,k),(cp(im1,i,j,k),im1=1,19) 
      end do 
      write(6,*) 'stopping in linerelax' 
      stop 
      return 
      END                                           
