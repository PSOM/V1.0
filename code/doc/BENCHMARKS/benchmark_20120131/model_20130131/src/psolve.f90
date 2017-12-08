subroutine psolve(p,dtime,nstep,edt) 
  USE header, only : NI,NJ,NK,rc_kind
  !     --------------                                                    
  integer iter,i,j,k,nstep,l,imax,jmax,kmax 
  REAL(kind=rc_kind) :: cp(19,NI,NJ,NK),rhs(NI,NJ,NK),                   &
       &     p(0:NI+1,0:NJ+1,0:NK+1),res,                                 &
       &     dtime,tol,maxres,pstar,rlx,edt                               
  !                                                                       
  !     10^-6 is the tolerance on (u_x +v_y +eps*w_z) and the tolerance on
  !     residual r is  edt*(u_x +v_y +eps*w_z).                           
  !     c      tol= 1.d-6*edt                                             
  tol= 1.d-7*edt 
  !     compute fine grid coefficients and rhs on fine grid               
  call cpfine(dtime,cp,rhs) 
  !     write(6,*) 'cpfine called from psolve'                            
  !     c      call checkcp(1,NI,NJ,NK,cp)                                
  !     if (nstep.eq.3) call writout(rhs,NI,NJ,NK)                        
  !     call writout(cp,rhs,p,NI,NJ,NK)                                   
  !     **********************************                                
  !     write(6,*) 'enter rlx'                                            
  !     read(5,*) rlx                                                     
  !     rlx= 1.75                                                         
  rlx= 1.82 
  !     do 1000 iter=1,100                                                
  do 1000 iter=1,200 
     !     $DOACROSSSHARE(cp,p,rlx,rhs),                                     
     !     $&LOCAL(i,j,k,pstar)                                              
     do 110 k=1,NK 
        do 120 j=1,NJ 
           i=1 
           pstar= ( cp(2,i,j,k)*p(i+1,j,k)                          &
     &              +cp(3,i,j,k)*p(i,j+1,k)                             &
     &              +cp(4,i,j,k)*p(NI,j,k)                              &
     &              +cp(5,i,j,k)*p(i,j-1,k)                             &
     &              +cp(6,i,j,k)*p(i,j,k+1)                             &
     &              +cp(7,i,j,k)*p(i,j,k-1)                             &
     &              +cp(8,i,j,k)*p(NI,j+1,k)                            &
     &              +cp(9,i,j,k)*p(NI,j-1,k)                            &
     &              +cp(10,i,j,k)*p(i+1,j-1,k)                          &
     &              +cp(11,i,j,k)*p(i+1,j+1,k)                          &
     &              +cp(12,i,j,k)*p(NI,j,k-1)                           &
     &              +cp(13,i,j,k)*p(i+1,j,k-1)                          &
     &              +cp(14,i,j,k)*p(i+1,j,k+1)                          &
     &              +cp(15,i,j,k)*p(NI,j,k+1)                           &
     &              +cp(16,i,j,k)*p(i,j-1,k-1)                          &
     &              +cp(17,i,j,k)*p(i,j+1,k-1)                          &
     &              +cp(18,i,j,k)*p(i,j+1,k+1)                          &
     &              +cp(19,i,j,k)*p(i,j-1,k+1)                          &
     &              - rhs(i,j,k) )/(-cp(1,i,j,k))                       
               p(i,j,k)= (1.d0 -rlx)*p(i,j,k) +rlx*pstar 
               do 130 i=2,NI-1 
                  pstar= ( cp(2,i,j,k)*p(i+1,j,k)                       &
     &                 +cp(3,i,j,k)*p(i,j+1,k)                          &
     &                 +cp(4,i,j,k)*p(i-1,j,k)                          &
     &                 +cp(5,i,j,k)*p(i,j-1,k)                          &
     &                 +cp(6,i,j,k)*p(i,j,k+1)                          &
     &                 +cp(7,i,j,k)*p(i,j,k-1)                          &
     &                 +cp(8,i,j,k)*p(i-1,j+1,k)                        &
     &                 +cp(9,i,j,k)*p(i-1,j-1,k)                        &
     &                 +cp(10,i,j,k)*p(i+1,j-1,k)                       &
     &                 +cp(11,i,j,k)*p(i+1,j+1,k)                       &
     &                 +cp(12,i,j,k)*p(i-1,j,k-1)                       &
     &                 +cp(13,i,j,k)*p(i+1,j,k-1)                       &
     &                 +cp(14,i,j,k)*p(i+1,j,k+1)                       &
     &                 +cp(15,i,j,k)*p(i-1,j,k+1)                       &
     &                 +cp(16,i,j,k)*p(i,j-1,k-1)                       &
     &                 +cp(17,i,j,k)*p(i,j+1,k-1)                       &
     &                 +cp(18,i,j,k)*p(i,j+1,k+1)                       &
     &                 +cp(19,i,j,k)*p(i,j-1,k+1)                       &
     &                 - rhs(i,j,k) )/(-cp(1,i,j,k))                    
                  p(i,j,k)= (1.d0 -rlx)*p(i,j,k) +rlx*pstar 
  130          continue 
               i=NI 
               pstar= ( cp(2,i,j,k)*p(1,j,k)                            &
     &              +cp(3,i,j,k)*p(i,j+1,k)                             &
     &              +cp(4,i,j,k)*p(i-1,j,k)                             &
     &              +cp(5,i,j,k)*p(i,j-1,k)                             &
     &              +cp(6,i,j,k)*p(i,j,k+1)                             &
     &              +cp(7,i,j,k)*p(i,j,k-1)                             &
     &              +cp(8,i,j,k)*p(i-1,j+1,k)                           &
     &              +cp(9,i,j,k)*p(i-1,j-1,k)                           &
     &              +cp(10,i,j,k)*p(1,j-1,k)                            &
     &              +cp(11,i,j,k)*p(1,j+1,k)                            &
     &              +cp(12,i,j,k)*p(i-1,j,k-1)                          &
     &              +cp(13,i,j,k)*p(1,j,k-1)                            &
     &              +cp(14,i,j,k)*p(1,j,k+1)                            &
     &              +cp(15,i,j,k)*p(i-1,j,k+1)                          &
     &              +cp(16,i,j,k)*p(i,j-1,k-1)                          &
     &              +cp(17,i,j,k)*p(i,j+1,k-1)                          &
     &              +cp(18,i,j,k)*p(i,j+1,k+1)                          &
     &              +cp(19,i,j,k)*p(i,j-1,k+1)                          &
     &              - rhs(i,j,k) )/(-cp(1,i,j,k))                       
               p(i,j,k)= (1.d0 -rlx)*p(i,j,k) +rlx*pstar 
  120       continue 
  110    continue 
!                                                                       
         do 140 l=1,3 
            call  mgpfill(dtime,p) 
  140    continue 
         maxres= 0.d0 
!     $DOACROSSSHARE(cp,p,rhs,maxres,imax,jmax,kmax),                   
!     $&LOCAL(i,j,k,res)                                                
         do 210 k=1,NK 
            do 220 j=1,NJ 
!     res = b - A x'                                                    
               i=1 
               res= rhs(i,j,k)-                                         &
     &              (cp(1,i,j,k)*p(i,j,k)                               &
     &              +cp(2,i,j,k)*p(i+1,j,k)                             &
     &              +cp(3,i,j,k)*p(i,j+1,k)                             &
     &              +cp(4,i,j,k)*p(NI,j,k)                              &
     &              +cp(5,i,j,k)*p(i,j-1,k)                             &
     &              +cp(6,i,j,k)*p(i,j,k+1)                             &
     &              +cp(7,i,j,k)*p(i,j,k-1)                             &
     &              +cp(8,i,j,k)*p(NI,j+1,k)                            &
     &              +cp(9,i,j,k)*p(NI,j-1,k)                            &
     &              +cp(10,i,j,k)*p(i+1,j-1,k)                          &
     &              +cp(11,i,j,k)*p(i+1,j+1,k)                          &
     &              +cp(12,i,j,k)*p(NI,j,k-1)                           &
     &              +cp(13,i,j,k)*p(i+1,j,k-1)                          &
     &              +cp(14,i,j,k)*p(i+1,j,k+1)                          &
     &              +cp(15,i,j,k)*p(NI,j,k+1)                           &
     &              +cp(16,i,j,k)*p(i,j-1,k-1)                          &
     &              +cp(17,i,j,k)*p(i,j+1,k-1)                          &
     &              +cp(18,i,j,k)*p(i,j+1,k+1)                          &
     &              +cp(19,i,j,k)*p(i,j-1,k+1))                         
               if (abs(res).gt.maxres) then 
                  maxres= abs(res) 
                  imax= i 
                  jmax= j 
                  kmax= k 
               end if 
                                                                        
               do 230 i=2,NI-1 
!     res = b - A x'                                                    
                  res= rhs(i,j,k)-                                      &
     &                 (cp(1,i,j,k)*p(i,j,k)                            &
     &                 +cp(2,i,j,k)*p(i+1,j,k)                          &
     &                 +cp(3,i,j,k)*p(i,j+1,k)                          &
     &                 +cp(4,i,j,k)*p(i-1,j,k)                          &
     &                 +cp(5,i,j,k)*p(i,j-1,k)                          &
     &                 +cp(6,i,j,k)*p(i,j,k+1)                          &
     &                 +cp(7,i,j,k)*p(i,j,k-1)                          &
     &                 +cp(8,i,j,k)*p(i-1,j+1,k)                        &
     &                 +cp(9,i,j,k)*p(i-1,j-1,k)                        &
     &                 +cp(10,i,j,k)*p(i+1,j-1,k)                       &
     &                 +cp(11,i,j,k)*p(i+1,j+1,k)                       &
     &                 +cp(12,i,j,k)*p(i-1,j,k-1)                       &
     &                 +cp(13,i,j,k)*p(i+1,j,k-1)                       &
     &                 +cp(14,i,j,k)*p(i+1,j,k+1)                       &
     &                 +cp(15,i,j,k)*p(i-1,j,k+1)                       &
     &                 +cp(16,i,j,k)*p(i,j-1,k-1)                       &
     &                 +cp(17,i,j,k)*p(i,j+1,k-1)                       &
     &                 +cp(18,i,j,k)*p(i,j+1,k+1)                       &
     &                 +cp(19,i,j,k)*p(i,j-1,k+1))                      
                  if (abs(res).gt.maxres) then 
                     maxres= abs(res) 
                     imax= i 
                     jmax= j 
                     kmax= k 
                  end if 
  230          continue 
               i=NI 
!     res = b - A x'                                                    
               res= rhs(i,j,k)-                                         &
     &              (cp(1,i,j,k)*p(i,j,k)                               &
     &              +cp(2,i,j,k)*p(1,j,k)                               &
     &              +cp(3,i,j,k)*p(i,j+1,k)                             &
     &              +cp(4,i,j,k)*p(i-1,j,k)                             &
     &              +cp(5,i,j,k)*p(i,j-1,k)                             &
     &              +cp(6,i,j,k)*p(i,j,k+1)                             &
     &              +cp(7,i,j,k)*p(i,j,k-1)                             &
     &              +cp(8,i,j,k)*p(i-1,j+1,k)                           &
     &              +cp(9,i,j,k)*p(i-1,j-1,k)                           &
     &              +cp(10,i,j,k)*p(1,j-1,k)                            &
     &              +cp(11,i,j,k)*p(1,j+1,k)                            &
     &              +cp(12,i,j,k)*p(i-1,j,k-1)                          &
     &              +cp(13,i,j,k)*p(1,j,k-1)                            &
     &              +cp(14,i,j,k)*p(1,j,k+1)                            &
     &              +cp(15,i,j,k)*p(i-1,j,k+1)                          &
     &              +cp(16,i,j,k)*p(i,j-1,k-1)                          &
     &              +cp(17,i,j,k)*p(i,j+1,k-1)                          &
     &              +cp(18,i,j,k)*p(i,j+1,k+1)                          &
     &              +cp(19,i,j,k)*p(i,j-1,k+1))                         
               if (abs(res).gt.maxres) then 
                  maxres= abs(res) 
                  imax= i 
                  jmax= j 
                  kmax= k 
               end if 
  220       continue 
  210    continue 
!                                                                       
         if (maxres.gt.3000.d0) then 
            print*, 'psolve.f90: STOP. res too large, i,j,k,maxres(in psolve)=', &
     &           maxres,imax,jmax,kmax                                  
            stop 
         end if 
         print*, 'psolve.f90: ', iter,maxres,imax,jmax,kmax 
!     do l=1,3                                                          
!     call  mgpfill(dtime,p)                                            
!     end do                                                            
         if (maxres.lt.tol) goto 13 
 1000 continue 
   13 maxres= maxres/edt 
!     maxres/edt is the value of (u_x +v_y +ep*w_Z)                     
      print*, "psolve.f90", 'iter=',iter,maxres 
!     do l=1,3                                                          
!     call  mgpfill(dtime,p)                                            
!     end do                                                            
!                                                                       
!     c      call output(0,0)                                           
!     c      write(6,*) 'output called from psolve'                     
!      do 50 k=0,NK+1                                                   
!         do 50 j=0,NJ+1                                                
!            do 50 i=0,NI+1                                             
!               write(200,*) p(i,j,k)                                   
! 50   continue                                                         
!c      close(200)                                                      
!     stop                                                              
            return 
      END                                           
