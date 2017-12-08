subroutine hsolve(h,oldh,hdt,dtime) 
!     --------------                                                    
  USE header, only : NI,NJ,rc_kind
!     modified for periodicew boundaries                                
!     ----------------------------------                                
      integer iter,i,j,kount,l 
      REAL(kind=rc_kind) :: ch(9,NI,NJ),rhs(NI,NJ),                          &
     &     oldh(0:NI+1,0:NJ+1),h(0:NI+1,0:NJ+1),res,                    &
     &     hdt(0:NI+1,0:NJ+1),dti,dtime,tol,maxres,hstar,rlx            
!     &     dti,dtime,tol,maxres,hstar,rlx                              
!                                                                       
!c      tol= 1.d-8                                                      
!     tol of 1.d-7 doesn't seem to change the solution. 1.d-6 would be o
      tol= 1.d-12 
!      tol= 1.d-18                                                      
      dti= 1.d0/dtime 
      kount=0 
!                                                                       
!     maxres= 0.d0                                                      
!     res= 0.d0                                                         
!     iter=0                                                            
!     write(6,*) 'in hsolve',tol,dti,maxres,res,iter                    
!                                                                       
!     save old value of h                                               
!                                                                       
!     compute fine grid coefficients and rhs on fine grid               
      call chfine(dtime,ch,rhs) 
!      call checkch(1,NI,NJ,ch)                                         
!     **********************************                                
!      write(6,*) 'enter rlx'                                           
!      read(5,*) rlx                                                    
!      do j=1,NJ                                                        
!         write(100,*) 'j=',j                                           
!         do i=1,NI                                                     
!            write(100,*) (ch(l,i,j),l=1,9),rhs(i,j)                    
!         end do                                                        
!      end do                                                           
!      stop                                                             
                                                                        
      rlx= 1.72 
  999 do 1000 iter=1,100 
!$DOACROSSSHARE(ch,h,rlx,rhs),                                          
!$&LOCAL(i,j,hstar)                                                     
         do 120 j=1,NJ 
            i=1 
            hstar= ( ch(2,i,j)*h(i+1,j)                                 &
     &           +ch(3,i,j)*h(NI,j)                                     &
     &           +ch(4,i,j)*h(i,j+1)                                    &
     &           +ch(5,i,j)*h(i,j-1)                                    &
     &           +ch(6,i,j)*h(i+1,j+1)                                  &
     &           +ch(7,i,j)*h(NI,j+1)                                   &
     &           +ch(8,i,j)*h(i+1,j-1)                                  &
     &           +ch(9,i,j)*h(NI,j-1)                                   &
     &           - rhs(i,j) )/(-ch(1,i,j))                              
            h(i,j)= (1.d0 -rlx)*h(i,j) +rlx*hstar 
            do 130 i=2,NI-1 
               hstar= ( ch(2,i,j)*h(i+1,j)                              &
     &              +ch(3,i,j)*h(i-1,j)                                 &
     &              +ch(4,i,j)*h(i,j+1)                                 &
     &              +ch(5,i,j)*h(i,j-1)                                 &
     &              +ch(6,i,j)*h(i+1,j+1)                               &
     &              +ch(7,i,j)*h(i-1,j+1)                               &
     &              +ch(8,i,j)*h(i+1,j-1)                               &
     &              +ch(9,i,j)*h(i-1,j-1)                               &
     &              - rhs(i,j) )/(-ch(1,i,j))                           
               h(i,j)= (1.d0 -rlx)*h(i,j) +rlx*hstar 
  130       continue 
            i=NI 
            hstar= ( ch(2,i,j)*h(1,j)                                   &
     &           +ch(3,i,j)*h(i-1,j)                                    &
     &           +ch(4,i,j)*h(i,j+1)                                    &
     &           +ch(5,i,j)*h(i,j-1)                                    &
     &           +ch(6,i,j)*h(1,j+1)                                    &
     &           +ch(7,i,j)*h(i-1,j+1)                                  &
     &           +ch(8,i,j)*h(1,j-1)                                    &
     &           +ch(9,i,j)*h(i-1,j-1)                                  &
     &           - rhs(i,j) )/(-ch(1,i,j))                              
            h(i,j)= (1.d0 -rlx)*h(i,j) +rlx*hstar 
  120    continue 
!                                                                       
         do 140 l=1,3 
            call hfill(dtime,h) 
  140    continue 
                                                                        
         maxres= 0.d0 
!$DOACROSSSHARE(ch,h,rhs,maxres),                                       
!$&LOCAL(i,j,res)                                                       
         do 220 j=1,NJ 
            i=1 
            res= rhs(i,j)-                                              &
     &           (ch(1,i,j)*h(i,j)                                      &
     &           +ch(2,i,j)*h(i+1,j)                                    &
     &           +ch(3,i,j)*h(NI,j)                                     &
     &           +ch(4,i,j)*h(i,j+1)                                    &
     &           +ch(5,i,j)*h(i,j-1)                                    &
     &           +ch(6,i,j)*h(i+1,j+1)                                  &
     &           +ch(7,i,j)*h(NI,j+1)                                   &
     &           +ch(8,i,j)*h(i+1,j-1)                                  &
     &           +ch(9,i,j)*h(NI,j-1))                                  
            if (abs(res).gt.maxres) then 
               maxres= abs(res) 
            end if 
            do 230 i=2,NI-1 
!     res = b - A x'                                                    
               res= rhs(i,j)-                                           &
     &              (ch(1,i,j)*h(i,j)                                   &
     &              +ch(2,i,j)*h(i+1,j)                                 &
     &              +ch(3,i,j)*h(i-1,j)                                 &
     &              +ch(4,i,j)*h(i,j+1)                                 &
     &              +ch(5,i,j)*h(i,j-1)                                 &
     &              +ch(6,i,j)*h(i+1,j+1)                               &
     &              +ch(7,i,j)*h(i-1,j+1)                               &
     &              +ch(8,i,j)*h(i+1,j-1)                               &
     &              +ch(9,i,j)*h(i-1,j-1))                              
               if (abs(res).gt.maxres) then 
                  maxres= abs(res) 
               end if 
  230       continue 
            i=NI 
            res= rhs(i,j)-                                              &
     &           (ch(1,i,j)*h(i,j)                                      &
     &           +ch(2,i,j)*h(1,j)                                      &
     &           +ch(3,i,j)*h(i-1,j)                                    &
     &           +ch(4,i,j)*h(i,j+1)                                    &
     &           +ch(5,i,j)*h(i,j-1)                                    &
     &           +ch(6,i,j)*h(1,j+1)                                    &
     &           +ch(7,i,j)*h(i-1,j+1)                                  &
     &           +ch(8,i,j)*h(1,j-1)                                    &
     &           +ch(9,i,j)*h(i-1,j-1))                                 
            if (abs(res).gt.maxres) then 
               maxres= abs(res) 
            end if 
  220    continue 
!                                                                       
!         write(6,*) 'iter,maxres',iter,maxres                          
         if (maxres.gt.3000.d0) then 
            print*, 'hsolve.f90: STOP. res too large, i,j,maxres=',              &
     &           i,j,maxres                                             
            stop 
         end if 
!     may have to call hfill a few times if g12 is non-zero             
!     if g12 is zero, we may call hfill after convergence is            
!     reached or may not call it at all.                                
!+         do 240 l=1,3                                                 
!+            call  hfill(dtime,h)                                      
!+ 240     continue                                                     
!         call mprove(h,ch,rhs,dtime)                                   
         if (maxres.lt.tol) goto 13 
 1000 continue 
!   13 print*,  'hsolve.f90: iter=',iter,maxres 
 13 continue
!      call  hfill(dtime,h)                                             
!+      kount=kount+1                                                   
!                                                                       
!     Bring the Coriolis terms up to the (n+1) time level. New sifc, sjf
!+      if (kount.ge.3) goto 301                                        
!+      call uvinterm(h,dtime)                                          
!+      call coriolis(1)                                                
!+      call srcface(1)                                                 
!+      call newrhs(dtime,rhs)                                          
!+      goto 999                                                        
!                                                                       
!c     hdt is computed in vhydro for i=1,NI,j=1,NJ... but for outer poin
!c     hdt is actually (dh'/dt')*(H/D)                                  
! 301  do 310 j=0,NJ+1                                                  
!         hdt(0,j)= dti*( h(0,j) -oldh(0,j) )*HDL                       
!         hdt(NI+1,j)= dti*( h(NI+1,j) -oldh(NI+1,j) )*HDL              
! 310  continue                                                         
!      do 320 i=1,NI                                                    
!         hdt(i,0)= dti*( h(i,0) -oldh(i,0))*HDL                        
!         hdt(i,NJ+1)= dti*( h(i,NJ+1) -oldh(i,NJ+1))*HDL               
! 320  continue                                                         
!                                                                       
      do l=1,1 
         call mprove(h,ch,rhs,dtime) 
      end do 
                                                                        
      return 
                                                                        
!!$!      do j=16,18                                                       
!!$         write(*,*) 'j =',j 
!!$         write(*,*) 'j =',j 
!!$         write(*,*) 'j =',j 
!!$         do i=1,NI 
!!$            write(*,*) 'rhs',rhs(i,17),(ch(l,i,17),l=1,2) 
!!$            write(*,*) 'ch',(ch(l,i,17),l=3,5) 
!!$         end do 
!!$         do i=1,NI 
!!$            write(*,*) h(i,16),h(i,17),h(i,18) 
!!$         end do 
!!$!      end do                                                           
      stop 
                                                                        
  301 return 
      END                                           
