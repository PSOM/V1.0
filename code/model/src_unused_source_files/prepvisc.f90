subroutine prepvisc(n,duml,hnusc,gxi,geta) 
  !---------------------------------------------------------------        
  USE header
  !     Called from Viscous to Evaluate the viscous terms                 
  integer i,j,k,n 
      REAL(kind=rc_kind) :: vnusc,hnusc,                                     &
     &     gk3nu(NI,NJ,0:NK), duml(0:NI+1,0:NJ+1,0:NK+1,0:1),            &
     &     gxi(0:NI,NJ,NK),geta(NI,0:NJ,NK),gsig(NI,NJ,0:NK)            
!     g33 on the k face has to been re-evaluated                        
!                                                                       
      do 50 k=1,NK 
         do 60 j=1,NJ 
            do 70 i=0,NI 
               gxi(i,j,k)= (gi(i,j,k,1)*(duml(i+1,j,k,n)-duml(i,j,k,n))   &
     &              +0.25*( gi(i,j,k,2)*(duml(i+1,j+1,k,n)+duml(i,j+1,k,n)&
     &              -duml(i+1,j-1,k,n)-duml(i,j-1,k,n)) + gi3(i,j,k)*(    &
     &              duml(i+1,j,k+1,n)+duml(i,j,k+1,n) -duml(i+1,j,k-1,n)   &
     &              -duml(i,j,k-1,n)) ) )*hnusc                          
   70       continue 
   60    continue 
   50 continue 
!                                                                       
      do 150 k=1,NK 
         do 160 j=0,NJ 
            do 170 i=1,NI 
               geta(i,j,k)= (gj(i,j,k,2)*(duml(i,j+1,k,n)-duml(i,j,k,n))  &
     &              +0.25*( gj(i,j,k,1)*(duml(i+1,j+1,k,n)+duml(i+1,j,k,n)&
     &              -duml(i-1,j+1,k,n)-duml(i-1,j,k,n)) + gj3(i,j,k)*(    &
     &              duml(i,j+1,k+1,n)+duml(i,j,k+1,n) -duml(i,j+1,k-1,n)   &
     &              -duml(i,j,k-1,n)) ) )*hnusc                          
  170       continue 
  160     continue 
  150  continue 
!                                                                       
!      do 250 k=0,NK                                                    
!         do 260 j=1,NJ                                                 
!            do 270 i=1,NI                                              
!               gsig(i,j,k)= gk3nu(i,j,k)*(duml(i,j,k+1,n)-duml(i,j,k,n)) 
!     &              +0.25*hnusc*( gk(i,j,k,1)*(duml(i+1,j,k+1,n)        
!     &              +duml(i+1,j,k,n)-duml(i-1,j,k+1,n)-duml(i-1,j,k,n))   
!     &              + gk(i,j,k,2)*(duml(i,j+1,k+1,n)+duml(i,j+1,k,n)     
!     &              -duml(i,j-1,k+1,n) -duml(i,j-1,k,n)) )               
! 270        continue                                                   
! 260      continue                                                     
! 250   continue                                                        
!                                                                       
       return 
      END                                           
