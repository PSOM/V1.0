subroutine pbc(cpf,fn,dtimel) 
!-------------------------------------------------------------------c   
  USE header
!     evaluates the boundaries in psolve c fills in                     
!     cp0,cpim1,cpip1,cpjm1,cpjp1,cpkm1,cpkp1,fn for the boundaries     
!     3-D Poisson solver - uses SOR                                     
!     Written for the control volume formulation                        
!     Uses the 19 point stencil (in 3d).                                
!                                                                       
!     fn is the contribution to the source terms on the RHS.            
!     Whenever possible we substitue the the velocities themselves      
!     (at the boundaries). We cannot substitute the u,w vels at the entr
!     Our eqn is                                                        
!     cpim1*p(i-1) +cpip1*p(i+1) +cpjm1p(j-1) ... -cp0*p(i,j,k)= f      
      integer i,j,k,l 
      REAL(kind=rc_kind) :: dtimel,edt,cpf(19,NI,NJ,NK),                      &
     &     fn(NI,NJ,NK)                                                 
!                                                                       
!     delinv= 1/delta, delta= D/L                                       
      edt= EPS/dtimel 
!                                                                       
!     *****************************                                     
!     do the edges separately                                           
      do 40 k=1,NK 
         do 50 j=1,NJ 
            do 60 i=1,NI 
!c        if ((i.eq.1).or.(i.eq.NI).or.(j.eq.1).or. changed for priodice
               if ((j.eq.1).or.                                         &
     &              (j.eq.NJ).or.(k.eq.1).or.(k.eq.NK)) then            
                  fn(i,j,k)= 0.d0 
                  do l=1,19 
                     cpf(l,i,j,k)= 0.d0 
                  end do 
!     add values for the east face if (i.ne.NI or i.eq.NI)              
                  fn(i,j,k)= fn(i,j,k) + edt*cxf(i,j,k)                 &
     &                 -sifc(i,j,k)                                     
                  cpf(1,i,j,k)= cpf(1,i,j,k) - gqi(i,j,k,1) 
                  cpf(2,i,j,k)= cpf(2,i,j,k) +gqi(i,j,k,1) 
                  cpf(3,i,j,k)= cpf(3,i,j,k) +0.25*gqi(i,j,k,2) 
                  cpf(5,i,j,k)= cpf(5,i,j,k) -0.25*gqi(i,j,k,2) 
                  cpf(6,i,j,k)= cpf(6,i,j,k) +0.25*gqi3(i,j,k) 
                  cpf(7,i,j,k)= cpf(7,i,j,k) -0.25*gqi3(i,j,k) 
                  cpf(10,i,j,k)= cpf(10,i,j,k) -0.25*gqi(i,j,k,2) 
                  cpf(11,i,j,k)= cpf(11,i,j,k) +0.25*gqi(i,j,k,2) 
                  cpf(13,i,j,k)= cpf(13,i,j,k) -0.25*gqi3(i,j,k) 
                  cpf(14,i,j,k)= cpf(14,i,j,k) +0.25*gqi3(i,j,k) 
!                                                                       
!     add values for the west face if (i.ne.1 or eq. 1)                 
                  fn(i,j,k)=  fn(i,j,k) -edt*cxf(i-1,j,k)               &
     &                 +sifc(i-1,j,k)                                   
                  cpf(1,i,j,k)= cpf(1,i,j,k) -gqi(i-1,j,k,1) 
                  cpf(3,i,j,k)= cpf(3,i,j,k) -0.25*gqi(i-1,j,k,2) 
                  cpf(4,i,j,k)= cpf(4,i,j,k) +gqi(i-1,j,k,1) 
                  cpf(5,i,j,k)= cpf(5,i,j,k) +0.25*gqi(i-1,j,k,2) 
                  cpf(6,i,j,k)= cpf(6,i,j,k) -0.25*gqi3(i-1,j,k) 
                  cpf(7,i,j,k)= cpf(7,i,j,k) +0.25*gqi3(i-1,j,k) 
                  cpf(8,i,j,k)= cpf(8,i,j,k) -0.25*gqi(i-1,j,k,2) 
                  cpf(9,i,j,k)= cpf(9,i,j,k) +0.25*gqi(i-1,j,k,2) 
                  cpf(12,i,j,k)= cpf(12,i,j,k) +0.25*gqi3(i-1,j,k) 
                  cpf(15,i,j,k)= cpf(15,i,j,k) -0.25*gqi3(i-1,j,k) 
                  if (j.eq.NJ) then 
                     fn(i,j,k)= fn(i,j,k) + edt*vfbcn(i,k) 
                  else 
!     add values for the north face if (j.ne.NJ)                        
                     fn(i,j,k)=  fn(i,j,k) +edt*cyf(i,j,k)              &
     &                    -sjfc(i,j,k)                                  
                     cpf(1,i,j,k)= cpf(1,i,j,k) -gqj(i,j,k,2) 
                     cpf(2,i,j,k)= cpf(2,i,j,k) +0.25*gqj(i,j,k,1) 
                     cpf(3,i,j,k)= cpf(3,i,j,k) +gqj(i,j,k,2) 
                     cpf(4,i,j,k)= cpf(4,i,j,k) -0.25*gqj(i,j,k,1) 
                     cpf(6,i,j,k)= cpf(6,i,j,k) +0.25*gqj3(i,j,k) 
                     cpf(7,i,j,k)= cpf(7,i,j,k) -0.25*gqj3(i,j,k) 
                     cpf(8,i,j,k)= cpf(8,i,j,k) -0.25*gqj(i,j,k,1) 
                     cpf(11,i,j,k)= cpf(11,i,j,k) +0.25*gqj(i,j,k,1) 
                     cpf(17,i,j,k)= cpf(17,i,j,k) -0.25*gqj3(i,j,k) 
                     cpf(18,i,j,k)= cpf(18,i,j,k) +0.25*gqj3(i,j,k) 
                  endif 
                  if (j.eq.1) then 
                     fn(i,j,k)= fn(i,j,k) - edt*vfbcs(i,k) 
                  else 
!     add values for the south face if (j.ne.1)                         
                     fn(i,j,k)=  fn(i,j,k) -edt*cyf(i,j-1,k)            &
     &                    +sjfc(i,j-1,k)                                
                     cpf(1,i,j,k)= cpf(1,i,j,k) -gqj(i,j-1,k,2) 
                     cpf(2,i,j,k)= cpf(2,i,j,k) -0.25*gqj(i,j-1,k,1) 
                     cpf(4,i,j,k)= cpf(4,i,j,k) +0.25*gqj(i,j-1,k,1) 
                     cpf(5,i,j,k)= cpf(5,i,j,k) +gqj(i,j-1,k,2) 
                     cpf(6,i,j,k)= cpf(6,i,j,k) -0.25*gqj3(i,j-1,k) 
                     cpf(7,i,j,k)= cpf(7,i,j,k) +0.25*gqj3(i,j-1,k) 
                     cpf(9,i,j,k)= cpf(9,i,j,k) +0.25*gqj(i,j-1,k,1) 
                     cpf(10,i,j,k)= cpf(10,i,j,k) -0.25*gqj(i,j-1,k,1) 
                     cpf(16,i,j,k)= cpf(16,i,j,k) +0.25*gqj3(i,j-1,k) 
                     cpf(19,i,j,k)= cpf(19,i,j,k) -0.25*gqj3(i,j-1,k) 
                  end if 
                  if  (k.eq.NK) then 
                     fn(i,j,k)= fn(i,j,k) + edt*czf(i,j,k) -skfc(i,j,k) 
                     cpf(1,i,j,k)= cpf(1,i,j,k) -2.d0*gqk(i,j,k,3) 
                  else 
!     add values for the top face if (k.ne.NK)                          
!c                     temp= fn(i,j,k)                                  
                     fn(i,j,k)=  fn(i,j,k) +edt*czf(i,j,k)              &
     &                    -skfc(i,j,k)                                  
                     cpf(1,i,j,k)= cpf(1,i,j,k) -gqk(i,j,k,3) 
                     cpf(2,i,j,k)= cpf(2,i,j,k) +0.25*gqk(i,j,k,1) 
                     cpf(3,i,j,k)= cpf(3,i,j,k) +0.25*gqk(i,j,k,2) 
                     cpf(4,i,j,k)= cpf(4,i,j,k) -0.25*gqk(i,j,k,1) 
                     cpf(5,i,j,k)= cpf(5,i,j,k) -0.25*gqk(i,j,k,2) 
                     cpf(6,i,j,k)= cpf(6,i,j,k) +gqk(i,j,k,3) 
                     cpf(14,i,j,k)= cpf(14,i,j,k) +0.25*gqk(i,j,k,1) 
                     cpf(15,i,j,k)= cpf(15,i,j,k) -0.25*gqk(i,j,k,1) 
                     cpf(18,i,j,k)= cpf(18,i,j,k) +0.25*gqk(i,j,k,2) 
                     cpf(19,i,j,k)= cpf(19,i,j,k) -0.25*gqk(i,j,k,2) 
                  endif 
                  if (k.eq.1) then 
                     fn(i,j,k)= fn(i,j,k) - edt*wfbcb(i,j) 
!c    since wf(k=0) is zero, we need not write this line                
                  else 
!                  if (k.ne.1) then                                     
!     add values for the bottom face if (k.ne.1)                        
!c                     temp= fn(i,j,k)                                  
                     fn(i,j,k)=  fn(i,j,k) -edt*czf(i,j,k-1)            &
     &                    +skfc(i,j,k-1)                                
!c                     if (k.eq.NK) write(400,*) temp,fn(i,j,k)         
                     cpf(1,i,j,k)= cpf(1,i,j,k) - gqk(i,j,k-1,3) 
                     cpf(2,i,j,k)= cpf(2,i,j,k) -0.25*gqk(i,j,k-1,1) 
                     cpf(3,i,j,k)= cpf(3,i,j,k) -0.25*gqk(i,j,k-1,2) 
                     cpf(4,i,j,k)= cpf(4,i,j,k) +0.25*gqk(i,j,k-1,1) 
                     cpf(5,i,j,k)= cpf(5,i,j,k) +0.25*gqk(i,j,k-1,2) 
                     cpf(7,i,j,k)= cpf(7,i,j,k) +gqk(i,j,k-1,3) 
                     cpf(12,i,j,k)= cpf(12,i,j,k) +0.25*gqk(i,j,k-1,1) 
                     cpf(13,i,j,k)= cpf(13,i,j,k) -0.25*gqk(i,j,k-1,1) 
                     cpf(16,i,j,k)= cpf(16,i,j,k) +0.25*gqk(i,j,k-1,2) 
                     cpf(17,i,j,k)= cpf(17,i,j,k) -0.25*gqk(i,j,k-1,2) 
                  endif 
               endif 
   60       continue 
   50    continue 
   40 continue 
!                                                                       
!c      do k=0,NK                                                       
!c         write(300,*) 'k=',k                                          
!c         do j=1,NJ                                                    
!c            write(300,*) 'j=',j                                       
!c            do i=1,NI                                                 
!c               write(300,*) skfc(i,j,k),wf(i,j,k)                     
!c            enddo                                                     
!c         enddo                                                        
!c      end do                                                          
      return 
      END                                           
