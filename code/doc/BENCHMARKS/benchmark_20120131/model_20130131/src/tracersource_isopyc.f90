subroutine tracersource(n,dtime) 
  !     -----------------------                                           
  USE header
                                                                        
!     tracer sources                                                    
!                                                                       
  integer  i,j,k,n,it,m 
  REAL(kind=rc_kind) :: dtime,fac(ntr) 
  REAL(kind=rc_kind) :: trinit,facmortal,phyconvert 
!                                                                       
!     consump is the uptake (or if negative, then addition) or tracer   
!     in milimoles/m^3/one time step                                    
                                                                        
!     tau = 1,32 days                                                   
!     nondimtau= 3*86400*UL/LEN                                         
  fac(1)= dtime/(1.d0*86400.d0*UL/LEN) 
  fac(2)= dtime/(2.d0*86400.d0*UL/LEN) 
  fac(3)= dtime/(4.d0*86400.d0*UL/LEN) 
  fac(4)= dtime/(8.d0*86400.d0*UL/LEN) 
  fac(5)= dtime/(16.d0*86400.d0*UL/LEN) 
  fac(6)= dtime/(32.d0*86400.d0*UL/LEN) 
                                                                    
!     mortality time scale greater by a factor of 10                    
  facmortal = dtime/(30.d0*86400.d0*UL/LEN) 
!=      fac(3)= dtime/(12.d0*86400.d0*UL/LEN)                           
                                                                        
   do k=0,NK+1 
     do j=0,NJ+1 
       do i=0,NI+1 
         if (rho(i,j,k).lt.1025.2) then 
           trinit= 1.d0 
          else 
           trinit= 0.d0 
         end if 
                                                                     
         !     RESTORE ALL TRACERS                                               
         do it=1,ntr 
           consump(i,j,k,it)= fac(it)*(T(it,i,j,k,n)- trinit)                          
           T(it,i,j,k,n)= T(it,i,j,k,n) - consump(i,j,k,it) 
         end do 
       end do 
     end do 
   end do 
                                                                     
   do m=0,1 
     do k=0,NK+1 
       do j=0,NJ+1 
         do it=1,ntr 
           t(it,0,j,k,m)= t(it,NI,j,k,m) 
           t(it,NI+1,j,k,m)= t(it,1,j,k,m) 
          end do 
       end do 
     end do 
   end do 
                                                                        


      return 




!     salnity forcing                                                   
!     ---------------------------                                       
      alpha = 0.7*864000.0*UL/LEN 
      dtfac = dtime/alpha 
      sfresh  = 34.0 -S0 
      vol = Jac(NI/2,1,NK) 
                                                                        
                                                                        
      do j=1,NJ 
         if (yc(j).gt.yfront) then 
            shelfdep= D(1,j-1) 
            goto 101 
         end if 
      end do 
                                                                        
                                                                        
  101 do j=0,NJ+1 
!         if (yc(j).le.yfront) then                                     
!            facy =  (yfront - yc(j))/yfront                            
                                                                        
         do  i=0,NI+1 
            do k=0,NK 
!               fac = (1.0 - zf(i,j,k)/shelfdep ) *facy                 
!=               fac = (1.0 - zf(i,j,k)/shelfdep )*vol/Jac(i,j,k)       
                                                                        
               s(i,j,k,n)= s(i,j,k,n) 
!c     &              -dtfac*fac*(s(i,j,k,n)-sfresh)                    
            end do 
         end do 
      end do 
                                                                        
      return 
                                                                        
      END                                           
