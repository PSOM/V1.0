subroutine mgpfill(dtime,pf) 
!     -----------------------------------                               
  USE header
!     modified for periodicew bc                                        
!     --------------------------                                        
!     In this version of mgpfill, the dpdsig term is left out           
!     at the vertical faces because g13 g23 are larger than             
!     g11 and g12.                                                      
!     call mgpfill(dtime,p(loco(m)))                                    
!     computes the pressure at points outside the impermeable boundaries
!     by using central difference to specify the pressure gradient      
!     at the boundary which is known from the condition v_normal=0.     
!     We use uf,vf,wf at n=0  as the BC on uf,vf,wf.  This is ok for    
!     the case where the BCs are not changing witth time.               
!     We should ideally use one-sided differencing near the edges and co
!     so as to not use these values which we are filling just by averagi
!     Then we should also modify cpfine so as to not access these values
!     Right now we will not bother to do this.                          
!     If the BCs are changing with time we must modify this subroutine. 
      integer i,j,k 
      REAL(kind=rc_kind) :: dtime,edt,dpdxi,dpdeta,dpdsig,pxi,peta,psig,     &
     &     pf(0:NI+1,0:NJ+1,0:NK+1),gin                                 
!                                                                       
      if (dtime.gt.1.d-16) then 
         edt= EPS/dtime 
      else 
!         write(6,*) 'dtime 0, set edt to 0 in mgpfill'                 
         edt = 0.0 
      end if 
                                                                        
                                                                        
!      delinv= 1/delta, delta= D/L                                      
!                                                                       
      if (qpr.ne.0) then 
                                                                        
!     face k=0                                                          
!     --------                                                          
      k=0 
      do 50 i=1,NI 
         do  50 j=1,NJ 
            gin= 1.d0/gqk(i,j,k,3) 
!     note that quantities are not defined at i=0,NI+1,j=0,NJ+1, and    
!     so we revert to filling pf inside the the domain only.            
            dpdxi= 0.25*( pf(i+1,j,k+1) +pf(i+1,j,k)                    &
     &           -pf(i-1,j,k+1) -pf(i-1,j,k) )                          
            dpdeta= 0.25*( pf(i,j+1,k+1) +pf(i,j+1,k)                   &
     &           -pf(i,j-1,k+1) -pf(i,j-1,k) )                          
            psig= gin*(-skfc(i,j,k)                                     &
     &           - gqk(i,j,k,1)*dpdxi - gqk(i,j,k,2)*dpdeta +           &
     &           edt*(czf(i,j,k)-wfbcb(i,j)) )                          
!     since  : wf(i,j,k)= czf(i,j,k) -dte*delta*(pz +skfc(i,j,k))       
            pf(i,j,0)= pf(i,j,1) - psig 
   50 continue 
                                                                        
!     bottom edges                                                      
      j=0 
      do i=1,NI 
         gin= 1.d0/gqk(i,j,0,3) 
         dpdxi= 0.25*( pf(i+1,j,k+1) +pf(i+1,j,k)                       &
     &        -pf(i-1,j,k+1) -pf(i-1,j,k) )                             
         dpdeta= 0.5*( pf(i,j+1,k+1) +pf(i,j+1,k)                       &
     &        -pf(i,j,k+1) -pf(i,j,k) )                                 
         psig= gin*(-skfc(i,j,k)                                        &
     &        - gqk(i,j,k,1)*dpdxi - gqk(i,j,k,2)*dpdeta )              
         pf(i,j,0)= pf(i,j,1) - psig 
      end do 
      j=NJ+1 
      do i=1,NI 
         gin= 1.d0/gqk(i,j,0,3) 
         dpdxi= 0.25*( pf(i+1,j,k+1) +pf(i+1,j,k)                       &
     &        -pf(i-1,j,k+1) -pf(i-1,j,k) )                             
         dpdeta= 0.5*( pf(i,j,k+1) +pf(i,j,k)                           &
     &        -pf(i,j-1,k+1) -pf(i,j-1,k) )                             
         psig= gin*(-skfc(i,j,k)                                        &
     &        - gqk(i,j,k,1)*dpdxi - gqk(i,j,k,2)*dpdeta )              
         pf(i,j,0)= pf(i,j,1) - psig 
      end do 
!                                                                       
!     face k=NK                                                         
!     --------                                                          
      k=NK 
      do 60 i=1,NI 
         do  60 j=0,NJ+1 
            pf(i,j,NK+1)= -pf(i,j,NK) 
!            pf(i,j,NK+1)= 0.                                           
   60 continue 
!                                                                       
!      pf(0,0,0)= (pf(1,0,0)+pf(0,1,0)+pf(0,0,1))/3.d0                  
!      pf(NI+1,0,0)= (pf(NI,0,0)+pf(NI+1,1,0)+pf(NI+1,0,1))/3.d0        
!      pf(NI+1,NJ+1,0)=(pf(NI,NJ+1,0)+pf(NI+1,NJ,0)+pf(NI+1,NJ+1,1))/3.d
!      pf(0,NJ+1,0)= (pf(1,NJ+1,0)+pf(0,NJ,0)+pf(0,NJ+1,1))/3.d0        
                                                                        
                                                                        
!     faces j=0,j=NJ                                                    
!     --------------                                                    
      j=0 
      do 30 i=1,NI 
         do 31 k=1,NK 
            gin= 1.d0/gqj(i,j,k,2) 
            dpdxi= 0.25*( pf(i+1,j+1,k) +pf(i+1,j,k)                    &
     &           -pf(i-1,j+1,k) -pf(i-1,j,k) )                          
            dpdsig= 0.25*( pf(i,j+1,k+1) +pf(i,j,k+1)                   &
     &           -pf(i,j+1,k-1) -pf(i,j,k-1) )                          
            peta= gin*(                                                 &
     &           - gqj(i,j,k,1)*dpdxi - gqj3(i,j,k)*dpdsig              &
     &           +(edt*(cyf(i,j,k)-vfbcs(i,k))-sjfc(i,j,k)) )           
!c     &           - gqj(i,j,k,1)*dpdxi                                 
            pf(i,0,k)= pf(i,1,k) - peta 
   31    continue 
   30 continue 
!                                                                       
      j=NJ 
      do 40 i=1,NI 
         do 41 k=1,NK 
            gin= 1.d0/gqj(i,j,k,2) 
            dpdxi= 0.25*( pf(i+1,j+1,k) +pf(i+1,j,k)                    &
     &           -pf(i-1,j+1,k) -pf(i-1,j,k) )                          
            dpdsig= 0.25*( pf(i,j+1,k+1) +pf(i,j,k+1)                   &
     &           -pf(i,j+1,k-1) -pf(i,j,k-1) )                          
            peta= gin*(                                                 &
     &           - gqj(i,j,k,1)*dpdxi - gqj3(i,j,k)*dpdsig              &
     &           +(edt*(cyf(i,j,k)-vfbcn(i,k))-sjfc(i,j,k)) )           
!c     &           - gqj(i,j,k,1)*dpdxi                                 
            pf(i,NJ+1,k)= pf(i,NJ,k) + peta 
   41    continue 
   40 continue 
!                                                                       
                                                                        
      else 
!     face k=NK                                                         
!     --------                                                          
      k=NK 
      do i=1,NI 
         do j=0,NJ+1 
            pf(i,j,NK+1)= -pf(i,j,NK) 
         end do 
      end do 
      end if 
!     the above is done only if qpr is non-zero                         
!                                                                       
!     faces i=0,i=NI : periodicew bcs                                   
!     -------------------------------                                   
      i=0 
      do 10 j=0,NJ+1 
         do 11 k=0,NK+1 
            pf(0,j,k)= pf(NI,j,k) 
            pf(NI+1,j,k)= pf(1,j,k) 
   11    continue 
   10 continue 
!                                                                       
!                                                                       
!                                                                       
      return 
      END                                           
