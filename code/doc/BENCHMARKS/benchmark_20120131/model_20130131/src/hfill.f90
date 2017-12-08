subroutine hfill(dtimel,hf) 
!     -----------------------------------                               
  USE header
!     MODIFIED FOR PERIODICEW                                           
!     -----------------------                                           
!     call hfill(dtime,h)                                               
!     computes h at points outside the boundaries                       
!     by using central difference to specify the h gradient             
!     at the boundary which is known from the condition v_normal=0.     
!     We use uf,vf,wf at n=0  as the BC on uf,vf,wf.  This is ok for    
!     the case where the BCs are not changing witth time.               
!     We make the approx that tilde-u_nonrmal = uf_normal at the boundar
!     We should ideally use one-sided differencing near the edges and co
!     so as to not use these values which we are filling just by averagi
!     Then we should also modify cpfine so as to not access these values
!     Right now we will not bother to do this -  but will call this rout
!     a few times (alternating with hsolve) to iterate when g12 is non-z
!     If the BCs are changing with time we must modify this subroutine. 
      integer i,j,k,l 
      REAL(kind=rc_kind) :: dtimel,edt,sumcxf,sumuf,sumsif,sumcyf,sumvf,      &
     &     sumsjf,hxi,heta,ginv,hf(0:NI+1,0:NJ+1),cons,dnk,gkn,         &
     &     gradh,sumhxn,sumhyn,sumgi(2),sumgj(2)                        
!                                                                       
      edt= EPS/dtimel 
      dnk= dble(NK) 
!      cons= (1.d0 -kappah)*dnk*gpr                                     
!      gkn= gpr*kappah*dnk                                              
      cons= (1.d0 -kappah)*gpr 
      gkn= gpr*kappah 
!                                                                       
!     faces j=0,j=NJ                                                    
!     --------------                                                    
      do 30 i=1,NI 
!         ginv= 1.d0/gj(i,0,2)                                          
!         if (i.eq.1) then                                              
!            hxi= 0.5d0*(hf(i+1,0)+hf(i+1,1) -hf(i,0)-hf(i,1))          
!         else if (i.eq.NI) then                                        
!            hxi= 0.5d0*(hf(i,0)+hf(i,1) -hf(i-1,0)-hf(i-1,1))          
!         else                                                          
            hxi= 0.25*(hf(i+1,0)+hf(i+1,1) -hf(i-1,0)-hf(i-1,1)) 
!         end if                                                        
         sumcyf= 0.d0 
         sumvf= 0.d0 
         sumsjf= 0.d0 
         sumhyn= 0.d0 
         do l=1,2 
            sumgj(l)= 0.d0 
         end do 
         do 35 k=1,NK 
            sumcyf= sumcyf +cyf(i,0,k) 
            sumvf= sumvf +vfbcs(i,k) 
            sumsjf= sumsjf +sjfc(i,0,k) 
            sumhyn= sumhyn +hyn(i,0,k) 
            do l=1,2 
               sumgj(l)= sumgj(l) + gj(i,0,k,l) 
            end do 
   35    continue 
         gradh= (edt*(sumcyf -sumvf) -sumsjf -cons*sumhyn)/gkn 
         hf(i,0)= hf(i,1) + (sumgj(1)*hxi -gradh)/sumgj(2) 
   30 continue 
                                                                        
      do 40 i=1,NI 
!         ginv= 1.d0/gj(i,NJ,2)                                         
!         if (i.eq.1) then                                              
!            hxi= 0.5d0*(hf(i+1,NJ)+hf(i+1,NJ+1) -hf(i,NJ)-hf(i,NJ+1))  
!         else if (i.eq.NI) then                                        
!            hxi= 0.5d0*(hf(i,NJ)+hf(i,NJ+1) -hf(i-1,NJ)-hf(i-1,NJ+1))  
!         else                                                          
            hxi= 0.25*(hf(i+1,NJ)+hf(i+1,NJ+1)-hf(i-1,NJ)-hf(i-1,NJ+1)) 
!         end if                                                        
         sumcyf= 0.d0 
         sumvf= 0.d0 
         sumsjf= 0.d0 
         sumhyn= 0.d0 
         do l=1,2 
            sumgj(l)= 0.d0 
         end do 
         do 45 k=1,NK 
            sumcyf= sumcyf +cyf(i,NJ,k) 
            sumvf= sumvf +vfbcn(i,k) 
            sumsjf= sumsjf +sjfc(i,NJ,k) 
            sumhyn= sumhyn +hyn(i,NJ,k) 
            do l=1,2 
               sumgj(l)= sumgj(l) +gj(i,NJ,k,l) 
            end do 
   45    continue 
         gradh= (edt*(sumcyf -sumvf) -sumsjf -cons*sumhyn)/gkn 
         hf(i,NJ+1)= hf(i,NJ) + (-sumgj(1)*hxi +gradh)/sumgj(2) 
   40 continue 
                                                                        
!                                                                       
!     faces i=0,i=NI                                                    
!     --------------                                                    
      do j=0,NJ+1 
         hf(0,j)= hf(NI,j) 
         hf(NI+1,j)= hf(1,j) 
      end do 
!                                                                       
      return 
      END                                           
