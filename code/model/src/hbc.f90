subroutine hbc(chf,fn,dtimel) 
!-------------------------------------------------------------------c   
  USE header
!     Modified for periodicew boundaries                                
!     ----------------------------------                                
!     evaluates the boundaries in psolve c fills in                     
!     ch0,chim1,chip1,chjm1,chjp1,chkm1,chkp1,fn for the boundaries     
!     3-D Poisson solver - uses SOR                                     
!     Written for the control volume formulation                        
!     Uses the 19 point stencil (in 3d).                                
!                                                                       
!     fn is the contribution to the source terms on the RHS.            
!     Whenever possible we substitue the the velocities themselves      
!     (at the boundaries). We cannot substitute the u,w vels at the entr
!     Our eqn is                                                        
!     chim1*p(i-1) +chip1*p(i+1) +chjm1p(j-1) ... -ch0*p(i,j,k)= f      
   integer i,j,k,l,im1 
   REAL(kind=rc_kind) :: dtimel,edt,chf(9,NI,NJ),fn(NI,NJ),dtinv,sumuf,sumvf,sumsif,sumsjf,                   &
  &     sumcxf,sumcyf,gprinv,eg,const,kaph1,sumhxn,sumgi(2),sumhyn,sumgj(2)                              

   kaph1= 1.d0 -kappah 
   edt= EPS/dtimel 
   eg= edt/(gpr*kappah) 
   gprinv= 1.d0/(gpr*kappah) 
   dtinv= HDL/(dtimel*kappah) 
   !     dnk= dfloat(NK)                                                   
   const= kaphinv -1.d0 

   !     For i=2...NI-1, j=2...NJ-1, k=2,NK-1  it is straightforward       

   do i=1,NI 
     do j=1,NJ 
       if ((j.eq.1).or.(j.eq.NJ)) then 
         fn(i,j)= eg*( kaphinv*wfbcb(i,j)-J2d(i,j)*oldh(i,j)*dtinv )                         
         chf(1,i,j)= -eg*J2d(i,j)*dtinv 
         do l=2,9 
            chf(l,i,j)= 0.d0 
         end do 
         sumsif = 0.d0 
         sumcxf = 0.d0 
         sumuf  = 0.d0 
         sumhxn = 0.d0 
         do l=1,2 
            sumgi(l)= 0.d0 
         end do 
         do k=1,NK 
           sumsif = sumsif + sifc(i,j,k) 
           sumcxf = sumcxf + cxf(i,j,k) 
           sumuf  = sumuf  + uf(i,j,k) 
           sumhxn = sumhxn + hxn(i,j,k) 
           do l=1,2 
             sumgi(l)= sumgi(l) + gi(i,j,k,l) 
           end do 
         end do 
         fn(i,j)= fn(i,j) + eg*(sumcxf +const*sumuf) - gprinv*sumsif - const*sumhxn                        
         chf(1,i,j)= chf(1,i,j) -      sumgi(1) 
         chf(2,i,j)= chf(2,i,j) +      sumgi(1) 
         chf(4,i,j)= chf(4,i,j) + 0.25*sumgi(2) 
         chf(5,i,j)= chf(5,i,j) - 0.25*sumgi(2) 
         chf(6,i,j)= chf(6,i,j) + 0.25*sumgi(2) 
         chf(8,i,j)= chf(8,i,j) - 0.25*sumgi(2) 
                                                                  
         sumsif = 0.d0 
         sumcxf = 0.d0 
         sumuf  = 0.d0 
         sumhxn = 0.d0 
         do l=1,2 
            sumgi(l)= 0.d0 
         end do 
         if (i.eq.1) then 
           im1=NI 
          else 
           im1= i-1 
         endif 
         do k=1,NK 
           sumsif = sumsif + sifc(im1,j,k) 
           sumcxf = sumcxf + cxf(im1,j,k) 
           sumuf  = sumuf  + uf(im1,j,k) 
           sumhxn = sumhxn + hxn(im1,j,k) 
           do l=1,2 
             sumgi(l)= sumgi(l) +gi(im1,j,k,l) 
           end do 
         enddo
         fn(i,j)= fn(i,j) - eg*(sumcxf +const*sumuf) + gprinv*sumsif + const*sumhxn                      
         chf(1,i,j)= chf(1,i,j) -      sumgi(1) 
         chf(3,i,j)= chf(3,i,j) +      sumgi(1) 
         chf(4,i,j)= chf(4,i,j) - 0.25*sumgi(2) 
         chf(5,i,j)= chf(5,i,j) + 0.25*sumgi(2) 
         chf(7,i,j)= chf(7,i,j) - 0.25*sumgi(2) 
         chf(9,i,j)= chf(9,i,j) + 0.25*sumgi(2) 
                                                                   
         if (j.eq.NJ) then 
            sumvf= 0.d0 
            do k=1,NK 
               sumvf= sumvf +vfbcn(i,k) 
            enddo
            fn(i,j)= fn(i,j) +eg*sumvf*kaphinv 
         else 
            sumsjf= 0.d0 
            sumcyf= 0.d0 
            sumvf= 0.d0 
            sumhyn= 0.d0 
            do l=1,2 
              sumgj(l)= 0.d0 
            end do 
            do k=1,NK 
              sumsjf = sumsjf + sjfc(i,j,k) 
              sumcyf = sumcyf + cyf(i,j,k) 
              sumvf  =  sumvf + vf(i,j,k) 
              sumhyn = sumhyn + hyn(i,j,k) 
              do l=1,2 
                sumgj(l)= sumgj(l) +gj(i,j,k,l) 
              end do 
            enddo
            fn(i,j)= fn(i,j) + eg*(sumcyf +const*sumvf) - gprinv*sumsjf - const*sumhyn                   
            chf(1,i,j)= chf(1,i,j) -      sumgj(2) 
            chf(2,i,j)= chf(2,i,j) + 0.25*sumgj(1) 
            chf(3,i,j)= chf(3,i,j) - 0.25*sumgj(1) 
            chf(4,i,j)= chf(4,i,j) +      sumgj(2) 
            chf(6,i,j)= chf(6,i,j) + 0.25*sumgj(1) 
            chf(7,i,j)= chf(7,i,j) - 0.25*sumgj(1) 
         endif 

         if (j.eq.1) then 
            sumvf= 0.d0 
            do k=1,NK 
               sumvf= sumvf +vfbcs(i,k) 
            enddo  
            fn(i,j)= fn(i,j) -eg*sumvf*kaphinv 
           else 
            sumsjf= 0.d0 
            sumcyf= 0.d0 
            sumvf= 0.d0 
            sumhyn= 0.d0 
            do l=1,2 
               sumgj(l)= 0.d0 
            end do 
            do k=1,NK 
              sumsjf= sumsjf + sjfc(i,j-1,k) 
              sumcyf= sumcyf + cyf(i,j-1,k) 
              sumvf = sumvf  + vf(i,j-1,k) 
              sumhyn= sumhyn + hyn(i,j-1,k) 
              do l=1,2 
                 sumgj(l)= sumgj(l) +gj(i,j-1,k,l) 
              end do 
            end do 
            fn(i,j)= fn(i,j) - eg*(sumcyf +const*sumvf) + gprinv*sumsjf + const*sumhyn                   
            chf(1,i,j)= chf(1,i,j) -      sumgj(2) 
            chf(2,i,j)= chf(2,i,j) - 0.25*sumgj(1) 
            chf(3,i,j)= chf(3,i,j) + 0.25*sumgj(1) 
            chf(5,i,j)= chf(5,i,j) +      sumgj(2) 
            chf(8,i,j)= chf(8,i,j) - 0.25*sumgj(1) 
            chf(9,i,j)= chf(9,i,j) + 0.25*sumgj(1) 
         endif 
       endif 
     enddo
   enddo
!                                                                       
      return 
      END                                           
