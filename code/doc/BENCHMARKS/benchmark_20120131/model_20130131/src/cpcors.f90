      subroutine cpcors(nxm,nym,nzm,cpf,cpc) 
!     ---------------------------------------                           
      USE header
!     call cpcors(m,nx(m),ny(m),nz(m),cp(loccp(m)))                     
!     Computes the coefficients in the elliptic operator by averaging   
!     the coefficients from the finer mesh.  Since the grid spacing     
!     is doubled the coefficients should be divided by 4, but instead   
!     we will multiply the rhs by 4.                                    
!     Uses the 19 point stencil (in 3d).                                
!                                                                       
      integer i,j,k,ic,jc,kc,nxm,nym,nzm 
      REAL(kind=rc_kind) ::                                                  &
     &     cpc(19,nxm,nym,nzm),cpf(19,2*nxm,2*nym,2*nzm)                
!                                                                       
      do 10 kc=1,nzm 
         k= 2*kc 
         do 20 jc=1,nym 
            j= 2*jc 
            do 30 ic=1,nxm 
               i= 2*ic 
               cpc(1,ic,jc,kc)= 0.125*(cpf(1,i,j,k)+cpf(1,i-1,j,k)+     &
     &              cpf(1,i,j-1,k) +cpf(1,i,j,k-1) +cpf(1,i-1,j-1,k)+   &
     &              cpf(1,i-1,j,k-1) +cpf(1,i,j-1,k-1) +                &
     &              cpf(1,i-1,j-1,k-1) )                                
               cpc(2,ic,jc,kc)= 0.25*(cpf(2,i,j,k)+cpf(2,i,j-1,k)+      &
     &              cpf(2,i,j,k-1) +cpf(2,i,j-1,k-1))                   
               cpc(3,ic,jc,kc)= 0.25*(cpf(3,i,j,k)+cpf(3,i-1,j,k)+      &
     &              cpf(3,i-1,j,k-1) +cpf(3,i,j,k-1) )                  
               cpc(4,ic,jc,kc)= 0.25*(cpf(4,i-1,j,k) +cpf(4,i-1,j-1,k)+ &
     &              cpf(4,i-1,j,k-1) +cpf(4,i-1,j-1,k-1) )              
               cpc(5,ic,jc,kc)= 0.25*(cpf(5,i,j-1,k) +cpf(5,i-1,j-1,k)+ &
     &              cpf(5,i,j-1,k-1) +cpf(5,i-1,j-1,k-1) )              
               cpc(6,ic,jc,kc)= 0.25*(cpf(6,i,j,k) +cpf(6,i-1,j,k)+     &
     &              cpf(6,i,j-1,k) +cpf(6,i-1,j-1,k) )                  
               cpc(7,ic,jc,kc)= 0.25*(cpf(7,i,j,k-1) +cpf(7,i-1,j,k-1)+ &
     &              cpf(7,i,j-1,k-1) +cpf(7,i-1,j-1,k-1) )              
               cpc(8,ic,jc,kc)= 0.5*(cpf(8,i-1,j,k) +                   &
     &              cpf(8,i-1,j,k-1))                                   
               cpc(9,ic,jc,kc)= 0.5*(cpf(9,i-1,j-1,k) +                 &
     &              cpf(9,i-1,j-1,k-1))                                 
               cpc(10,ic,jc,kc)= 0.5*(cpf(10,i,j-1,k) +                 &
     &              cpf(10,i,j-1,k-1))                                  
               cpc(11,ic,jc,kc)= 0.5*(cpf(11,i,j,k) +cpf(11,i,j,k-1)) 
               cpc(12,ic,jc,kc)= 0.5*(cpf(12,i-1,j,k-1) +               &
     &              cpf(12,i-1,j-1,k-1))                                
               cpc(13,ic,jc,kc)= 0.5*(cpf(13,i,j,k-1) +                 &
     &              cpf(13,i,j-1,k-1))                                  
               cpc(14,ic,jc,kc)= 0.5*(cpf(14,i,j,k) +cpf(14,i,j-1,k)) 
               cpc(15,ic,jc,kc)= 0.5*(cpf(15,i-1,j,k) +                 &
     &              cpf(15,i-1,j-1,k))                                  
               cpc(16,ic,jc,kc)= 0.5*(cpf(16,i,j-1,k-1) +               &
     &              cpf(16,i-1,j-1,k-1))                                
               cpc(17,ic,jc,kc)= 0.5*(cpf(17,i,j,k-1) +                 &
     &              cpf(17,i-1,j,k-1))                                  
               cpc(18,ic,jc,kc)= 0.5*(cpf(18,i,j,k) +cpf(18,i-1,j,k)) 
               cpc(19,ic,jc,kc)= 0.5*(cpf(19,i,j-1,k) +                 &
     &              cpf(19,i-1,j-1,k))                                  
   30       continue 
   20    continue 
   10 continue 
!                                                                       
      return 
      END                                           
