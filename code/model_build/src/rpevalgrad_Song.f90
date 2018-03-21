      subroutine rpevalgrad_Song(n) 
!      subroutine rpflatbtm(n)                                          
!     ---------------------------------------------                     
  USE header
!  USE rpgrads
!     This routine evaluates rp and its gradients for the case          
!     where d sigma/dx, d sigma/dy are 0, which is approx true for      
!     the flat-bottommed domain.  When topography is present use        
!     the routine  rpevalgrad.f  instead. rpevalgrad.f is slow, and     
!     hence we have written this for the special case of no topography. 
!                                                                       
!     CHANGED FOR THE INCOMPRESSIBLE CASE where rho=potrho(s,T,p=atmos) 
!     Modified for the free surface case. This routine now evaluates    
!     r or the component of the pressure associated with the deviation  
!     of density from the mean.                                         

      integer i,j,k,n 
      REAL(kind=rc_kind) :: z,hpd,zt,zb,                                      &
     &     temp,rg,const,dep,ht,R02,rpxi,rpeta                          
!      REAL(kind=rc_kind) :: findz                                           
!      REAL(kind=rc_kind) :: drpx(NI,NJ,NK),drpy(NI,NJ,NK),                   &
!     &     grpifc(0:NI,NJ,NK),grpjfc(NI,0:NJ,NK)                        
!      common/rpgrads/drpx,drpy,grpifc,grpjfc 
!      if (.NOT(rect)) then 
      if (rect.eqv.(.false.)) then 
!        need to modify grpifc,grpjfc to contain cross diff terms       
         write(6,*) 'modify grpifc,grpjfc' 
         stop 
      end if 
!     Given the salinity and temperature field we evaluate the pot densi
!     field by calling potdens. rho is the dimensionless density \rho'  
!     rho= rho(s,T),                                                    
!                                                                       
!      const= G*gpr/P1                                                  
      const= 10.d0*gpr/P1 
      call findzall 
      call evalrho(rho,n) 
!                                                                       
      do 10 i=0,NI+1 
         do 20 j=0,NJ+1 
            hpd= h(i,j)*HDL +D(i,j) 
            dep= D(i,j) 
            ht= h(i,j)*HDL 
!     for k=NK (firt for the toppmost layer of cells)                   
!     z is the dimensional z and is therefore multiplied by DL          
            k= NK 
            z= DL*zc(i,j,k) 
            zt= DL*zf(i,j,k) 
            zb= DL*zf(i,j,k-1) 
            rg= (rho(i,j,k) - R0)*const 
!            rg= (rhreal - R0)*const *R0/rhreal -non-Boussinesq for this
            rp(i,j,k)= rg*(zt-z) 
            temp= rg*(z-zb) 
!                                                                       
!****************                                                       
!            rp(i,j,k)= 0.d0                                            
!*******************                                                    
!     for the rest of the column                                        
!     remember that findz will not work for the k=0 cell where sig is -v
!     We will simply fill rp at k=0 using extrapolation - since this    
!     is what we use to fill s and T in any case.                       
!            do 30 k=NK-1,0,-1                                          
            do 30 k=NK-1,1,-1 
!     z is the dimensional z and is therefore multiplied by DL          
!               z= DL*findz(dfloat(k)-0.5,dztop,dep,ht,NK)              
               z= DL*zc(i,j,k) 
               zt= zb 
               zb= DL*zf(i,j,k-1) 
               rg= (rho(i,j,k) - R0)*const 

!     rg= (rhreal - R0)*const *R0/rhreal-non-Boussinesq for this term   
               rp(i,j,k)= rg*(zt-z) +temp +rp(i,j,k+1) 
               temp= rg*(z-zb) 
   30       continue 
!-            rp(i,j,0)= 2.d0*rp(i,j,1) -rp(i,j,2)                      
   20    continue 
   10 continue 
!                                                                       

!     Boundary conditions                                               
      do k=1,NK 
         do i=1,NI 
            rp(i,0,k)= rp(i,1,k) 
            rp(i,NJ+1,k)= rp(i,NJ,k) 
         end do 
!     periodic-ew boundaries                                            
         do j=0,NJ+1 
            rp(NI+1,j,k)= rp(1,j,k) 
            rp(0,j,k)= rp(NI,j,k) 
         end do 
      end do 
!                                                                       
      do 40 k=1,NK 
         do 50 j=1,NJ 
            do 60 i=1,NI 
               rpxi= 0.5d0*(rp(i+1,j,k) -rp(i-1,j,k)) 
               rpeta= 0.5d0*(rp(i,j+1,k) -rp(i,j-1,k)) 
               drpx(i,j,k)= rpxi*ux(i,j) + rpeta*vx(i,j) 
               drpy(i,j,k)= rpxi*uy(i,j) + rpeta*vy(i,j) 
   60       continue 
   50    continue 
         do 70 j=1,NJ 
            do 75 i=0,NI 
               grpifc(i,j,k)= (rp(i+1,j,k) -rp(i,j,k))*gi(i,j,k,1) 
   75       continue 
   70    continue 
         do 80 j=0,NJ 
            do 85 i=1,NI 
               grpjfc(i,j,k)= (rp(i,j+1,k) -rp(i,j,k))*gj(i,j,k,2) 
   85       continue 
   80    continue 
   40 continue 
!                                                                       
      return 
      END                                           
