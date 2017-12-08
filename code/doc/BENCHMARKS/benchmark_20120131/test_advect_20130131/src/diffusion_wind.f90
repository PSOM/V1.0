subroutine diffusion(var,vardif,m,step,iv_compute_kzl) 
  !     ---------------------------------------------                     
  USE header
  !     Wind Stress can be specified here                                 
  !     Kz*du/dz = -Tx/rho,   Kz*dv/dz = -Ty/rho                          
  !     use level m                                                       
  !     computes d2s/dz2 at the cell centers.                             
  !     fricu and fricv are u,v diffusion terms (in header)               
  !     that need to be saved (if m.eq.0) for n2budget)                   

  implicit none 
  REAL(kind=rc_kind) :: dSin,posneg,facreverse
  integer i,j,k,m,step,nday
  integer iv_compute_kzl 
  REAL(kind=rc_kind) :: zdep,yw,yset,ymid,timday,stress,ustar,Ekdepth,ff                                      
                                                                        

  REAL(kind=rc_kind), dimension(0:NI+1,0:NJ+1, 0:NK+1) :: var 
  REAL(kind=rc_kind) :: vardif(NI,NJ,NK)

  REAL(kind=rc_kind) :: dvardzfc(0:NK)
  REAL(kind=rc_kind) :: Kdudzb,Kdvdzb,Kdudzt,Kdvdzt,wzkth,fact,Cp     
  REAL(kind=rc_kind) :: fac,facb,dfactor,Kzmax,facy,rhoinv,diff,stressy
  REAL(kind=rc_kind) :: wgt,day,ts,KzmaxTr,ztransit,zextent,thy                                         

  PARAMETER (Kzmax= 1.d-3, KzmaxTr=1.d-3) !     The viscosity is Kz*Kzmax                                         

  timday= dble(step)*dtf*TL/86400.d0 
  ts=  dble(step)/4000.d0 

  stressmax=  0.1d0                                                
  stressy= 0.d0 
                   
                                                    
  facb = RR*DL      ! Linear Drag                                                       
  ! facb = RR*DL*UL ! Quadratic drag                                                    

  fac= 1.d0/(UL*DL*delta) 
  fact = DL/UL 
                                                                        
  do j=1,NJ 
    do i=1,NI 

      if(iv_compute_kzl==1) then

        !-------------------
        ! COMPUTATION OF KZ
        !-------------------
        stress= sqrt(stressx(j)*stressx(j)+ stressy*stressy) 
        ustar = sqrt( stress /R0 ) 
        ff= ffc(i,j)*FPAR 
        Ekdepth= 0.4d0*ustar/ff 
        zextent= 0.5d0*Ekdepth !     Set zextent
        ztransit = -Ekdepth 
        do k=NK,0,-1 
          Kz(i,j,k)= 1.d0 
          zdep = zf(i,j,k)*DL 
          thy = (1.d0 +tanh(((zdep-ztransit)/zextent)*PI))*0.5d0                                              
          Kz(i,j,k)= max(0.01d0,thy) 
        end do
 
      endif

      !-------------------
      ! ADDITIONAL COMPUTATIONS
      !-------------------
                                                                      
      !     Quadratic drag                                                    
      !=            Kdudzb= facb*u(i,j,1,m)*abs(u(i,j,1,m))                  
      !=            Kdvdzb= facb*v(i,j,1,m)*abs(v(i,j,1,m))                  
      !     Linear Drag                                                       
      !Kdudzb= facb*u(i,j,1,m) 
      !Kdvdzb= facb*v(i,j,1,m) 
      !                                                              
      !rhoinv = 1.d0/rho(i,j,NK) 
      !Kdudzt= stressx(j)*rhoinv*fact 
      !Kdvdzt= stressy*rhoinv*fact 


      !-------------------
      ! COMPUTATION OF THE VARIATIONS OF THE VARIABLE ALONG Z
      !-------------------

      do k=1,NK-1 
        wzkth= wz(i,j,k) 
        dvardzfc(k)= wzkth*(var(i,j,k+1)-var(i,j,k)) 
      end do 
      dvardzfc(0)= 0.d0 
      dvardzfc(NK)= 0.0 


      !-------------------
      ! COMPUTATION OF THE FLUX DIVERGENCE
      !-------------------

      k=1 
!     ---
      dfactor=  fac*Jac(i,j,k)*wz(i,j,k) 
      vardif(i,j,k)= dfactor*KzmaxTr*(dvardzfc(k)*Kz(i,j,k)- dvardzfc(k-1)*Kz(i,j,k-1) )       

      do k=2,NK-1 
!     -----------
        dfactor=  fac*KzmaxTr*Jac(i,j,k)*wz(i,j,k) 
        vardif(i,j,k)= dfactor*(dvardzfc(k)*Kz(i,j,k)- dvardzfc(k-1)*Kz(i,j,k-1) )    
      end do 
                                                                       
      k=NK 
!     ----                                                        
      dfactor=  fac*Jac(i,j,k)*wz(i,j,k) 
      vardif(i,j,k)= dfactor*KzmaxTr*(dvardzfc(k)*Kz(i,j,k)- dvardzfc(k-1)*Kz(i,j,k-1) )       

    end do ! i
  end do ! j
                                                                        
return 
END                                           
                                                                        
                                                                        
      subroutine viscosity(dudz,dvdz,drdz,i,j) 
!     ---------------------------------------------                     
        USE header
!     Compute a Richardson number based vertical viscosity              
!     coefficient Kz using the final velocities u,v and density         
!     with the n=0 subscript                                            
!     Kz is situated at k cell faces                                    
!     The algorithm is from Rev Geophys., p 373, 1994. Large et al.     
!                                                                       
!     Assumes that rho is evaluated just prior to this call             
                                                                        
      integer i,j,k 
                                                                        
      integer n1 
      parameter (n1=3) 
      REAL(kind=rc_kind) :: fac,bvfreq,grho,RiCr,Ri,vshear 
      REAL(kind=rc_kind) :: dudz(0:NK),dvdz(0:NK),drdz(0:NK) 
      parameter (RiCr= 0.7d0) 
                                                                        
      grho= 9.81/R0 
!     fac= DL/(UL*UL)                                                   
      DLinv = 1.0/DL 
      fac = UL*UL/(DL*DL) 
!                                                                       
!     Set kz(k=0)= 1. It is required only for the s,T equations, since  
!     for the momentum equations,  K*du/dz= ru.                         
!     Kz is not needed at k=0 and NK since stress bcs are used.         
      Kz(i,j,0)= 0.d0 
      Kz(i,j,NK)= 0.d0 
!      if  (i.eq.16) write(100,*) 'j = ', j                             
      do k=1,NK-1 
         bvfreq=  -grho*drdz(k)*DLinv 
!     BVfreq is N**2 and is in s^-2 if re-dim by DLinv                  
         vshear= (dudz(k)*dudz(k) + dvdz(k)*dvdz(k))*fac 
         if (vshear.eq.0.d0) then 
!         if (vshear.le.1.d-12) then                                    
            Ri= 100.d0 
         else 
            Ri= bvfreq/vshear 
         end if 
!     unstable density profile => Ri < 0                                
         if (Ri.le.0.d0) then 
            Kz(i,j,k)= 1.d0 
         else if (Ri.lt.RiCr) then 
            Kz(i,j,k)= (1.d0 - (Ri*Ri)/(RiCr*RiCr))**n1 
         else 
            Kz(i,j,k)= 0.d0 
         endif 
!         if ((i.eq.NI/2).and.(j.eq.NJ/2)) write(6,*) Ri ,bvfreq,vshear,
!     &        Kz(i,j,k)                                                
      end do 
                                                                        
!     The value of Kz to be used is  Kz(k) *Kzmax                       
                                                                        
      return 
      END                                           
