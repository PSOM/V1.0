FUNCTION  potdens(sl,Tl)
  !-----------------------------------------------------------------------
  !     computes the POTENTIAL density from the equation of state         
  !     Uses the 1 atmos. International eqn of state for sea water (1980) 
  !     and the high press. international eqn of state for sea water (19
  !     (Unesco technical Papers in Marine Science 38)                    
  !                                                                       
  !     s is the practical salinity ie. (mass of solute/mass of water)*10^
  !     T is temp in  deg C, p is pressure in bar                         
  !     r0= rho0(s,T,p=0), rw= density of standard mean ocean water       
  !                
      USE header, ONLY : rc_kind                                                       
      implicit none
      REAL(kind=rc_kind) :: potdens
      REAL(kind=rc_kind) ::  sl,Tl,p,rho0,r0,rw,Aw,Bw,Kw,a,b,K0,K,s15
      !                                                                       
      s15= sl*sqrt(sl) 
      !     One atmos. International eqn of state for sea water (1980)        
      !     ----------------------------------------------------------        
      rw= 999.842594d0 +  6.793952d-2*Tl - 9.095290d-3*Tl**2 +            &
           1.001685d-4*Tl**3 -1.120083d-6*Tl**4 +                         &
           6.56332d-9*Tl**5                                               
      potdens= rw + (8.24493d-1 - 4.0899d-3*Tl + 7.6438d-5*Tl**2          &
           -8.2467d-7*Tl**3 + 5.3875d-9*Tl**4 )*sl +                       &
           (-5.72466d-3 + 1.0227d-4*Tl -1.6546d-6*Tl**2)*s15 +            &
           4.8314d-4*sl*sl                                                 
!                                                                       
!     High press. international eqn of state for sea water (1980)       
!     ----------------------------------------------------------        
!-      Kw= 19652.21d0 + 148.4206*T - 2.327105d0*T**2 +1.360477d-2*T**3 
!-     &    - 5.155288d-5*T**4                                          
!-      Aw= 3.239908  + 1.43713d-3*T + 1.16092d-4*T**2 - 5.77905d-7*T**3
!-      Bw= 8.50935d-5 - 6.12293d-6*T + 5.2787d-8*T**2                  
!-c                                                                     
!-      K0= Kw + (54.6746 -0.603459*T + 1.09987d-2*T**2 -               
!-     &    6.167d-5*T**3)*s + (7.944d-2 + 1.6483d-2*T                  
!-     &    - 5.3009d-4*T**2)*s15                                       
!-      a= Aw + (2.2838d-3 - 1.0981d-5*T -1.6078d-6*T**2)*s +           
!-     &   1.91075d-4*s15                                               
!-      b= Bw + (-9.9348d-7 + 2.0816d-8*T + 9.1697d-10*T**2)*s          
!-c                                                                     
!-      K= K0 + a*p + b*p*p                                             
!-c                                                                     
!-      rho0= r0/(1.d0-p/K)                                             
!-c     ------------------                                              
      return 
    end FUNCTION potdens
