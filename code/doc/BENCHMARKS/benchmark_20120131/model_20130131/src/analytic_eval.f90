real(kind=rc_kind) function analytic_eval(mld, z, nml, n1, ndep, w) 
  !      -----------------------------------------------------            
  !      DESCRIPTION                                                      
                                                                        
  !      Evaluates background density profile according                   
  !      to analytic buoyancy frequency function.                         
                                                                        
  !      __________________                                               
  !      ------------------                                               
  !      SET VARIABLES                                                    
     use header, only : rc_kind

  implicit none                                                                        
  REAL(kind=rc_kind), parameter ::PI=3.14159265358979323846
  REAL(kind=rc_kind) ::  mld, z, nml, n1, ndep, w, a, b, c, e, d,f,g,h,erf
  !      __________________                                               
  !      ------------------                                               
  !      EVALUATE FUNCTION                                                
                                                                        
       h = (3.14159265358979323846**(0.5)) 
       g = (mld + 3*w - z)/w 
       a = z*ndep 
       b = w*(log(tanh(g) + 1) + (2*z)/w) 
       c = ndep/2 - nml/2 
       d = h*n1*erf(g) 
       e = (2/w) 
                                                                        
       f = a - (b * c) - (d / e) 
                                                                        
                                                                        
       analytic_eval = (1026.0d0/9.8d0) * f + 1026.0 
!      ___________________                                              
                                                                        
       return 
      END                                           
