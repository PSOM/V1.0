      SUBROUTINE DGTSL (N, C, D, E, B, INFO) 
!***BEGINPROLOGUE  DGTSL                                                
!***PURPOSESolve a tridiagonal linear system.                           
!***LIBRARYSLATEC (LINPACK)                                             
!***CATEGORYD2A2A                                                       
!***TYPEDOUBLE PRECISION (SGTSL-S, DGTSL-D, CGTSL-C)                    
!***KEYWORDSLINEAR ALGEBRA, LINPACK, MATRIX, SOLVE, TRIDIAGONAL         
!***AUTHORDongarra, J., (ANL)                                           
!***DESCRIPTION                                                         
!                                                                       
!     DGTSL given a general tridiagonal matrix and a right hand         
!     side will find the solution.                                      
!                                                                       
!     On Entry                                                          
!                                                                       
!     N       INTEGER                                                   
!     is the order of the tridiagonal matrix.                           
!                                                                       
!     C       DOUBLE PRECISION(N)                                       
!     is the subdiagonal of the tridiagonal matrix.                     
!     C(2) through C(N) should contain the subdiagonal.                 
!     On output C is destroyed.                                         
!                                                                       
!     D       DOUBLE PRECISION(N)                                       
!     is the diagonal of the tridiagonal matrix.                        
!     On output D is destroyed.                                         
!                                                                       
!     E       DOUBLE PRECISION(N)                                       
!     is the superdiagonal of the tridiagonal matrix.                   
!     E(1) through E(N-1) should contain the superdiagonal.             
!     On output E is destroyed.                                         
!                                                                       
!     B       DOUBLE PRECISION(N)                                       
!     is the right hand side vector.                                    
!                                                                       
!     On Return                                                         
!                                                                       
!     B       is the solution vector.                                   
!                                                                       
!     INFO    INTEGER                                                   
!     = 0 normal value.                                                 
!     = K if the K-th element of the diagonal becomes                   
!     exactly zero.  The subroutine returns when                        
!     this is detected.                                                 
!                                                                       
!***REFERENCESJ. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.       
!     Stewart, LINPACK Users' Guide, SIAM, 1979.                        
!***ROUTINESCALLED  (NONE)                                              
!***REVISIONHISTORY  (YYMMDD)                                           
!     780814  DATE WRITTEN                                              
!     890531  Changed all specific intrinsics to generic.  (WRB)        
!     890831  Modified array declarations.  (WRB)                       
!     890831  REVISION DATE from Version 3.2                            
!     891214  Prologue converted to Version 4.0 format.  (BAB)          
!     900326  Removed duplicate information from DESCRIPTION section.   
!     (WRB)                                                             
!     920501  Reformatted the REFERENCES section.  (WRB)                
!***ENDPROLOGUE  DGTSL                                                  
      INTEGER N,INFO 
      DOUBLE PRECISION C(*),D(*),E(*),B(*) 
!                                                                       
      INTEGER K,KB,KP1,NM1,NM2 
      DOUBLE PRECISION T 
!***FIRSTEXECUTABLE STATEMENT  DGTSL                                    
      INFO = 0 
      C(1) = D(1) 
      NM1 = N - 1 
      IF (NM1 .LT. 1) GO TO 40 
      D(1) = E(1) 
      E(1) = 0.0D0 
      E(N) = 0.0D0 
!                                                                       
      DO 30 K = 1, NM1 
         KP1 = K + 1 
!                                                                       
!     FIND THE LARGEST OF THE TWO ROWS                                  
!                                                                       
         IF (ABS(C(KP1)) .LT. ABS(C(K))) GO TO 10 
!                                                                       
!     INTERCHANGE ROW                                                   
!                                                                       
         T = C(KP1) 
         C(KP1) = C(K) 
         C(K) = T 
         T = D(KP1) 
         D(KP1) = D(K) 
         D(K) = T 
         T = E(KP1) 
         E(KP1) = E(K) 
         E(K) = T 
         T = B(KP1) 
         B(KP1) = B(K) 
         B(K) = T 
   10    CONTINUE 
!                                                                       
!     ZERO ELEMENTS                                                     
!                                                                       
         IF (C(K) .NE. 0.0D0) GO TO 20 
         INFO = K 
         GO TO 100 
   20    CONTINUE 
         T = -C(KP1)/C(K) 
         C(KP1) = D(KP1) + T*D(K) 
         D(KP1) = E(KP1) + T*E(K) 
         E(KP1) = 0.0D0 
         B(KP1) = B(KP1) + T*B(K) 
   30 END DO 
   40 CONTINUE 
      IF (C(N) .NE. 0.0D0) GO TO 50 
      INFO = N 
      GO TO 90 
   50 CONTINUE 
!                                                                       
!     BACK SOLVE                                                        
!                                                                       
      NM2 = N - 2 
      B(N) = B(N)/C(N) 
      IF (N .EQ. 1) GO TO 80 
      B(NM1) = (B(NM1) - D(NM1)*B(N))/C(NM1) 
      IF (NM2 .LT. 1) GO TO 70 
      DO 60 KB = 1, NM2 
         K = NM2 - KB + 1 
         B(K) = (B(K) - D(K)*B(K+1) - E(K)*B(K+2))/C(K) 
   60 END DO 
   70 CONTINUE 
   80 CONTINUE 
   90 CONTINUE 
  100 CONTINUE 
!                                                                       
      RETURN 
      END                                           
