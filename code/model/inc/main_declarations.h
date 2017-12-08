

! Various declarations. Would need some clearing up... 
  INTEGER :: step,n,i,j,k,nblock,ik,inew,iold
  INTEGER ::m
 
  REAL(kind=rc_kind) ::  fdiv,ctrdiv,hsum    
  REAL(kind=rc_kind) ::  tarray(2), tim 
  REAL(kind=rc_kind) ::  pcorr(maxout)
  character(len=10) :: stepchar                                       

