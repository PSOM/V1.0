subroutine  sigma2z(x,y,sigma,z)
  !     -------------------------------------------------------           
  USE header, only: D,HDL,dztop,NK,h,pfac,rc_kind
  !     finds the value of z (non-dim by DL) given the value of sigma.    
  !     ht and dep are non-dim by DL.                                     
  !     At every (i,j), the column is divided into NK equal-depth cells.  
  REAL(kind=rc_kind) :: sigma,dnkm1,hpd,xfac,x,y,di,dj,h1,h2,z                     
  integer :: i,j
  !    pfac is the stretching in z. higher pfac gives more points near surf.

  pfac= 2.0d0 
  !      pfac= 5.d0  !c=NK32,dztop=0.5m                                   
  !==      pfac= 4.d0  !c=NK32,dztop=0.5m USED FOR SEVERAL MLI runs       
  !=      pfac=3.d0  !c used with NK=32 for first set of runs             
  dnkm1= dble(NK-1) 

  i = INT(x)
  j = INT(y)
  di = x - i
  dj = y - j
        
        if (sigma>=dnkm1) then
           !In the surface layer
           h1 = (h(i+1,j)-h(i,j))*di + h(i,j)
           h2 = (h(i+1,j+1) - h(i,j+1))*di + h(i,j+1)
           hpd = (h2 - h1) * dj + h1                                              
           hpd= hpd*HDL +dztop 
           z = ( sigma - dnkm1 ) * hpd - dztop
        else
           !below the surface layer
           xfac = (dnkm1 - sigma ) / dnkm1 
           z    = (exp(pfac*xfac)-1.d0)* &
                  (D(i,j) + dztop) / (exp(pfac)-1d0) - dztop                               
        endif
  
  return 
end subroutine sigma2z
