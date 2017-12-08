subroutine wind_stress(udif,vdif,step) 
  !     ---------------------------------------------                     
  USE header

  implicit none 
  integer i,j,k,m,step,n
  integer, parameter :: nw=1460
                                                                        
  REAL(kind=rc_kind) :: udif(NI,NJ,NK),vdif(NI,NJ,NK)

  REAL(kind=rc_kind) :: Kdudzt,Kdvdzt,fact,fac,rhoinv

  REAL(kind=rc_kind) :: xstressmax, ystressmax, edgeshelf, ywind, yice, edgeice
  REAL(kind=rc_kind) :: yrday,startday,del,d1,wt2,wt1,ycenter,ymax
  REAL(kind=rc_kind) :: ice(NI,NJ), yrdaywind(nw), stress(2,nw)

  real(kind=rc_kind) :: yw , yset, y0
  integer            :: iyw, iy0

  real :: stressmax

  udif=0.;vdif=0.;


!*********************************
! COMPUTATION OF THE WIND STRESS

  stress_top_x=0.;stress_top_y=0; 
  stressmax= 0.d0  
  
      if (step.eq.1) then
         open(unit=40,file='./BoB/2007windstress.in')
         do n=1,nw
            read(40,*) yrdaywind(n),stress(1,n),stress(2,n)
         end do
         close(40)
      end if

      startday=1.d0
      yrday= dble(step)*dtf*TL/86400.d0 + startday
      do n=1,nw-1
         if ((yrday.gt.yrdaywind(n)).and.(yrday.le.yrdaywind(n+1))) then
            del= yrdaywind(n+1)-yrdaywind(n)
            d1= yrday-yrdaywind(n)
            wt2= d1/del
            wt1=1.d0 - wt2
            xstressmax= stress(1,n)*wt1 + stress(1,n+1)*wt2
            ystressmax= stress(2,n)*wt1 + stress(2,n+1)*wt2
            goto 13
         end if
      end do


13      ycenter = 0.5*(yc(NJ)+yc(1)) 
      ymax = yc(NJ-10)-ycenter 
      edgeshelf=0.03
      ywind= 25.0 ! 25 km from coast
!!!!  KEEP ICE FIXED
      yice=(yc(NJ)-25.d0)

! reduce tightnness by 3/4
      edgeice=0.03

      do j=1,NJ
         do i=1,NI
            ice(i,j)= 0.5*(tanh(edgeice*(yc(j)-yice)*PI) + 1.d0)
         end do
      end do

      do j=1,NJ
         if (yc(j).lt.(2.*ywind)) then
            do i=1,NI
            stress_top_x(i,j)= xstressmax*0.5*(tanh(edgeshelf*(yc(j)-ywind)*PI)+1.d0)
            stress_top_y(i,j)= ystressmax*0.5*(tanh(edgeshelf*(yc(j)-ywind)*PI)+1.d0)
            end do
         else
            do i=1,NI
            stress_top_x(i,j)= xstressmax*(1.d0-ice(i,j))
            stress_top_y(i,j)= ystressmax*(1.d0-ice(i,j))
            end do
         end if
      end do

!************************************
! COMPUTATION OF THE SOURCE TERM


  fac= 1.d0/(UL*DL*delta) 
  fact = DL/UL 

  do j=1,NJ 
    do i=1,NI 
      stress_top(i,j) = sqrt(stress_top_x(i,j)*stress_top_x(i,j)+ stress_top_y(i,j)*stress_top_y(i,j))
      !stress_top(i,j)=0
      rhoinv = 1.d0/rho(i,j,NK) 
      !rhoinv = 1.d0/R0 
      Kdudzt= stress_top_x(i,j)*rhoinv*fact 
      Kdvdzt= stress_top_y(i,j)*rhoinv*fact 

      udif(i,j,NK)= fac*Jac(i,j,NK)*wz(i,j,NK)*Kdudzt   
      vdif(i,j,NK)= fac*Jac(i,j,NK)*wz(i,j,NK)*Kdvdzt   
    
    end do ! i
  end do ! j
                                                                        
return 
END                                           
                                                                        
                                                                        
