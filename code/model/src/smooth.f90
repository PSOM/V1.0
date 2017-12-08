subroutine smooth 
#include "cppdefs.h"
!     ------------------------------------------------                  
  USE header
!     smooths out D and computes new Ddx,Ddy                            
!                                                                       
      integer i,j,ncyc 
      REAL(kind=rc_kind) ::  Ddu,Ddv 
!c     FLAT BOTTOM                                                      
!      do j=0,NJ+1                                                      
!         do i=0,NI+1                                                   
!           D(i,j)= 1.5                                                 
!           Ddx(i,j)= 0.d0                                              
!           Ddy(i,j)= 0.d0                                              
!        end do                                                         
!      end do                                                           
!      goto 999                                                         
!---------------------------                                            
!     If smoothing of topography is required                            
!c      do 100 ncyc= 1,4                                                
!c         do 10 j=0,NJ+1                                               
!c            do 10 i=1,NI                                              
!c               D(i,j)= 0.25*(D(i+1,j)+D(i-1,j)) +0.5d0*D(i,j)         
!c 10      continue                                                     
!c         do 20 j=1,NJ                                                 
!c            do 20 i=0,NI+1                                            
!c               D(i,j)= 0.25*(D(i,j+1)+D(i,j-1)) +0.5d0*D(i,j)         
!c 20      continue                                                     
!c 100  continue                                                        
!                                                                       
      do j=1,NJ 
         do i=1,NI 
            Ddu= 0.5d0*(D(i+1,j) -D(i-1,j)) 
            Ddv= 0.5d0*(D(i,j+1) -D(i,j-1)) 
            Ddx(i,j)= Ddu*ux(i,j) +Ddv*vx(i,j) 
            Ddy(i,j)= Ddu*uy(i,j) +Ddv*vy(i,j) 
        enddo
      enddo
      i=0 
      do 40 j=1,NJ 
         Ddu= D(i+1,j) -D(i,j) 
         Ddv= 0.5d0*(D(i,j+1) -D(i,j-1)) 
         Ddx(i,j)= Ddu*ux(i,j) +Ddv*vx(i,j) 
         Ddy(i,j)= Ddu*uy(i,j) +Ddv*vy(i,j) 
   40 continue 
      i=NI+1 
      do 50 j=1,NJ 
         Ddu= D(i,j) -D(i-1,j) 
         Ddv= 0.5d0*(D(i,j+1) -D(i,j-1)) 
         Ddx(i,j)= Ddu*ux(i,j) +Ddv*vx(i,j) 
         Ddy(i,j)= Ddu*uy(i,j) +Ddv*vy(i,j) 
   50 continue 
      j=0 
      do 60 i=1,NI 
         Ddu= 0.5d0*(D(i+1,j) -D(i-1,j)) 
         Ddv= D(i,j+1) -D(i,j) 
         Ddx(i,j)= Ddu*ux(i,j) +Ddv*vx(i,j) 
         Ddy(i,j)= Ddu*uy(i,j) +Ddv*vy(i,j) 
   60 continue 
      j=NJ+1 
      do 70 i=1,NI 
         Ddu= 0.5d0*(D(i+1,j) -D(i-1,j)) 
         Ddv= D(i,j) -D(i,j-1) 
         Ddx(i,j)= Ddu*ux(i,j) +Ddv*vx(i,j) 
         Ddy(i,j)= Ddu*uy(i,j) +Ddv*vy(i,j) 
   70 continue 
!                                                                       
      Ddx(0,0)= 0.5d0*(Ddx(1,0)+Ddx(0,1)) 
      Ddy(0,0)= 0.5d0*(Ddy(1,0)+Ddy(0,1)) 
      Ddx(NI+1,0)= 0.5d0*(Ddx(NI,0)+Ddx(NI+1,1)) 
      Ddy(NI+1,0)= 0.5d0*(Ddy(NI,0)+Ddy(NI+1,1)) 
      Ddx(NI+1,NJ+1)= 0.5d0*(Ddx(NI,NJ+1)+Ddx(NI+1,NJ)) 
      Ddy(NI+1,NJ+1)= 0.5d0*(Ddy(NI,NJ+1)+Ddy(NI+1,NJ)) 
      Ddx(0,NJ+1)= 0.5d0*(Ddx(1,NJ+1)+Ddx(0,NJ)) 
      Ddy(0,NJ+1)= 0.5d0*(Ddy(1,NJ+1)+Ddy(0,NJ)) 
!                                   
#ifdef file_output                                    
      open (unit=70, file=TRIM(dirout)//'Ddxdy.dat') 
      do j=0,NJ+1 
         write(70,*) 'j= ',j 
         do i=0,NI+1 
            write(70,*) D(i,j),Ddx(i,j),Ddy(i,j) 
        end do 
      end do 
      close(70) 
#endif
                                                                        
      return 
      END                                           
