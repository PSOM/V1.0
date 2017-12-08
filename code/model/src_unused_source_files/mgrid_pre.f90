subroutine mgrid_pre()
use header, only : NI,NJ,NK,ngrid,maxout,maxint,int1
!     ------------------                                                
!     multigrid solver for the elliptic equation Del2 p= f              
!     Run pgrogram "preprocess" to get the value of maxdim              
!     input :                                                           
!     -----                                                             
!     maxdim : max dimension on the one-dim array - includes outer point
!     maxint: max dimension of the one-dim array for variables at interi
!              points of grid                                           
!     ngrid : number of grid levels                                     
!     nx(m),ny(m),nz(m) = number of int grid points at level m. m=1...ng
!     m=1 is the finest grid, m=ngrid refers to the coarsest grid.      
!     ntout(m) is the total number of grid points (storage locations) at
!               grid level m - including the outer ficticious points    
!     ntint(m) is the total number of interior grid points (storage loca
!               at grid level m                                         
!     nbc is the total number boundary location values                  
!      integer ngrid,m,maxdim,maxint,nbc,ntemp 
      integer m,nbc,ntemp 
!      parameter (ngrid=4)
!=      parameter (ngrid=5) 
!      parameter (ngrid=4)                                              
      integer nx(ngrid),ny(ngrid),nz(ngrid),ntout(ngrid),ntint(ngrid) 
!      open (unit=90, file='mg.in') 
!      read(*,*) nx(1),ny(1),nz(1)
      nx(1)=NI
      ny(1)=NJ
      nz(1)=NK
      print*, 'Number of grid points on fine grid: nx,ny,nz',        &
     &     nx(1),ny(1),nz(1)                                            
!      close(90) 
      ntint(1)= nx(1)*ny(1)*nz(1) 
      ntout(1)= (nx(1)+2)*(ny(1)+2)*(nz(1)+2) 
      do m=2,ngrid 
         if (mod(nx(m-1),2).ne.0) goto 20 
         if (mod(ny(m-1),2).ne.0) goto 20 
         if (mod(nz(m-1),2).ne.0) goto 20 
         nx(m)= nx(m-1)/2 
         ny(m)= ny(m-1)/2 
         nz(m)= nz(m-1)/2 
         ntint(m)= nx(m)*ny(m)*nz(m) 
         ntout(m)= (nx(m)+2)*(ny(m)+2)*(nz(m)+2) 
      end do 
      nbc= 0 
      do m=1,ngrid 
         nbc= nbc + 2*( nx(m)*ny(m) +nx(m)*nz(m) + ny(m)*nz(m)) 
         ntemp= 2*( nx(m)*ny(m) +nx(m)*nz(m) + ny(m)*nz(m)) 
         print*, 'm,ntint,ntout,nbc(m)',m,ntint(m),ntout(m),ntemp 
      end do 
!                                                                       
      maxout= 0 
      do m=1,ngrid 
         maxout= maxout + ntout(m) 
         maxint= maxint + ntint(m) 
      enddo 
        int1=NI*NJ*NK
      print*, 'maxout = ',maxout, ' maxint = ',maxint,'  nbc =',nbc, 'int1=',int1
      goto 21 
   20 print*, 'cannot coarsen this grid as specified' 
   21 stop 
 end subroutine  mgrid_pre              
