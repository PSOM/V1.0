SUBROUTINE mgrid(p,dtime,edt,cfcdiv)
  USE header, ONLY : NI,NJ,NK,maxout,maxint,int1,ngrid,rc_kind,user1
  !     --------------                                                    
  !     multigrid solver for the elliptic equation Del2 p= f              
  !     Run pgrogram "preprocess" to get the value of maxout              
  !     input :                                                           
  !     -----                                                             
  !     call the subroutine mgrid with p_initial on fine grid             
  !     maxout : max dimension on the one-dim array - includes outer point
  !     maxint: max dimension of the one-dim array for variables at interi
  !              points of grid                                           
  !     int1 :  number of internal points on the  fine grid               
  !     ngrid : number of grid levels                                     
  !     nx(m),ny(m),nz(m) = number of int grid points at level m. m=1...ng
  !     m=1 is the finest grid, m=ngrid refers to the coarsest grid.      
  !     ntout(m) is the total number of grid points (storage locations) at
  !               grid level m - including the outer ficticious points    
  !     ntint(m) is the total number of interior grid points (storage loca
  !               at grid level m                                         
  !     loco(m) is the storage location of a scalar variable in a one-dim 
  !             of dimension  Summation(m=1,ngrid) {ntout(m)}             
  !     loci(m) is the storage location of a scalar specified only at inte
  !             grid points                                               
  !     loccp(m) is the storage location for cp on the m-th level         
  !     res(i)  temporary stores the residual before it interpolated onto 
  !             the coarse grid                                           
  !     dxf,dyf,dzf  are the grid spacings on the fine grid               
  !     integer :: ngrid,m,maxout,maxint,int1,l,ncycle,ll
  INTEGER :: m,l,ncycle,ll
  INTEGER, PARAMETER :: nu1=2,nu2=1,noc=30

  INTEGER ntout(ngrid),ntint(ngrid),nx(ngrid),ny(ngrid),nz(ngrid),  &
       &     loco(ngrid),loci(ngrid),loccp(ngrid)                         
  REAL(kind=rc_kind) ::  cp(19*maxint),rhs(maxint),p(maxout),             &
       &     res(int1),dtime,tol,maxres,edt,ratio,oldres,cfcdiv
  REAL(kind=rc_kind) :: tol_r
  !     edt= EPS/dtime. If the tolerance is on (u_x +v_y +eps*w_z)        
  !     then  the tolerance on the residual r= edt*(u_x+ v_y + eps*w_z)   
  !     is equal to edt*(the specified tolerance)                         
  !
  !commented out by jbw
  !      open (unit=90, file='mg.in')
  !      read(90,*) nx(1),ny(1),nz(1),noc,nu1,nu2,tol
  !      close(90)

!  IF(ABS(user1-8)<0.05) then
!    tol_r=1d-10
!   else
   tol_r=1d-13
!  ENDIF

  nx(1)=NI
  ny(1)=NJ
  nz(1)=NK
  tol=tol_r
  IF (cfcdiv.LE.tol) THEN 
     DO m=1,maxout 
        p(m)= 0.d0 
     END DO
     maxres = edt*cfcdiv 
     ncycle = 0 
     GOTO 101 
  END IF

  !     redefine the tolerance as tol*edt                                 
  tol= edt*tol 
  ntint(1)= nx(1)*ny(1)*nz(1) 
  ntout(1)= (nx(1)+2)*(ny(1)+2)*(nz(1)+2) 
  !                                                                       
  loco(1)= 1 
  loci(1)= 1 
  loccp(1)= 1 
  DO m=2,ngrid 
     !check this at the beginning of the program instead of doing for every step. -jbw
     !    IF (MOD(nx(m-1),2).NE.0) GOTO 20 
     !    IF (MOD(ny(m-1),2).NE.0) GOTO 20 
     !    IF (MOD(nz(m-1),2).NE.0) GOTO 20 
     nx(m)= nx(m-1)/2 
     ny(m)= ny(m-1)/2 
     nz(m)= nz(m-1)/2 
     ntint(m)= nx(m)*ny(m)*nz(m) 
     ntout(m)= (nx(m)+2)*(ny(m)+2)*(nz(m)+2) 
     loco(m)= loco(m-1) + ntout(m-1) 
     loci(m)= loci(m-1) + ntint(m-1) 
     !     for the 19 point stencil                                          
     loccp(m)= loccp(m-1) + 19*ntint(m-1) 
  END DO
  !                                                                       
  !     compute fine grid coefficients and rhs on fine grid               
  CALL cpfine(dtime,cp,rhs) 
  !      call checkcp(1,nx(1),ny(1),nz(1),cp(loccp(1)))                   
  !      write(6,*) 'cpfine checked'                                      
  !     **********************************                                
  DO m=2,ngrid 
     CALL cpcors(nx(m),ny(m),nz(m),cp(loccp(m-1)),cp(loccp(m))) 
     !         call checkcp(m,nx(m),ny(m),nz(m),cp(loccp(m)))                
  END DO
  !                                                                       
  DO 1000 ncycle=1,noc 
     !+         if (mod(ncycle,3).eq.0) then                                 
     !+c     then we update the rhs                                          
     !+            call intermq(dtime)                                       
     !+            call remcor                                               
     !+            call newsrc                                               
     !+            call exitp                                                
     !+            call newqfn(dtime,rhs)                                    
     !+         endif                                                        
     !w         write(6,*) 'cycle=',ncycle                                   
     !     initialize the coarse grid values of p                            
     DO l=loco(2),maxout 
        p(l)= 0.d0 
     END DO
     !     DOWN CYCLE                                                        
     !     for m=1                                                           
     m=1 
     !     -----------------------                                           
     !wc     write out starting resid                                        
     !w      call resid(m,nx(m),ny(m),nz(m),cp(loccp(m)),                    
     !w     &        p(loco(m)),rhs(loci(m)),res,maxres)                     
     !w      write(6,*) 'starting resid=',maxres                             
     !     ---------------------------------                                 
     DO 11 l=1,nu1 
        CALL linerelax(nx(m),ny(m),nz(m),cp(loccp(m)),p(loco(m)),rhs(loci(m)))                                             
        !     *******try and put mgpfill here ******                            
        CALL mgpfill(dtime,p) 
        !         call mgpfill(dtime,p)                                         
11   ENDDO
     !c      call mgpfill(dtime,p)                                           
     CALL resid(m,nx(m),ny(m),nz(m),cp(loccp(m)),                      &
          &        p(loco(m)),rhs(loci(m)),res,maxres)                       
     !w      if (ncycle.ne.1) then                                           
     !w         ratio = maxres/oldres                                        
     !w         write(6,*) ncycle,'maxres',maxres,'conv ratio=',ratio        
     !w      endif                                                           
     !w      oldres= maxres                                                  
     !w      write(6,*) 'm= ',m, 'maxres=',maxres                            
     IF (maxres.LT.tol) GOTO 101 
     CALL restrict(nx(m),ny(m),nz(m),res,rhs(loci(m+1))) 
     !      do 100 m=1,ngrid-1                                               
     DO 100 m=2,ngrid-1 
        DO 21 l=1,nu1 
           CALL linerelax(nx(m),ny(m),nz(m),cp(loccp(m)),p(loco(m)),   &
                &           rhs(loci(m)))                                          
21      ENDDO
        CALL resid(m,nx(m),ny(m),nz(m),cp(loccp(m)),                   &
             &        p(loco(m)),rhs(loci(m)),res,maxres)                       
        !w         write(6,*) 'm= ',m, 'maxres=',maxres                         
        CALL restrict(nx(m),ny(m),nz(m),res,rhs(loci(m+1))) 
100  ENDDO
     !                                                                       
     !     UP CYCLE                                                          
     DO 200 m= ngrid,2,-1 
        IF (m.EQ.ngrid) THEN 
25         CALL sor(nx(m),ny(m),nz(m),cp(loccp(m)),p(loco(m)),rhs(loci(m)))                                          
        ELSE 
           DO 41 l=1,nu2 
              CALL linerelax(nx(m),ny(m),nz(m),cp(loccp(m)),p(loco(m)),rhs(loci(m)))                            
41         ENDDO
        ENDIF
        !     this call to resid is not necessary - temporary for res check     
        !w         call resid(m,nx(m),ny(m),nz(m),cp(loccp(m)),                 
        !w     &        p(loco(m)),rhs(loci(m)),res,maxres)                     
        !w         write(6,*) 'm= ',m, 'maxres=',maxres                         
        CALL efill(nx(m),ny(m),nz(m),p(loco(m))) 
        CALL prolong(nx(m),ny(m),nz(m),p(loco(m)),p(loco(m-1)) ) 
200  ENDDO
     !      do l=1,3                                                         
     !         call  mgpfill(dtime,p)                                        
     !      end do                                                           
1000 ENDDO
  !*      write(6,*) '100 mg V cycles, not converged'                     
  !     added Nov 17,1993                                                 
  m=1 
  DO 111 l=1,nu1 
     CALL linerelax(nx(m),ny(m),nz(m),cp(loccp(m)),p(loco(m)),      &
          &        rhs(loci(m)))                                             
     !     *******try and put mgpfill here ******                            
     CALL  mgpfill(dtime,p) 
111  enddo
     !                                                                       
     m=1 
     CALL resid(m,nx(m),ny(m),nz(m),cp(loccp(m)),                      &
          &     p(loco(m)),rhs(loci(m)),res,maxres)                          
     !w      write(6,*) 'max-res=',maxres                                    
     !     normliz can be called only for the case of dirichlet bcs.         
101  CONTINUE
     ! 101  call normliz(nx(1),ny(1),nz(1),p)                                
     !      do l=1,3                                                         
     !           call  mgpfill(dtime,p)                                      
     !      end do                                                           
     !     maxres/edt is the value of (u_x +v_y +ep*w_Z)                     
     maxres= maxres/edt 
!     WRITE(6,*) 'mg-cycle ', ncycle,'  maxres=',maxres 
     !     this call to resid is not necessary - temporary for res check     
     !      call resid(m,nx(m),ny(m),nz(m),cp(loccp(m)),                     
     !     &     p(loco(m)),rhs(loci(m)),res,maxres)                         
     !     GOTO 201 
     !20   WRITE(6,*) 'cannot coarsen this grid as specified' 

201  RETURN 


   END SUBROUTINE mgrid
