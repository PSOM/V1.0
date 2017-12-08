
SUBROUTINE momentum(pcorr,step) 
  !----------------------------------------------------                   
  USE header
!  use relaxation
  INTEGER step,i,j,k 
  REAL(kind=rc_kind) ::  dtim,fdiv,cfcdiv,ctrdiv,edt,pcorr(maxout),temp 
  REAL(kind=rc_kind) :: tarray(2),tim 
  !     --------------------------------------------------------          
  !     initialized with u(t=0) and uf(t=0), ie. u(n=0),uf(n=0)           
  !     ------------------------------------------------------            
  !     get ready for the 1st step                                        
  !     evalrp and srcface were callled at the finish or the previous step
  !     or in init                                                        
  !     save old value of h  in array oldh
  oldh= h

  !     1st R-K step                                                      
  !     -----------------                                                 
  dtim= dtf/3.d0 
  !     compute new zc,zf and wz,wx,wy and Jac (since grid has moved)     


  CALL findzall 
  CALL sigma 
  !      call outflowbc(0,dtim)                                           
  !     advecn(m,n,dtim) means use level m and write s,T to level n       
  CALL advecn(0,1,dtim,step) 
  !     call intpol with the time level at which uf is known              
  !     besides interpolating the convective terms, intpol computes u*    
  CALL intpol 
  !     we call the following 3 after advecn, in which case               
  !     evalrp would use the new value of s,T(n+1), but old value of h(n) 


 if(.NOT.(use_Shchepetkin)) then ! calculates baroclinic press gradients
    CALL rpevalgrad(1) 
   else
    CALL rpevalgrad(1) 
    CALL rpevalgrad_Sche(1);
! #include "test_Sche.f90"
    drpx(1:NI,1:NJ,1:NK)=ru4_Sche(1:NI,1:NJ,1:NK);  drpy(1:NI,1:NJ,1:NK)=rv4_Sche(1:NI,1:NJ,1:NK); grpifc(0:NI,1:NJ,1:NK)=ru2_Sche(0:NI,1:NJ,1:NK); grpjfc(1:NI,0:NJ,1:NK)=rv2_Sche(1:NI,0:NJ,1:NK);
 endif

  !      call rpflatbtm(1)                                                
  !      call viscous(0)                                                  
  CALL coriolis(0) 
  CALL srcface(0,step)      ! puts coriolis terms on faces
  CALL hsolve(h,oldh,hdt,dtim)        ! solve for h
  CALL calcskfc                       ! puts on faces
  CALL vhydro(dtim)                    ! calculates hydrostatic vel (predictor)
  CALL cfdiv(cfcdiv)   

  ! v,w, present beacuse cannot bal both v and vf in geostrophic perfectly
  !     in calling newcor, n should be the same as in coriolis            
  CALL newcor(dtim,0)                ! update to coriolis -not being used
  CALL newsrc                        ! puts on faces
  !     compute q (correction to NH press p)                              
  !     cpuf must be called before mgrid which calls vface and overwrites 
  edt= EPS/dtim 
  CALL mgrid(pcorr,dtim,edt,cfcdiv)  ! solves the non-hydrostatic pressure correction  

  !      call pfilter(pcorr)                                              
  !c      call linesolve(p,dtim,step,edt)                                 
  !      call psolve(p,dtim,step,edt)                                     
  !c      call pnsolve(p,dtim,step,edt)                                   

  CALL vface(pcorr,dtim)           ! final face velocities
  !     call vcenter with the value of n that we wish to fill             
  CALL vcenter(pcorr,dtim,1)       ! final grid center velocities
  !     computes the final vel (at n+1 -th time step) at the cell centers 
  !     Now correct the NH pressure                                       
  IF (fnhhy.NE.0) CALL pcorrect(pcorr)    ! correct the non-hyd pressure
  CALL facediv(dtim,fdiv)              !check divergence
  CALL cdiv(dtim,ctrdiv,1)            
  !      write(6,*)  'rk-1 cfcdiv ',cfcdiv,' fdiv ',fdiv,' cdiv ',ctrdiv 

  !     2nd R-K step                                                      
  !     -----------------                                                 
  dtim= 0.5d0*dtf 
  !     compute new zc,zf and wz,wx,wy and Jac (since grid has moved)     
  CALL findzall 
  CALL sigma 
  !     advecn(m,n,dtim) means use level m and write s,T to level n       
  CALL advecn(1,1,dtim,step) 
  !     call intpol with the time level at which uf is known              
  !     besides interpolating the convective terms, intpol computes u*    
  CALL intpol 
  !     we call the following 3 after advecn, in which case               
  !     evalrp would use the new value of s,T(n+1), but old value of h(n) 
  !      tim= dtime(tarray)                                               

 if(.NOT.(use_Shchepetkin)) then ! calculates baroclinic press gradients
    CALL rpevalgrad(1) 
   else
    CALL rpevalgrad(1) 
    CALL rpevalgrad_Sche(1);
! #include "test_Sche.f90"
    drpx(1:NI,1:NJ,1:NK)=ru4_Sche(1:NI,1:NJ,1:NK);  drpy(1:NI,1:NJ,1:NK)=rv4_Sche(1:NI,1:NJ,1:NK); grpifc(0:NI,1:NJ,1:NK)=ru2_Sche(0:NI,1:NJ,1:NK); grpjfc(1:NI,0:NJ,1:NK)=rv2_Sche(1:NI,0:NJ,1:NK);
 endif

  CALL coriolis(1) 
  CALL srcface(1,step) 

  CALL hsolve(h,oldh,hdt,dtim) 
  CALL calcskfc 
  CALL vhydro(dtim) 
  CALL cfdiv(cfcdiv) 
  !     in calling newcor, n should be the same as in coriolis.           

  CALL newcor(dtim,1) 
  CALL newsrc 
  edt= EPS/dtim 
  CALL mgrid(pcorr,dtim,edt,cfcdiv) 
  !     ********************                                              
  CALL vface(pcorr,dtim) 
  !     call vcenter with the value of n that we wish to fill             
  CALL vcenter(pcorr,dtim,1) 
  !     computes the final vel (at n+1 -th time step) at the cell centers 
  !     Now correct the NH pressure                                       
  IF (fnhhy.NE.0) CALL pcorrect(pcorr) 
  CALL facediv(dtim,fdiv) 
  CALL cdiv(dtim,ctrdiv,1) 

  !     3rd R-K step                                                      
  !     -----------------                                                 
  dtim= dtf 
  !     compute new zc,zf and wz,wx,wy and Jac (since grid has moved)     
  CALL findzall 
  CALL sigma 
  !     advecn(m,n,dtim) means use level m and write s,T to level n       
  CALL advecn(1,0,dtim,step) 
  !     call intpol with the time level at which uf is known              
  !     besides interpolating the convective terms, intpol computes u*    
  CALL intpol 


 if(.NOT.(use_Shchepetkin)) then ! calculates baroclinic press gradients
    CALL rpevalgrad(0) 
   else
    CALL rpevalgrad(0) 
    CALL rpevalgrad_Sche(0);
! #include "test_Sche.f90"
    drpx(1:NI,1:NJ,1:NK)=ru4_Sche(1:NI,1:NJ,1:NK);  drpy(1:NI,1:NJ,1:NK)=rv4_Sche(1:NI,1:NJ,1:NK); grpifc(0:NI,1:NJ,1:NK)=ru2_Sche(0:NI,1:NJ,1:NK); grpjfc(1:NI,0:NJ,1:NK)=rv2_Sche(1:NI,0:NJ,1:NK);
 endif

  CALL coriolis(1) 
  CALL srcface(1,step) 

  CALL hsolve(h,oldh,hdt,dtim) 
  CALL calcskfc 
  !c      call hnsolve(h,oldh,hdt,dtim)                                   
  CALL vhydro(dtim) 
  CALL cfdiv(cfcdiv) 
  CALL newcor(dtim,1) 
  CALL newsrc 
  edt= EPS/dtim 
  CALL mgrid(pcorr,dtim,edt,cfcdiv) 
  !     *********************                                             
  CALL vface(pcorr,dtim) 

  !     computes the vel at the cell faces                                

  !     call vcenter with the value of n that we wish to fill             
  CALL vcenter(pcorr,dtim,0) 
  !     computes the final vel (at n+1 -th time step) at the cell centers 
  !     Now correct the NH pressure (remains 0 if hy)                     
  IF (fnhhy.NE.0) CALL pcorrect(pcorr) 
  CALL facediv(dtim,fdiv) 
  CALL cdiv(dtim,ctrdiv,0) 
  !     Convective adjustment carried out every 9 time steps.             
  !      if (mod(step,9).eq.0) then                                       
  CALL evalrho(rho,0) 
  CALL conadjust(step,0) 
  !      end if                                                           
  CALL calcn2 
  !     This was eliminated for non reactive tracer - ONR expts.          
  ! CALL tracersource(0,dtim)   
  !     call eddydivergence(dtim) - needs cbar
  !                                                                       
!  CALL sponge(0,dtf)

END SUBROUTINE momentum
