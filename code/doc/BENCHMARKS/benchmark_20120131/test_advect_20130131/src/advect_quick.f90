subroutine advecn(m,n,dtimel,step) 
!     ---------------------------------------------                     
USE header
!     convecn(m,n,dtime) means use level m and write s,T to level n     
!     computes Cx,Cy,Cz at the cell centers.                            
!     uflx,vflx,wflx,sflx,Tflx are the fluxes defined at                
!     cell faces (using QUICK).                                         
!                                                                       

INTEGER :: i,i2,j,k,n,m,step

REAL(kind=rc_kind), PARAMETER  :: Eighth=0.125d0,Half=0.5
REAL(kind=rc_kind), PARAMETER :: courinv=2.d0

REAL(kind=rc_kind) :: absvel,left,right,ctr,fc,upos,uneg 

REAL(kind=rc_kind) :: dtimel

REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1) :: var      ! var=(s,T,u,v,w,Tr(it,:,:,:)                    
REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1) :: uvarx    ! uvarx  is the source term, divergence of the advective fluxes                   
REAL(kind=rc_kind), dimension(    0:NI  ,0:NJ  , 0:NK  ) :: varflx   ! varflx is a temporary variable to compute fluxes at the faces
REAL(kind=rc_kind), dimension(    1:NI  ,1:NJ  , 1:NK  ) :: vardif   ! vardif is the source term from diabatic processes      
REAL(kind=rc_kind) :: dvarx(0:NI+1),dvary(0:NJ),dvarz(0:NK+1)        ! dvarX is d(var)/dX at the faces (linear difference).

INTEGER :: selvar
INTEGER :: iv_compute_kz

do selvar=1,5+ntr
  if(selvar==1) var=T(:,:,:,m)
  if(selvar==2) var=s(:,:,:,m)
  if(selvar==3) var=u(:,:,:,m)
  if(selvar==4) var=v(:,:,:,m)
  if(selvar==5) var=w(:,:,:,m)
  do it=1,ntr
    if(selvar==5+it) var=Tr(it,:,:,:,m)
  enddo

                       
  !------------------------------------------
  ! computation of the source term in x.
                                                                                                     
  do k=1,NK 
    do j=1,NJ 

      dvarx(0)=var(1,j,k)-var(NI,j,k)
      do i=1,NI-1 
        dvarx(i)= var(i+1,j,k) - var(i,j,k) 
      end do 
      dvarx(NI)=dvarx(0)
      dvarx(NI+1)=dvarx(1)

      var(NI+1,j,k)= var(1,j,k) 
      var(0,j,k)= var(NI,j,k) 
       
      do i=1,NI
        if (uf(i,j,k).gt.0.d0) then 
          left= Eighth*(dvarx(i) -dvarx(i-1)) 
          ctr= Half* (var(i+1,j,k) +var(i,j,k)) 
          fc= ctr -left 
          call ultim(var(i-1,j,k),var(i,j,k),var(i+1,j,k),dvarx(i),dvarx(i-1),courinv,fc)         
          varflx(i,j,k)= uf(i,j,k)*fc 
         else 
          if(i==NI) then                  
            i2=0;
           else
            i2=i;
          endif
          right= Eighth*(dvarx(i2+1) -dvarx(i2)) 
          ctr= Half* (var(i2+1,j,k) +var(i,j,k)) 
          fc= ctr - right 
          call ultim(var(i2+2,j,k),var(i2+1,j,k),var(i,j,k),dvarx(i),dvarx(i2+1),courinv,fc)           
          varflx(i,j,k)= uf(i,j,k)*fc 
        end if 
      end do
      varflx(0,j,k)= varflx(NI,j,k) 
    end do 
  end do 

  do i=1,NI 
    uvarx(i,1:NJ,1:NK)=varflx(i,1:NJ,1:NK) -varflx(i-1,1:NJ,1:NK) 
  end do 

  ! computation of the source term in x done.
  !------------------------------------------


                                                       
  !------------------------------------------
  ! computation of the source term in y.
                                                               
  do k=1,NK 
    do i=1,NI
 
      dvary(0)= 0.d0; 
      dvary(NJ)= 0.d0; 
      do j=1,NJ-1 
        dvary(j)= var(i,j+1,k) - var(i,j,k) 
      end do 

      var(i,0,k)= var(i,1,k) 
      var(i,NJ+1,k)= var(i,NJ,k) 

      do j=1,NJ-1 
        if (vf(i,j,k).ge.0.d0) then 
          left= Eighth*(dvary(j) -dvary(j-1)) 
          ctr= Half* (var(i,j+1,k) +var(i,j,k)) 
          fc= ctr -left 
          call ultim(var(i,j-1,k),var(i,j,k),var(i,j+1,k),dvary(j),dvary(j-1),courinv,fc)         
          varflx(i,j,k)= vf(i,j,k)*fc 
        else 
         right= Eighth*(dvary(j+1) -dvary(j)) 
         ctr= Half* (var(i,j+1,k) +var(i,j,k)) 
         fc= ctr - right 
         call ultim(var(i,j+2,k),var(i,j+1,k),var(i,j,k),dvary(j),dvary(j+1),courinv,fc)           
         varflx(i,j,k)= vf(i,j,k)*fc 
       end if 
     end do 
     varflx(i,0,k)=0.d0;
     varflx(i,NJ,k)=0.d0;
    end do 
  end do 
      

  do j=1,NJ 
    uvarx(1:NI,j,1:NK)= (varflx(1:NI,j,1:NK) -varflx(1:NI,j-1,1:NK) ) + uvarx(1:NI,j,1:NK) 
  end do 
            

  ! computation of the source term in y done.
  !------------------------------------------



  !------------------------------------------
  ! computation of the source term in z.
                                                           
  do j=1,NJ 
    do i=1,NI 

      var(i,j,0)= var(i,j,1) 

      dvarz(0)= 0.d0 
      do k=1,NK-1 
         dvarz(k)= var(i,j,k+1) - var(i,j,k) 
      end do 
      dvarz(NK:NK+1)= 0.d0  !     Top boundary - Zero gradient                                      

      var(i,j,NK+1)= var(i,j,NK) 
                                                                  
      do k=1,NK 
        if (wf(i,j,k).ge.0.d0) then 
          left= Eighth*(dvarz(k) -dvarz(k-1)) 
          ctr= Half* (var(i,j,k+1) +var(i,j,k)) 
          fc= ctr -left 
          call ultim(var(i,j,k-1),var(i,j,k),var(i,j,k+1),dvarz(k),dvarz(k-1),courinv,fc)         
          varflx(i,j,k)= wf(i,j,k)*fc 
         else 
          if (k.eq.NK) then 
            varflx(i,j,k)= wf(i,j,k)*Half*(var(i,j,k+1) +var(i,j,k))                    
           else
            right= Eighth*(dvarz(k+1) -dvarz(k)) 
            ctr= Half* (var(i,j,k+1) +var(i,j,k)) 
            fc= ctr - right 
            call ultim(var(i,j,k+2),var(i,j,k+1),var(i,j,k),dvarz(k),dvarz(k+1),courinv,fc)           
            varflx(i,j,k)= wf(i,j,k)*fc 
          end if 
        end if 
      end do 
      varflx(i,j,0)= 0.d0   !     Solid Boundary                                                    
    end do 
  end do 

  do k=1,NK 
    uvarx(1:NI,1:NJ,k)=  ( varflx(1:NI,1:NJ,k) -varflx(1:NI,1:NJ,k-1) ) + uvarx(1:NI,1:NJ,k) 
  enddo 

  ! computation of the source term in z done.
  !------------------------------------------


  if(selvar<5.5) then
    call viscous(var,vardif); uvarx(1:NI,1:NJ,1:NK)=uvarx(1:NI,1:NJ,1:NK)-vardif(1:NI,1:NJ,1:NK)
  endif
                     
  if(selvar<4.5) then
    if(selvar==1) then
       iv_compute_kz=1 
      else
       iv_compute_kz=0
    endif
    call diffusion(var,vardif,m,step,iv_compute_kz); uvarx(1:NI,1:NJ,1:NK)=uvarx(1:NI,1:NJ,1:NK)-vardif(1:NI,1:NJ,1:NK)
  endif

  
  if(selvar==1) T(1:NI,1:NJ,1:NK,n)=T(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
  if(selvar==2) s(1:NI,1:NJ,1:NK,n)=s(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
  if(selvar==3) cx(1:NI,1:NJ,1:NK) =u(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
  if(selvar==4) cy(1:NI,1:NJ,1:NK) =v(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
  if(selvar==5) cz(1:NI,1:NJ,1:NK) =w(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
  do it=1,ntr
    if(selvar==5+it) Tr(it,1:NI,1:NJ,1:NK,n) =Tr(it,1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
  enddo


enddo !selvar




return 
                                                                        
END                                           
      
!----------------------------------------------------------------------
! SUBROUTINE ultim                                                                  
!     Implementation of the ULTIMATE limiter based on upstrm,dnstrm,ctr 
!     points and d(phi)/dx at the up- and dnstrm locations              
!     Returns the "ulitimate limited" face value of the variable        
!     that should be multiplied by uf,vf or wf to give the flux.        
                                                                         
SUBROUTINE ultim(up,ct,dn,gdn,gup,courinv,fc) 
USE header, only : rc_kind
                                                                  
REAL(kind=rc_kind) :: up,ct,dn,gup,gdn,courinv,fc,ref,del,adel,acurv 
                                                                  
del= dn -up 
adel= abs(del) 
acurv= abs(gdn -gup) 
if (acurv.ge.adel)  then 
  fc= ct 
  return 
 else 
  ref= up + (ct -up)*courinv 
  if (del.gt.0) then 
    fc= max(ct,fc) 
    ref= min(ref,dn) 
    fc= min(fc,ref) 
  else 
    fc= min(ct,fc) 
    ref= max(ref,dn) 
    fc= max(ref,fc) 
  endif 
end if 
return 
END                                           
