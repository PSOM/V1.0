subroutine advection_and_mixing(m,n,dtimel,step) 
  !     ---------------------------------------------                     
  USE header
  !     convecn(m,n,dtime) means use level m and write s,T to level n     
  !     computes Cx,Cy,Cz at the cell centers.                            
  !     uflx,vflx,wflx,sflx,Tflx are the fluxes defined at                
  !     cell faces (using QUICK).                                         
  !                                                                       
  integer i,j,k,n,m,step
  !     ntr is the number of tracers in variable T (Defined in header.f)  
  REAL(kind=rc_kind) :: Eighth,Half,absvel,left,right,ctr,upos,uneg 
  REAL(kind=rc_kind) :: dtimel,uflx(0:NI,0:NJ,0:NK),                      &
       vflx(0:NI,0:NJ,0:NK),wflx(0:NI,0:NJ,0:NK),                   &
       sflx(0:NI,0:NJ,0:NK),Trflx(ntr,0:NI,0:NJ,0:NK),              &
       Tflx(0:NI,0:NJ,0:NK),uTx(0:NI+1,0:NJ+1,0:NK+1),              &
       usx(0:NI+1,0:NJ+1,0:NK+1),uTrx(ntr,0:NI+1,0:NJ+1,0:NK+1),    &
       uux(0:NI+1,0:NJ+1,0:NK+1),uvx(0:NI+1,0:NJ+1,0:NK+1),         &
       uwx(0:NI+1,0:NJ+1,0:NK+1),Tdif(NI,NJ,NK),                    &
       sdif(NI,NJ,NK),Trdif(ntr,NI,NJ,NK),udif(NI,NJ,NK),           &
       vdif(NI,NJ,NK),wdif(NI,NJ,NK),dtJ                            
  REAL(kind=rc_kind) :: dTrx(ntr,0:NI+1),dsx(0:NI+1),dux(0:NI+1),        &
       dTx(0:NI+1),dvx(0:NI+1),dwx(0:NI+1),dTry(ntr,0:NJ),          &
       dTy(0:NJ),dsy(0:NJ),duy(0:NJ),dvy(0:NJ),dwy(0:NJ),           &
       dTrz(ntr,0:NK+1),dsz(0:NK+1),dTz(0:NK+1),                    &
       duz(0:NK+1),dvz(0:NK+1),dwz(0:NK+1)                          
  REAL(kind=rc_kind) :: s_restore(0:NJ+1),alph_restore(0:NJ+1) 
  REAL(kind=rc_kind) :: courinv,fc 
                                                                        
  parameter (Eighth=0.125d0, Half=0.5, courinv=2.d0) 
                                                                        
      do k=1,NK 
         do j=1,NJ 
!     for periodic-ew boundaries                                        
!     dTx,dsx,.... at faces...
                                          
           !---------------
           ! i=0

            i=0 
            do it=1,ntr 
               dTrx(it,i)= Tr(it,1,j,k,m) - Tr(it,NI,j,k,m) 
            end do 
            dTx(i)= T(1,j,k,m) - T(NI,j,k,m) 
            dsx(i)= s(1,j,k,m) - s(NI,j,k,m) 
            dux(i)= u(1,j,k,m) - u(NI,j,k,m) 
            dvx(i)= v(1,j,k,m) - v(NI,j,k,m) 
            dwx(i)= w(1,j,k,m) - w(NI,j,k,m) 

           !---------------
           ! i=1,NI-1
           
            do i=1,NI-1 
               do it=1,ntr 
                  dTrx(it,i)= Tr(it,i+1,j,k,m) - Tr(it,i,j,k,m) 
               end do 
               dsx(i)= s(i+1,j,k,m) - s(i,j,k,m) 
               dTx(i)= T(i+1,j,k,m) - T(i,j,k,m) 
               dux(i)= u(i+1,j,k,m) - u(i,j,k,m) 
               dvx(i)= v(i+1,j,k,m) - v(i,j,k,m) 
               dwx(i)= w(i+1,j,k,m) - w(i,j,k,m) 
            end do 

           !---------------
           ! i=NI

!     specially for periodic-ew boundaries dsx,dTx... from 0:NI+1       

            dTrx(:,NI)  =dTrx(:,0); dsx(NI)  =dsx(0); dTx(NI)  =dTx(0); dux(NI)  =dux(0); dvx(NI)  =dvx(0); dwx(NI)  =dwx(0);
            dTrx(:,NI+1)=dTrx(:,1); dsx(NI+1)=dsx(1); dTx(NI+1)=dTx(1); dux(NI+1)=dux(1); dvx(NI+1)=dvx(1); dwx(NI+1)=dwx(1);
                                                                        
!     periodic-ew boundaries  : do i=1,NI  (instead of i=1,NI-1)        
            do it=1,ntr 
               Tr(it,NI+1,j,k,m)= Tr(it,1,j,k,m) 
            end do 
            s(NI+1,j,k,m)= s(1,j,k,m) 
            T(NI+1,j,k,m)= T(1,j,k,m) 
            u(NI+1,j,k,m)= u(1,j,k,m) 
            v(NI+1,j,k,m)= v(1,j,k,m) 
            w(NI+1,j,k,m)= w(1,j,k,m) 
                                                    
            do i=1,NI 
               goto 108 
               if (uf(i,j,k).gt.0.d0) then 
                                                                        
                  do it=1,ntr 
                     left= Eighth*(dTrx(it,i) -dTrx(it,i-1)) 
                     ctr= Half* (Tr(it,i+1,j,k,m) +Tr(it,i,j,k,m)) 
                     fc= ctr -left 
                     call ultim(Tr(it,i-1,j,k,m),Tr(it,i,j,k,m),        &
     &                    Tr(it,i+1,j,k,m),dTrx(it,i),dTrx(it,i-1),     &
     &                    courinv,fc)                                   
                     Trflx(it,i,j,k)= uf(i,j,k)*fc 
                  end do 
                                                                        
                  left= Eighth*(dsx(i) -dsx(i-1)) 
                  ctr= Half* (s(i+1,j,k,m) +s(i,j,k,m)) 
                  fc= ctr -left 
                  call ultim(s(i-1,j,k,m),s(i,j,k,m),                   &
     &                 s(i+1,j,k,m),dsx(i),dsx(i-1),courinv,fc)         
                  sflx(i,j,k)= uf(i,j,k)*fc 

                  left= Eighth*(dTx(i) -dTx(i-1)) 
                  ctr= Half* (T(i+1,j,k,m) +T(i,j,k,m)) 
                  fc= ctr -left 
                  call ultim(T(i-1,j,k,m),T(i,j,k,m),                   &
     &                 T(i+1,j,k,m),dTx(i),dTx(i-1),courinv,fc)         
                  Tflx(i,j,k)= uf(i,j,k)*fc 
                                                                        
                  left= Eighth*(dux(i) -dux(i-1)) 
                  ctr= Half* (u(i+1,j,k,m) +u(i,j,k,m)) 
                  fc= ctr -left 
                  call ultim(u(i-1,j,k,m),u(i,j,k,m),                   &
     &                 u(i+1,j,k,m),dux(i),dux(i-1),courinv,fc)         
                  uflx(i,j,k)= uf(i,j,k)*fc 
                                                                        
                  left= Eighth*(dvx(i) -dvx(i-1)) 
                  ctr= Half* (v(i+1,j,k,m) +v(i,j,k,m)) 
                  fc= ctr -left 
                  call ultim(v(i-1,j,k,m),v(i,j,k,m),                   &
     &                 v(i+1,j,k,m),dvx(i),dvx(i-1),courinv,fc)         
                  vflx(i,j,k)= uf(i,j,k)*fc 
                                                                        
                  left= Eighth*(dwx(i) -dwx(i-1)) 
                  ctr= Half* (w(i+1,j,k,m) +w(i,j,k,m)) 
                  fc= ctr -left 
                  call ultim(w(i-1,j,k,m),w(i,j,k,m),                   &
     &                 w(i+1,j,k,m),dwx(i),dwx(i-1),courinv,fc)         
                  wflx(i,j,k)= uf(i,j,k)*fc 
                                                                        
                                                                        
               else 
                                                                        
                  do it=1,ntr 
                     right= Eighth*(dTrx(it,i+1) -dTrx(it,i)) 
                     ctr= Half* (Tr(it,i+1,j,k,m) +Tr(it,i,j,k,m)) 
                     fc= ctr - right 
                     call ultim(Tr(it,i+2,j,k,m),Tr(it,i+1,j,k,m),      &
     &                    Tr(it,i,j,k,m),dTrx(it,i),dTrx(it,i+1),       &
     &                    courinv,fc)                                   
                     Trflx(it,i,j,k)= uf(i,j,k)*fc 
                  end do 
                                                                        
                  right= Eighth*(dsx(i+1) -dsx(i)) 
                  ctr= Half* (s(i+1,j,k,m) +s(i,j,k,m)) 
                  fc= ctr - right 
                  call ultim(s(i+2,j,k,m),s(i+1,j,k,m),                 &
     &                 s(i,j,k,m),dsx(i),dsx(i+1),courinv,fc)           
                  sflx(i,j,k)= uf(i,j,k)*fc 

                  right= Eighth*(dTx(i+1) -dTx(i)) 
                  ctr= Half* (T(i+1,j,k,m) +T(i,j,k,m)) 
                  fc= ctr - right 
                  call ultim(T(i+2,j,k,m),T(i+1,j,k,m),                 &
     &                 T(i,j,k,m),dTx(i),dTx(i+1),courinv,fc)           
                  Tflx(i,j,k)= uf(i,j,k)*fc 
                                                                        
                  right= Eighth*(dux(i+1) -dux(i)) 
                  ctr= Half* (u(i+1,j,k,m) +u(i,j,k,m)) 
                  fc= ctr - right 
                  call ultim(u(i+2,j,k,m),u(i+1,j,k,m),                 &
     &                 u(i,j,k,m),dux(i),dux(i+1),courinv,fc)           
                  uflx(i,j,k)= uf(i,j,k)*fc 
                                                                        
                  right= Eighth*(dvx(i+1) -dvx(i)) 
                  ctr= Half* (v(i+1,j,k,m) +v(i,j,k,m)) 
                  fc= ctr - right 
                  call ultim(v(i+2,j,k,m),v(i+1,j,k,m),                 &
     &                 v(i,j,k,m),dvx(i),dvx(i+1),courinv,fc)           
                  vflx(i,j,k)= uf(i,j,k)*fc 
                                                                        
                  right= Eighth*(dwx(i+1) -dwx(i)) 
                  ctr= Half* (w(i+1,j,k,m) +w(i,j,k,m)) 
                  fc= ctr - right 
                  call ultim(w(i+2,j,k,m),w(i+1,j,k,m),                 &
     &                 w(i,j,k,m),dwx(i),dwx(i+1),courinv,fc)           
                  wflx(i,j,k)= uf(i,j,k)*fc 
               end if 
                                                                        
!     THIS PART was without the limiter ultim in x-dir                  
!               goto 109                                              
  108          continue 
               absvel= abs(uf(i,j,k)) 
               upos= Half*(uf(i,j,k) +absvel) 
               uneg= Half*(uf(i,j,k) -absvel) 
                                                                        
               do it=1,ntr 
               left= Eighth*(dTrx(it,i) -dTrx(it,i-1)) 
               right= Eighth*(dTrx(it,i+1) -dTrx(it,i)) 
               ctr= Half* (Tr(it,i+1,j,k,m) +Tr(it,i,j,k,m)) 
               Trflx(it,i,j,k) = uf(i,j,k)*ctr - upos*left -uneg*right 
               end do 
                                                                        
               left= Eighth*(dsx(i) -dsx(i-1)) 
               right= Eighth*(dsx(i+1) -dsx(i)) 
               ctr= Half* (s(i+1,j,k,m) +s(i,j,k,m)) 
               sflx(i,j,k) = uf(i,j,k)*ctr - upos*left -uneg*right 

               left= Eighth*(dTx(i) -dTx(i-1)) 
               right= Eighth*(dTx(i+1) -dTx(i)) 
               ctr= Half* (T(i+1,j,k,m) +T(i,j,k,m)) 
               Tflx(i,j,k) = uf(i,j,k)*ctr - upos*left -uneg*right 
                                                                        
               left= Eighth*(dux(i) -dux(i-1)) 
               right= Eighth*(dux(i+1) -dux(i)) 
               ctr= Half* (u(i+1,j,k,m) +u(i,j,k,m)) 
               uflx(i,j,k) = uf(i,j,k)*ctr - upos*left -uneg*right 
                                                                        
               left= Eighth*(dvx(i) -dvx(i-1)) 
               right= Eighth*(dvx(i+1) -dvx(i)) 
               ctr= Half* (v(i+1,j,k,m) +v(i,j,k,m)) 
               vflx(i,j,k) = uf(i,j,k)*ctr - upos*left -uneg*right 
                                                                        
               left= Eighth*(dwx(i) -dwx(i-1)) 
               right= Eighth*(dwx(i+1) -dwx(i)) 
               ctr= Half* (w(i+1,j,k,m) +w(i,j,k,m)) 
               wflx(i,j,k) = uf(i,j,k)*ctr - upos*left -uneg*right 
                                                                        
  109       end do 
!     periodic-ew boundaries                                            
            do it=1,ntr 
               Trflx(it,0,j,k)= Trflx(it,NI,j,k) 
            end do 
            sflx(0,j,k)= sflx(NI,j,k) 
            Tflx(0,j,k)= Tflx(NI,j,k) 
            uflx(0,j,k)= uflx(NI,j,k) 
            vflx(0,j,k)= vflx(NI,j,k) 
            wflx(0,j,k)= wflx(NI,j,k) 
         end do 
      end do 
                                                                        
      do k=1,NK 
         do j=1,NJ 
            do i=1,NI 
               uux(i,j,k)=  uflx(i,j,k) -uflx(i-1,j,k) 
               uvx(i,j,k)=  vflx(i,j,k) -vflx(i-1,j,k) 
               uwx(i,j,k)=  wflx(i,j,k) -wflx(i-1,j,k) 
               usx(i,j,k)=  sflx(i,j,k) -sflx(i-1,j,k) 
               uTx(i,j,k)=  Tflx(i,j,k) -Tflx(i-1,j,k) 
               do it=1,ntr 
                  uTrx(it,i,j,k)=  Trflx(it,i,j,k) -Trflx(it,i-1,j,k) 
               end do 
            end do 
         end do 
      end do 
                                                                        
!     y-direction                                                       
!     -----------                                                       
      do k=1,NK 
         do i=1,NI 
!     dTy,dsy,.... at faces...                                          
            do j=0,NJ,NJ 
               do it=1,ntr 
                  dTry(it,j)= 0.d0 
               end do 
               dsy(j)= 0.d0 
               dTy(j)= 0.d0 
               duy(j)= 0.d0 
               dvy(j)= 0.d0 
               dwy(j)= 0.d0 
            end do 
            do j=1,NJ-1 
               do it=1,ntr 
                  dTry(it,j)= Tr(it,i,j+1,k,m) - Tr(it,i,j,k,m) 
               end do 
               dsy(j)= s(i,j+1,k,m) - s(i,j,k,m) 
               dTy(j)= T(i,j+1,k,m) - T(i,j,k,m) 
               duy(j)= u(i,j+1,k,m) - u(i,j,k,m) 
               dvy(j)= v(i,j+1,k,m) - v(i,j,k,m) 
               dwy(j)= w(i,j+1,k,m) - w(i,j,k,m) 
            end do 
            s(i,0,k,m)= s(i,1,k,m) 
            T(i,0,k,m)= T(i,1,k,m) 
            do it=1,ntr 
               Tr(it,i,0,k,m)= Tr(it,i,1,k,m) 
            end do 
            u(i,0,k,m)= u(i,1,k,m) 
            v(i,0,k,m)= v(i,1,k,m) 
            w(i,0,k,m)= w(i,1,k,m) 
                                                                        
            s(i,NJ+1,k,m)= s(i,NJ,k,m) 
            T(i,NJ+1,k,m)= T(i,NJ,k,m) 
            do it=1,ntr 
               Tr(it,i,NJ+1,k,m)= Tr(it,i,NJ,k,m) 
            end do 
            u(i,NJ+1,k,m)= u(i,NJ,k,m) 
            v(i,NJ+1,k,m)= v(i,NJ,k,m) 
            w(i,NJ+1,k,m)= w(i,NJ,k,m) 
                                                                        
            do j=1,NJ-1 
                                                                        
               if (vf(i,j,k).ge.0.d0) then 
                                                                        
                  do it=1,ntr 
                     left= Eighth*(dTry(it,j) -dTry(it,j-1)) 
                     ctr= Half* (Tr(it,i,j+1,k,m) +Tr(it,i,j,k,m)) 
                     fc= ctr -left 
                     call ultim(Tr(it,i,j-1,k,m),Tr(it,i,j,k,m),        &
     &                    Tr(it,i,j+1,k,m),dTry(it,j),dTry(it,j-1),     &
     &                    courinv,fc)                                   
                     Trflx(it,i,j,k)= vf(i,j,k)*fc 
                  end do 
                                                                        
                  left= Eighth*(dsy(j) -dsy(j-1)) 
                  ctr= Half* (s(i,j+1,k,m) +s(i,j,k,m)) 
                  fc= ctr -left 
                  call ultim(s(i,j-1,k,m),s(i,j,k,m),                   &
     &                 s(i,j+1,k,m),dsy(j),dsy(j-1),courinv,fc)         
                  sflx(i,j,k)= vf(i,j,k)*fc 

                  left= Eighth*(dTy(j) -dTy(j-1)) 
                  ctr= Half* (T(i,j+1,k,m) +T(i,j,k,m)) 
                  fc= ctr -left 
                  call ultim(T(i,j-1,k,m),T(i,j,k,m),                   &
     &                 T(i,j+1,k,m),dTy(j),dTy(j-1),courinv,fc)         
                  Tflx(i,j,k)= vf(i,j,k)*fc 
                                                                        
                  left= Eighth*(duy(j) -duy(j-1)) 
                  ctr= Half* (u(i,j+1,k,m) +u(i,j,k,m)) 
                  fc= ctr -left 
                  call ultim(u(i,j-1,k,m),u(i,j,k,m),                   &
     &                 u(i,j+1,k,m),duy(j),duy(j-1),courinv,fc)         
                  uflx(i,j,k)= vf(i,j,k)*fc 
                                                                        
                  left= Eighth*(dvy(j) -dvy(j-1)) 
                  ctr= Half* (v(i,j+1,k,m) +v(i,j,k,m)) 
                  fc= ctr -left 
                  call ultim(v(i,j-1,k,m),v(i,j,k,m),                   &
     &                 v(i,j+1,k,m),dvy(j),dvy(j-1),courinv,fc)         
                  vflx(i,j,k)= vf(i,j,k)*fc 
                                                                        
                  left= Eighth*(dwy(j) -dwy(j-1)) 
                  ctr= Half* (w(i,j+1,k,m) +w(i,j,k,m)) 
                  fc= ctr -left 
                  call ultim(w(i,j-1,k,m),w(i,j,k,m),                   &
     &                 w(i,j+1,k,m),dwy(j),dwy(j-1),courinv,fc)         
                  wflx(i,j,k)= vf(i,j,k)*fc 
                                                                        
               else 
                                                                        
                  do it=1,ntr 
                     right= Eighth*(dTry(it,j+1) -dTry(it,j)) 
                     ctr= Half* (Tr(it,i,j+1,k,m) +Tr(it,i,j,k,m)) 
                     fc= ctr - right 
                     call ultim(Tr(it,i,j+2,k,m),Tr(it,i,j+1,k,m),        &
     &                    Tr(it,i,j,k,m),dTry(it,j),dTry(it,j+1),          &
     &                    courinv,fc)                                   
                     Trflx(it,i,j,k)= vf(i,j,k)*fc 
                  end do 
                                                                        
                  right= Eighth*(dsy(j+1) -dsy(j)) 
                  ctr= Half* (s(i,j+1,k,m) +s(i,j,k,m)) 
                  fc= ctr - right 
                  call ultim(s(i,j+2,k,m),s(i,j+1,k,m),                 &
     &                 s(i,j,k,m),dsy(j),dsy(j+1),courinv,fc)           
                  sflx(i,j,k)= vf(i,j,k)*fc 
                                                                        
                  right= Eighth*(dTy(j+1) -dTy(j)) 
                  ctr= Half* (T(i,j+1,k,m) +T(i,j,k,m)) 
                  fc= ctr - right 
                  call ultim(T(i,j+2,k,m),T(i,j+1,k,m),                 &
     &                 T(i,j,k,m),dTy(j),dTy(j+1),courinv,fc)           
                  Tflx(i,j,k)= vf(i,j,k)*fc 
                                                                        
                  right= Eighth*(duy(j+1) -duy(j)) 
                  ctr= Half* (u(i,j+1,k,m) +u(i,j,k,m)) 
                  fc= ctr - right 
                  call ultim(u(i,j+2,k,m),u(i,j+1,k,m),                 &
     &                 u(i,j,k,m),duy(j),duy(j+1),courinv,fc)           
                  uflx(i,j,k)= vf(i,j,k)*fc 
                                                                        
                  right= Eighth*(dvy(j+1) -dvy(j)) 
                  ctr= Half* (v(i,j+1,k,m) +v(i,j,k,m)) 
                  fc= ctr - right 
                  call ultim(v(i,j+2,k,m),v(i,j+1,k,m),                 &
     &                 v(i,j,k,m),dvy(j),dvy(j+1),courinv,fc)           
                  vflx(i,j,k)= vf(i,j,k)*fc 
                                                                        
                  right= Eighth*(dwy(j+1) -dwy(j)) 
                  ctr= Half* (w(i,j+1,k,m) +w(i,j,k,m)) 
                  fc= ctr - right 
                  call ultim(w(i,j+2,k,m),w(i,j+1,k,m),                 &
     &                 w(i,j,k,m),dwy(j),dwy(j+1),courinv,fc)           
                  wflx(i,j,k)= vf(i,j,k)*fc 
                                                                        
               end if 
                                                                        
!               absvel= dabs(vf(i,j,k))                                 
!               upos= Half*(vf(i,j,k) +absvel)                          
!               uneg= Half*(vf(i,j,k) -absvel)                          
!                                                                       
!               left= Eighth*(dTy(j) -dTy(j-1))                         
!               right= Eighth*(dTy(j+1) -dTy(j))                        
!               ctr= Half* (T(i,j+1,k,m) +T(i,j,k,m))                   
!               Tflx(i,j,k) = vf(i,j,k)*ctr - upos*left -uneg*right     
!                                                                       
!               left= Eighth*(dsy(j) -dsy(j-1))                         
!               right= Eighth*(dsy(j+1) -dsy(j))                        
!               ctr= Half* (s(i,j+1,k,m) +s(i,j,k,m))                   
!               sflx(i,j,k) = vf(i,j,k)*ctr - upos*left -uneg*right     
!                                                                       
!               left= Eighth*(duy(j) -duy(j-1))                         
!               right= Eighth*(duy(j+1) -duy(j))                        
!               ctr= Half* (u(i,j+1,k,m) +u(i,j,k,m))                   
!               uflx(i,j,k) = vf(i,j,k)*ctr - upos*left -uneg*right     
!                                                                       
!               left= Eighth*(dvy(j) -dvy(j-1))                         
!               right= Eighth*(dvy(j+1) -dvy(j))                        
!               ctr= Half* (v(i,j+1,k,m) +v(i,j,k,m))                   
!               vflx(i,j,k) = vf(i,j,k)*ctr - upos*left -uneg*right     
!                                                                       
!               left= Eighth*(dwy(j) -dwy(j-1))                         
!               right= Eighth*(dwy(j+1) -dwy(j))                        
!               ctr= Half* (w(i,j+1,k,m) +w(i,j,k,m))                   
!               wflx(i,j,k) = vf(i,j,k)*ctr - upos*left -uneg*right     
                                                                        
            end do 
!     Solid Boundaries                                                  
            do j=0,NJ,NJ 
               do it=1,ntr 
                  Trflx(it,i,j,k)= 0.d0 
               end do 
               sflx(i,j,k)= 0.d0 
               Tflx(i,j,k)= 0.d0 
               uflx(i,j,k)= 0.d0 
               vflx(i,j,k)= 0.d0 
               wflx(i,j,k)= 0.d0 
            end do 
         end do 
      end do 
                                                                        
      do k=1,NK 
         do i=1,NI 
            do j=1,NJ 
               uux(i,j,k)= (uflx(i,j,k) -uflx(i,j-1,k) ) + uux(i,j,k) 
               uvx(i,j,k)= (vflx(i,j,k) -vflx(i,j-1,k) ) + uvx(i,j,k) 
               uwx(i,j,k)= (wflx(i,j,k) -wflx(i,j-1,k) ) + uwx(i,j,k) 
               usx(i,j,k)= (sflx(i,j,k) -sflx(i,j-1,k) ) + usx(i,j,k) 
               uTx(i,j,k)= (Tflx(i,j,k) -Tflx(i,j-1,k) ) + uTx(i,j,k) 
               do it=1,ntr 
                  uTrx(it,i,j,k)= (Trflx(it,i,j,k) -Trflx(it,i,j-1,k) )    &
     &                 + uTrx(it,i,j,k)                                  
               end do 
            end do 
         end do 
      end do 
                                                                        
!     z-direction                                                       
!     -----------                                                       
      do j=1,NJ 
         do i=1,NI 
!     dTz,dsz,.... at faces...                                          
!            do k=0,NK,NK                                               
            k=0 
            do it=1,ntr 
               dTrz(it,k)= 0.d0 
            end do 
            dsz(k)= 0.d0 
            dTz(k)= 0.d0 
            duz(k)= 0.d0 
            dvz(k)= 0.d0 
            dwz(k)= 0.d0 
!                                                                       
            s(i,j,k,m)= s(i,j,k+1,m) 
            T(i,j,k,m)= T(i,j,k+1,m) 
            do it=1,ntr 
               Tr(it,i,j,k,m)= Tr(it,i,j,k+1,m) 
            end do 
            u(i,j,k,m)= u(i,j,k+1,m) 
            v(i,j,k,m)= v(i,j,k+1,m) 
            w(i,j,k,m)= w(i,j,k+1,m) 
!                                                                       
            do k=1,NK-1 
               do it=1,ntr 
                  dTrz(it,k)= Tr(it,i,j,k+1,m) - Tr(it,i,j,k,m) 
               end do 
               dsz(k)= s(i,j,k+1,m) - s(i,j,k,m) 
               dTz(k)= T(i,j,k+1,m) - T(i,j,k,m) 
               duz(k)= u(i,j,k+1,m) - u(i,j,k,m) 
               dvz(k)= v(i,j,k+1,m) - v(i,j,k,m) 
               dwz(k)= w(i,j,k+1,m) - w(i,j,k,m) 
            end do 
!-c     Top boundary - linear extrapolation, constant gradient          
!-            do k=NK,NK+1                                              
!-               dTz(k)= dtz(k-1)                                       
!-               dsz(k)= dsz(k-1)                                       
!-               duz(k)= duz(k-1)                                       
!-               dvz(k)= dvz(k-1)                                       
!-               dwz(k)= dwz(k-1)                                       
!-            end do                                                    
!-c     Linear extrapolation                                            
!-            s(i,j,NK+1,m)= 2.d0*s(i,j,NK,m) - s(i,j,NK-1,m)           
!-            T(i,j,NK+1,m)= 2.d0*T(i,j,NK,m) - T(i,j,NK-1,m)           
!-            u(i,j,NK+1,m)= 2.d0*u(i,j,NK,m) - u(i,j,NK-1,m)           
!-            v(i,j,NK+1,m)= 2.d0*v(i,j,NK,m) - v(i,j,NK-1,m)           
!-            w(i,j,NK+1,m)= 2.d0*w(i,j,NK,m) - w(i,j,NK-1,m)           
!                                                                       
!     Top boundary - Zero gradient                                      
            do k=NK,NK+1 
               do it=1,ntr 
                  dTrz(it,k)= 0.d0 
               end do 
               dsz(k)= 0.d0 
               dTz(k)= 0.d0 
               duz(k)= 0.d0 
               dvz(k)= 0.d0 
               dwz(k)= 0.d0 
            end do 
!     Linear extrapolation                                              
            s(i,j,NK+1,m)= s(i,j,NK,m) 
            T(i,j,NK+1,m)= T(i,j,NK,m) 
            do it=1,ntr 
               Tr(it,i,j,NK+1,m)= Tr(it,i,j,NK,m) 
            end do 
            u(i,j,NK+1,m)= u(i,j,NK,m) 
            v(i,j,NK+1,m)= v(i,j,NK,m) 
            w(i,j,NK+1,m)= w(i,j,NK,m) 
                                                                        
            do k=1,NK 
                                                                        
               if (wf(i,j,k).ge.0.d0) then 
                                                                        
                  do it=1,ntr 
                     left= Eighth*(dTrz(it,k) -dTrz(it,k-1)) 
                     ctr= Half* (Tr(it,i,j,k+1,m) +Tr(it,i,j,k,m)) 
                     fc= ctr -left 
                     call ultim(Tr(it,i,j,k-1,m),Tr(it,i,j,k,m),          &
     &                    Tr(it,i,j,k+1,m),dTrz(it,k),dTrz(it,k-1),        &
     &                    courinv,fc)                                   
                     Trflx(it,i,j,k)= wf(i,j,k)*fc 
                  end do 
                                                                        
                  left= Eighth*(dsz(k) -dsz(k-1)) 
                  ctr= Half* (s(i,j,k+1,m) +s(i,j,k,m)) 
                  fc= ctr -left 
                  call ultim(s(i,j,k-1,m),s(i,j,k,m),                   &
     &                 s(i,j,k+1,m),dsz(k),dsz(k-1),courinv,fc)         
                  sflx(i,j,k)= wf(i,j,k)*fc 
                                                                        
                  left= Eighth*(dTz(k) -dTz(k-1)) 
                  ctr= Half* (T(i,j,k+1,m) +T(i,j,k,m)) 
                  fc= ctr -left 
                  call ultim(T(i,j,k-1,m),T(i,j,k,m),                   &
     &                 T(i,j,k+1,m),dTz(k),dTz(k-1),courinv,fc)         
                  Tflx(i,j,k)= wf(i,j,k)*fc 
                                                                        
                  left= Eighth*(duz(k) -duz(k-1)) 
                  ctr= Half* (u(i,j,k+1,m) +u(i,j,k,m)) 
                  fc= ctr -left 
                  call ultim(u(i,j,k-1,m),u(i,j,k,m),                   &
     &                 u(i,j,k+1,m),duz(k),duz(k-1),courinv,fc)         
                  uflx(i,j,k)= wf(i,j,k)*fc 
                                                                        
                  left= Eighth*(dvz(k) -dvz(k-1)) 
                  ctr= Half* (v(i,j,k+1,m) +v(i,j,k,m)) 
                  fc= ctr -left 
                  call ultim(v(i,j,k-1,m),v(i,j,k,m),                   &
     &                 v(i,j,k+1,m),dvz(k),dvz(k-1),courinv,fc)         
                  vflx(i,j,k)= wf(i,j,k)*fc 
                                                                        
                  left= Eighth*(dwz(k) -dwz(k-1)) 
                  ctr= Half* (w(i,j,k+1,m) +w(i,j,k,m)) 
                  fc= ctr -left 
                  call ultim(w(i,j,k-1,m),w(i,j,k,m),                   &
     &                 w(i,j,k+1,m),dwz(k),dwz(k-1),courinv,fc)         
                  wflx(i,j,k)= wf(i,j,k)*fc 
                                                                        
               else 
                                                                        
                  if (k.eq.NK) then 
                     do it=1,ntr 
                        Trflx(it,i,j,k)= wf(i,j,k)*Half*                 &
     &                    (Tr(it,i,j,k+1,m) +Tr(it,i,j,k,m))              
                     end do 
                     sflx(i,j,k)= wf(i,j,k)*Half*                       &
     &                    (s(i,j,k+1,m) +s(i,j,k,m))                    
                     Tflx(i,j,k)= wf(i,j,k)*Half*                       &
     &                    (T(i,j,k+1,m) +T(i,j,k,m))                    
                     uflx(i,j,k)= wf(i,j,k)*Half*                       &
     &                    (u(i,j,k+1,m) +u(i,j,k,m))                    
                     vflx(i,j,k)= wf(i,j,k)*Half*                       &
     &                    (v(i,j,k+1,m) +v(i,j,k,m))                    
                     wflx(i,j,k)= wf(i,j,k)*Half*                       &
     &                    (w(i,j,k+1,m) +w(i,j,k,m))                    
                     goto 209 
                  endif 
                                                                        
                  do it=1,ntr 
                     right= Eighth*(dTrz(it,k+1) -dTrz(it,k)) 
                     ctr= Half* (Tr(it,i,j,k+1,m) +Tr(it,i,j,k,m)) 
                     fc= ctr - right 
                     call ultim(Tr(it,i,j,k+2,m),Tr(it,i,j,k+1,m),        &
     &                    Tr(it,i,j,k,m),dTrz(it,k),dTrz(it,k+1),          &
     &                    courinv,fc)                                   
                     Trflx(it,i,j,k)= wf(i,j,k)*fc 
                  end do 
                                                                        
                  right= Eighth*(dsz(k+1) -dsz(k)) 
                  ctr= Half* (s(i,j,k+1,m) +s(i,j,k,m)) 
                  fc= ctr - right 
                  call ultim(s(i,j,k+2,m),s(i,j,k+1,m),                 &
     &                 s(i,j,k,m),dsz(k),dsz(k+1),courinv,fc)           
                  sflx(i,j,k)= wf(i,j,k)*fc 

                  right= Eighth*(dTz(k+1) -dTz(k)) 
                  ctr= Half* (T(i,j,k+1,m) +T(i,j,k,m)) 
                  fc= ctr - right 
                  call ultim(T(i,j,k+2,m),T(i,j,k+1,m),                 &
     &                 T(i,j,k,m),dTz(k),dTz(k+1),courinv,fc)           
                  Tflx(i,j,k)= wf(i,j,k)*fc 
                                                                        
                  right= Eighth*(duz(k+1) -duz(k)) 
                  ctr= Half* (u(i,j,k+1,m) +u(i,j,k,m)) 
                  fc= ctr - right 
                  call ultim(u(i,j,k+2,m),u(i,j,k+1,m),                 &
     &                 u(i,j,k,m),duz(k),duz(k+1),courinv,fc)           
                  uflx(i,j,k)= wf(i,j,k)*fc 
                                                                        
                  right= Eighth*(dvz(k+1) -dvz(k)) 
                  ctr= Half* (v(i,j,k+1,m) +v(i,j,k,m)) 
                  fc= ctr - right 
                  call ultim(v(i,j,k+2,m),v(i,j,k+1,m),                 &
     &                 v(i,j,k,m),dvz(k),dvz(k+1),courinv,fc)           
                  vflx(i,j,k)= wf(i,j,k)*fc 
                                                                        
                  right= Eighth*(dwz(k+1) -dwz(k)) 
                  ctr= Half* (w(i,j,k+1,m) +w(i,j,k,m)) 
                  fc= ctr - right 
                  call ultim(w(i,j,k+2,m),w(i,j,k+1,m),                 &
     &                 w(i,j,k,m),dwz(k),dwz(k+1),courinv,fc)           
                  wflx(i,j,k)= wf(i,j,k)*fc 
                                                                        
               end if 
                                                                        
!               absvel= dabs(wf(i,j,k))                                 
!               upos= Half*(wf(i,j,k) +absvel)                          
!               uneg= Half*(wf(i,j,k) -absvel)                          
!                                                                       
!               left= Eighth*(dTz(k) -dTz(k-1))                         
!               right= Eighth*(dTz(k+1) -dTz(k))                        
!               ctr= Half* (T(i,j,k+1,m) +T(i,j,k,m))                   
!               Tflx(i,j,k) = wf(i,j,k)*ctr - upos*left -uneg*right     
!                                                                       
!               left= Eighth*(dsz(k) -dsz(k-1))                         
!               right= Eighth*(dsz(k+1) -dsz(k))                        
!               ctr= Half* (s(i,j,k+1,m) +s(i,j,k,m))                   
!               sflx(i,j,k) = wf(i,j,k)*ctr - upos*left -uneg*right     
!                                                                       
!               left= Eighth*(duz(k) -duz(k-1))                         
!               right= Eighth*(duz(k+1) -duz(k))                        
!               ctr= Half* (u(i,j,k+1,m) +u(i,j,k,m))                   
!               uflx(i,j,k) = wf(i,j,k)*ctr - upos*left -uneg*right     
!                                                                       
!               left= Eighth*(dvz(k) -dvz(k-1))                         
!               right= Eighth*(dvz(k+1) -dvz(k))                        
!               ctr= Half* (v(i,j,k+1,m) +v(i,j,k,m))                   
!               vflx(i,j,k) = wf(i,j,k)*ctr - upos*left -uneg*right     
!                                                                       
!               left= Eighth*(dwz(k) -dwz(k-1))                         
!               right= Eighth*(dwz(k+1) -dwz(k))                        
!               ctr= Half* (w(i,j,k+1,m) +w(i,j,k,m))                   
!               wflx(i,j,k) = wf(i,j,k)*ctr - upos*left -uneg*right     
                                                                        
  209       end do 
!     Solid Boundary                                                    
            k=0 
            do it=1,ntr 
               Trflx(it,i,j,k)= 0.d0 
            end do 
            sflx(i,j,k)= 0.d0 
            Tflx(i,j,k)= 0.d0 
            uflx(i,j,k)= 0.d0 
            vflx(i,j,k)= 0.d0 
            wflx(i,j,k)= 0.d0 
!     Top surface (done already)                                        
!            k=NK                                                       
!            Tflx(i,j,k)= T(i,j,k,m)*wf(i,j,k)                          
!            sflx(i,j,k)= s(i,j,k,m)*wf(i,j,k)                          
!            uflx(i,j,k)= u(i,j,k,m)*wf(i,j,k)                          
!            vflx(i,j,k)= v(i,j,k,m)*wf(i,j,k)                          
!            wflx(i,j,k)= w(i,j,k,m)*wf(i,j,k)                          
         end do 
      end do 
                                                                        
!      call diffusion(sdif,Trdif,udif,vdif,m,step) !original
       call mixing_vertical(sdif,Tdif,Trdif,udif,vdif,m,step) !modified by amala, add Tr in
                                                                        
!=      call biharmonic(udif,vdif,m,step)                               
!1/25/05 moved     call viscous(udif,vdif,wdif,m)                       
!                                                                       
      do 250 k=1,NK 
         do 251 j=1,NJ 
            do 252 i=1,NI 
               uux(i,j,k)= uux(i,j,k) + (uflx(i,j,k) -uflx(i,j,k-1) )   &
     &              - udif(i,j,k)                                       
               uvx(i,j,k)= uvx(i,j,k) + (vflx(i,j,k) -vflx(i,j,k-1) )   &
     &              - vdif(i,j,k)                                       
               uwx(i,j,k)= uwx(i,j,k) + (wflx(i,j,k) -wflx(i,j,k-1) ) 
!     &              - wdif(i,j,k)                                      
               usx(i,j,k)= usx(i,j,k) + (sflx(i,j,k) -sflx(i,j,k-1) )   &
     &              - sdif(i,j,k)                                       
               uTx(i,j,k)= uTx(i,j,k) + (Tflx(i,j,k) -Tflx(i,j,k-1) )   &
     &              - Tdif(i,j,k)                                       
!     put back sdif to calc. advec of rho, non-dim by Tl                
               rhoadv(i,j,k)= (usx(i,j,k) +sdif(i,j,k))/Jac(i,j,k) 
               do it=1,ntr 
                  uTrx(it,i,j,k)= uTrx(it,i,j,k) + (Trflx(it,i,j,k)        &
     &                 -Trflx(it,i,j,k-1) )                              
!     NO VERTICAL DIFFUSION of PASSIVE TRACERS                          
!=     &                 - Trdif(i,j,k)                                  
               end do 
  252       continue 
  251    continue 
  250 continue 
      call mixing_horizontal(udif,vdif,wdif,sdif,Tdif,Trdif,m) 
      do k=1,NK
         do j=1,NJ 
            do i=1,NI 
               uux(i,j,k)= uux(i,j,k)-udif(i,j,k) 
               uvx(i,j,k)= uvx(i,j,k)-vdif(i,j,k) 
               uwx(i,j,k)= uwx(i,j,k)-wdif(i,j,k) 
               usx(i,j,k)= usx(i,j,k)-sdif(i,j,k) 
               uTx(i,j,k)= uTx(i,j,k)-Tdif(i,j,k) 
!     LATERAL DIFFUSION OF PASSIVE TRACER IS TURNED ON / OFF            
!=               do it=1,ntr                                            
!=                  uTrx(it,i,j,k)= uTrx(it,i,j,k)-Trdif(it,i,j,k)         
!=               end do                                                 
            end do 
         end do 
      end do 
!                                                                       
!                                                                       
!     *$*  ASSERT CONCURRENT CALL                                       
      do 300 j=1,NJ 
         do 301 i=1,NI 
            do 302 k=1,NK 
               dtJ= dtimel/Jac(i,j,k) 
               cx(i,j,k)= u(i,j,k,0) -dtJ*uux(i,j,k) 
               cy(i,j,k)= v(i,j,k,0) -dtJ*uvx(i,j,k) 
               cz(i,j,k)= w(i,j,k,0) -dtJ*uwx(i,j,k) 
               s(i,j,k,n)= s(i,j,k,0) -dtJ*usx(i,j,k) 
               T(i,j,k,n)= T(i,j,k,0) -dtJ*uTx(i,j,k) 
!==     &              - alph_restore(j)*(s(i,j,k,m)-s_restore(j))      
               do it=1,ntr 
                  Tr(it,i,j,k,n)= Tr(it,i,j,k,0) -dtJ*uTrx(it,i,j,k) 
               end do 
!     T_restore= 0.d0                                                   
  302       continue 
  301    continue 
  300 continue 
!                                                                       
!+      call remineralize(n) ( called from momentum)                    
!     Surface flux                                                      
!=      call surfaceflux(dtime,n)                                       
!                                                                       

 call restore_bndry(n)

      return 
                                                                        
!     This extrapolation wasn't there before and is not                 
!     needed because rho from outside the bndry is  not used            
!     Boundary conditions                                               
      do k=1,NK 
         do i=1,NI 
            s(i,0,k,n)= 3.d0*(s(i,1,k,n)-s(i,2,k,n))+s(i,3,k,n) 
            s(i,NJ+1,k,n)= 3.d0*(s(i,NJ,k,n)-s(i,NJ-1,k,n))             &
     &           +s(i,NJ-2,k,n)                                         
            T(i,0,k,n)= 3.d0*(T(i,1,k,n)-T(i,2,k,n))+T(i,3,k,n) 
            T(i,NJ+1,k,n)= 3.d0*(T(i,NJ,k,n)-T(i,NJ-1,k,n))             &
     &           +T(i,NJ-2,k,n)                                         
            do it=1,ntr 
               Tr(it,i,0,k,n)= 3.d0*(Tr(it,i,1,k,n)-Tr(it,i,2,k,n))        &
     &              +Tr(it,i,3,k,n)                                      
               Tr(it,i,NJ+1,k,n)= 3.d0*(Tr(it,i,NJ,k,n)-Tr(it,i,NJ-1,k,n)) &
     &              +Tr(it,i,NJ-2,k,n)                                   
            end do 
         end do 
      end do 
      do j=0,NJ+1 
         do i=0,NI+1 
            s(i,j,NK+1,n)= 3.d0*(s(i,j,NK,n)-s(i,j,NK-1,n))             &
     &           +s(i,j,NK-2,n)                                         
            s(i,j,0,n)= 3.d0*(s(i,j,1,n)-s(i,j,2,n))                    &
     &           +s(i,j,3,n)                                            
            T(i,j,NK+1,n)= 3.d0*(T(i,j,NK,n)-T(i,j,NK-1,n))             &
     &           +T(i,j,NK-2,n)                                         
            T(i,j,0,n)= 3.d0*(T(i,j,1,n)-T(i,j,2,n))                    &
     &           +T(i,j,3,n)                                            
            do it=1,ntr 
               Tr(it,i,j,NK+1,n)= 3.d0*(Tr(it,i,j,NK,n)-Tr(it,i,j,NK-1,n)) &
     &              +Tr(it,i,j,NK-2,n)                                   
               Tr(it,i,j,0,n)= 3.d0*(Tr(it,i,j,1,n)-Tr(it,i,j,2,n))        &
     &              +Tr(it,i,j,3,n)                                      
            end do 
         end do 
      end do 
!     periodic-ew boundaries                                            
      do k=0,NK+1 
         do j=0,NJ+1 
            s(NI+1,j,k,n)= s(1,j,k,n) 
            s(0,j,k,n)= s(NI,j,k,n) 
            T(NI+1,j,k,n)= T(1,j,k,n) 
            T(0,j,k,n)= T(NI,j,k,n) 
            do it=1,ntr 
               Tr(it,NI+1,j,k,n)= Tr(it,1,j,k,n) 
               Tr(it,0,j,k,n)= Tr(it,NI,j,k,n) 
            end do 
         end do 
      end do 
   
      s(:,0,:,n) = s(:,1,:,n)
      T(:,0,:,n) = T(:,1,:,n)
                                                                        
!     call sTbc                                                         
!==      call sTbc_periodicew(n)                                        
!                                                                       
      return 
      END                                           
