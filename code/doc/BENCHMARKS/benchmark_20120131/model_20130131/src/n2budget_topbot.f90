subroutine N2budget(step) 
  !     =======================                                           
  USE header
  !     Use the Gauss' Thm form - i.e. fluxes at  top and bottom only.    
  !     Calculate N2  (units: per s^2) and terms in EQ 6 in Thomas-Ferrari
  !     bouyancy b in (m per s^2)                                         
  !     vor in (per second)                                               
  !     PV in (per second^3)                                              
  !     zeta*b is vor*b. So zeta*b/(f H) is in (per s^2)                  
                                                                        
      implicit none 
      integer step,i,i2,j,k,m,n,jup,jdn,iup,idn,kup,kdn 
      INTEGER :: ktest(3)
      REAL(kind=rc_kind) :: depth,dz,dzt,dzb,dyinv,dxinv,dzinv,zbt,          &
     &     flocalinv,flocal,bx,by,bz,grav,adv,vol                       
      REAL(kind=rc_kind) :: avezf(NK),avezc(NK),avexyN2(NK),testdepth(3) 
      REAL(kind=rc_kind) :: b(NI,NJ,NK),Dbdt(NI,NJ,NK) 
!     pv is in header and is common                                     
      REAL(kind=rc_kind) :: pvx,pvy,pvz,gravR0,dia,                                      &
     &     gradBxF_top,gradBxF_bot,adv_top,adv_bot,dia_top,dia_bot      
                                                                        
!     WL= eps*delta*UL                                                  
      gravR0= gpr*10.d0/R0 
      n=0 
                                                                        
!     Horizontally average                                              
      do k=1,NK 
         avexyN2(k)= 0.d0 
         avezf(k)= 0.d0 
         avezc(k)= 0.d0 
         do j=1,NJ 
            do i=1,NI 
               avexyN2(k)= freqN2(i,j,k) +avexyN2(k) 
               avezf(k)= -zf(i,j,k) +avezf(k) 
               avezc(k)= -zc(i,j,k) +avezc(k) 
            end do 
         end do 
         avexyN2(k)= avexyN2(k)/dble(NI*NJ) 
         avezf(k)= DL*avezf(k)/dble(NI*NJ) 
         avezc(k)= DL*avezc(k)/dble(NI*NJ) 
      end do 
!     avezf and avezc are positive (in m)                               
                                                                        
!     Integrate vertically, down to 0.9*MLD, 0.75*MLD , and 0.4*MLD     
      testdepth(3)= 0.9d0*mldepth 
      testdepth(2)= 0.75d0*mldepth 
      testdepth(1)= 0.5d0*mldepth 
                                                                        
      do m=1,3 
         do k=NK,1,-1 
!     avezf and avezc are positive (refer to depth of grid level)       
            if (avezf(k-1).ge.testdepth(m))  then 
               dzb= avezf(k-1)- testdepth(m) 
               dzt= testdepth(m) -avezf(k) 
               if (dzt.le.dzb) then 
                  ktest(m)= k 
               else 
                  ktest(m)= k-1 
               end if 
               goto 101 
            end if 
         end do 
  101 end do 
!      do m=1,3                                                         
!         write(6,*) 'testdep',testdepth(m),ktest(m),avezf(ktest(m)),   
!     &        avezf(ktest(m)+1), avezf(ktest(m)-1)                     
!         testdepth(m)= avezf(ktest(m))                                 
!      end do                                                           
!      write(6,*) 'stop in n2budget'                                    
!      stop                                                             
                                                                        
!     Calculate mldN2                                                   
!     ---------------                                                   
      do m=1,3 
         depth = 0.d0 
         mldN2(m)= 0.d0 
         do k=NK,1,-1 
            if (k.ge.ktest(m)) then 
                                           ! avezf is positive (depth)  
               dz = -(avezf(k)-avezf(k-1)) 
               mldN2(m)= mldN2(m) + avexyN2(k)*dz 
               depth = depth + dz 
            end if 
         end do 
  201    mldN2(m)= mldN2(m)/depth 
                                                                        
         if (step.eq.0) then 
            mldN2init(m)= mldN2(m) 
         else 
            mldN2(m)= mldN2(m)-mldN2init(m) 
         end if 
      end do 
                                                                        
!     FRONTAL TERM                                                      
!     =============================                                     
!     Calculate zeta*b top and bot divided by fH in per s^2             
      call evalrho(rho,0) 
                                                                        
      call vort(0) 
                                                                        
!     Calculate b  in m/s^2                                             
!     -----------------------                                           
      do k=1,NK 
         do j=1,NJ 
            do i=1,NI 
               b(i,j,k)= -gravR0*(rho(i,j,k)-R0) 
            end do 
         end do 
      end do 
                                                                        
      do m=1,3 
         zbbot(m)= 0.d0 
      end do 
      zbt= 0.d0 
      do j=1,NJ 
         do i=1,NI 
            flocalinv= 1.d0/(ffc(i,j)*FPAR) 
            zbt = zbt - vor(i,j,NK)*b(i,j,NK)*flocalinv 
            do m=1,3 
               zbbot(m) = zbbot(m) - vor(i,j,ktest(m))*b(i,j,ktest(m))  &
     &              *flocalinv                                          
            end do 
         end do 
      end do 
      do m=1,3 
         zbtop(m) = zbt/ dble(NI*NJ) / testdepth(m) 
         zbbot(m) = zbbot(m)/dble(NI*NJ) / testdepth(m) 
      end do 
                                                                        
      if (step.eq.0) then 
         do m=1,3 
            zbtopinit(m)= zbtop(m) 
            zbbotinit(m)=zbbot(m) 
         end do 
      else 
         do m=1,3 
            zbtop(m)= zbtop(m) - zbtopinit(m) 
            zbbot(m)= zbbot(m) - zbbotinit(m) 
         end do 
      end if 
!     =======================================================           
!     ADVEC term  u \cdot  grad (q)                                     
!     Calculate Time_integ [ volume_integ (Div ( u q )) ]   in per s^2  
!     i.e.  Time_integ ( Volume Integ (u \cdot grad(q) ))               
      call pvcalc(0) 
!     PV is in per s^3, w is in m/s, testdepth in m, f in per s         
               
                                                         
      do k=1,NK 
         do j=1,NJ 
            do i=0,NI+1
               if(i==0) then
                 i2=NK
                else if(i==NI+1) then
                 i2=1
                else
                 i2=i
               endif
               pv(i,j,k)= pv1(i2,j,k)+pv2(i2,j,k) + pv3(i2,j,k) 
            end do 
         end do 
      end do 
                                                         
      do k=1,NK 
         do j=1,NJ 
            do i=1,NI 
               Dbdt(i,j,k)= -gravR0*((rho(i,j,k)-rhoprev(i,j,k))/dtf    &
     &              + rhoadv(i,j,k) )/TL                                
               rhoprev(i,j,k)= rho(i,j,k) 
            end do 
         end do 
      end do 
!     Integrate up to ktest(m) of deepest level, i.e. m=3               
                                                                        
      adv= 0.d0 
      dia = 0.d0 
      do m=1,3 
         friction(m)= 0.d0 
         diabatic(m)= 0.d0 
         advecpv(m)= 0.d0 
         frictop(m)= 0.d0 
         diatop(m)=0.d0 
         fricbot(m)= 0.d0 
         diabot(m)=0.d0 
      end do 
                                                                        
!     Integrate upto the deepest testdepth, i.e. testdepth(1)           
      do j=1,NJ 
         jup=j+1 
         jdn=j-1 
         dyinv= 1.d0/(2.d0*dy) 
         if (j.eq.NJ) then 
            jup=j 
            dyinv= 1.d0/dy 
         else if (j.eq.1) then 
            jdn=j 
            dyinv= 1.d0/dy 
         end if 
         do i=1,NI 
                                                                        
            iup=i+1 
            idn=i-1 
            dxinv= 1.d0/(2.d0*dx) 
!     Periodic-ew boundaries                                            
            if (i.eq.NI) then 
               iup=1 
            else if (i.eq.1) then 
               idn=NI 
            end if 
                                                                        
            flocalinv= 1.d0/(ffc(i,j)*FPAR) 
!     TOP                                                               
            k=NK 
!     ADVEC TERM                                                        
            adv_top= 0.d0 
                                                                        
!     FRIC TERM                                                         
!     CALCULATE bx,by,bz                                                
            bx= dxinv*(b(iup,j,k)-b(idn,j,k)) 
            by= dyinv*(b(i,jup,k)-b(i,jdn,k)) 
                                                                        
                                                                        
!===                                                                    
!                                                                       
            gradBxF_top=  bx*fricv(i,j,k) - by*fricu(i,j,k) 
!            if (step.eq.1) then                                        
!            if (i.eq.48) write(60,*) j,bx,by                           
!            if (i.eq.48) write(61,*)                                   
!     &           j,fricv(i,j,k),fricu(i,j,k),gradBxF_top               
!            endif                                                      
!==                                                                     
!     DIA TERM                                                          
            dia_top =  - (vor3(i,j,k)+ffc(i,j)*FPAR)*dbdt(i,j,k) 
                                                                        
                                                                        
            do m=1,3 
               k=ktest(m) 
                                                                        
               adv_bot= - flocalinv *WL*w(i,j,k,0)*pv(i,j,k) 
                                                                        
!     CALCULATE bx,by,bz                                                
               bx= dxinv*(b(iup,j,k)-b(idn,j,k)) 
               by= dyinv*(b(i,jup,k)-b(i,jdn,k)) 
               gradBxF_bot=  bx*fricv(i,j,k) - by*fricu(i,j,k) 
                                                                        
               dia_bot =  - (vor3(i,j,k)+ffc(i,j)*FPAR)*dbdt(i,j,k) 
               advecpv(m)=advecpv(m)+ dtf*TL* (adv_top-adv_bot)/        &
     &              testdepth(m)                                        
               friction(m) =friction(m) +dtf*TL*flocalinv*              &
     &              (gradBxF_top-gradBxF_bot)/                          &
     &              testdepth(m)                                        
               frictop(m) = frictop(m) +dtf*TL*flocalinv*               &
     &              gradBxF_top/testdepth(m)                            
               fricbot(m) = fricbot(m) +dtf*TL*flocalinv*               &
     &              gradBxF_bot/testdepth(m)                            
               diabatic(m) = diabatic(m)+ dtf*TL*flocalinv*             &
     &              (dia_top-dia_bot)/testdepth(m)                      
               diatop(m) = diatop(m)+ dtf*TL*flocalinv*                 &
     &              dia_top/testdepth(m)                                
               diabot(m) = diabot(m)+ dtf*TL*flocalinv*                 &
     &              dia_bot/testdepth(m)                                
!     write(6,*) m,friction(m),advecpv(m), diabatic(m)                  
                                                                        
            end do 
         end do 
      end do 
                                                                        
      do m=1,3 
         advecpv(m)= advecpv(m)/dble(NI*NJ) 
         friction(m)= friction(m)/dble(NI*NJ) 
         diabatic(m)= diabatic(m)/dble(NI*NJ) 
         frictop(m)= frictop(m)/dble(NI*NJ) 
         diatop(m)= diatop(m)/dble(NI*NJ) 
         fricbot(m)= fricbot(m)/dble(NI*NJ) 
         diabot(m)= diabot(m)/dble(NI*NJ) 
      end do 
      return 
      END                                           
