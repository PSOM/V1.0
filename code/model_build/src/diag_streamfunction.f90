subroutine diag_streamfunction(rhobar,n2bar,psiv,psiw,by,bybar,bzbar,     &
     vbbar,wbbar,psieddy,wpvbar,pvbar,Feddydiv,Freyndiv,          &
     cybar,czbar,vcbar,wcbar)                                     
  !     --------------                                                    
  USE header
!     Solve psi_yy + psi_zz = vor in 2-D                                
!     BC psi =0 on all boundaries                                       
                                                                        
      implicit NONE 
      integer i,j,k,n,jup,jdn,kup,kdn
      REAL(kind=rc_kind) :: fcor(NJ),wzbar(NJ,NK),vybar(NJ),vbar(NJ,NK),     &
     &     wbar(NJ,NK),psiw(NJ,NK),psiv(NJ,NK),bybar(NJ,NK), &
     &     bzbar(NJ,NK),bbar(NJ,NK),zcbar(NJ,NK),ubar(NJ,NK),           &
     &     b(NI,NJ,NK),by(NJ,NK),bz(NI,NJ,NK),                       &
     &     vbbar(NJ,NK),wbbar(NJ,NK),psieddy(NJ,NK),rhobar(NJ,NK),      &
     &     wpvbar(NJ,NK),pvbar(NJ,NK),                    &
     &     pv1bar(NJ,NK),pv2bar(NJ,NK),     &
     &     pv3bar(NJ,NK),n2bar(NJ,NK)                                   
      REAL(kind=rc_kind) :: Feddydiv(NJ,NK),Freyndiv(NJ,NK),                 &
     &     cybar(ntr,NJ,NK),czbar(ntr,NJ,NK),                           &
     &     cbar(ntr,NJ,NK),                                             &
     &     vcbar(ntr,NJ,NK),wcbar(ntr,NJ,NK),csfeddy(ntr,NJ,NK)         
      REAL(kind=rc_kind) :: psivprev,psiwprev,dyc,dz,vp,wp,bp,cp,par,        &
     &     rhomean,dyinv                                                
                                                                        
      n=0 
!     WL= eps*delta*UL                                                  
      call evalrho(rho,0) 
                                                                        
      call diag_pv(n) 
                                                                        
      rhomean=0.d0 
      do k=1,NK 
         do j=1,NJ 
            wzbar(j,k)= 0.d0 
            do i=1,NI 
               wzbar(j,k)= wzbar(j,k) +wz(i,j,k) 
               rhomean= rho(i,j,k) + rhomean 
            end do 
            wzbar(j,k)= wzbar(j,k)/dble(NI) 
         end do 
      end do 
      rhomean= rhomean/dble(NI*NJ*NK) 
                                                                        
!     Calculate b  in m/s^2                                             
      do k=1,NK 
         do j=1,NJ 
            do i=1,NI 
!=               b(i,j,k)= -gpr*10.d0*(rho(i,j,k)-rhomean)/rhomean      
               b(i,j,k)= -gpr*10.d0*(rho(i,j,k)-R0)/R0 
            end do 
         end do 
      end do 
!      calculate by on section i=NI/2
      i=NI/2
      do k=1,NK
         do j=1,NJ
            jup= j+1 
            jdn=j-1 
            if (j.eq.1) jdn=j 
            if (j.eq.NJ) jup= j 
            dyc= (yc(jup)-yc(jdn))*1.d3 
            by(j,k)= (b(i,jup,k)-b(i,jdn,k))/dyc 
         end do
      end do

!     Average zonally                                                   
!     =================                                                 
      do k=1,NK 
         do j=1,NJ 
                                                                        
            ubar(j,k) = 0.d0 
            vbar(j,k) = 0.d0 
            wbar(j,k) = 0.d0 
            bbar(j,k) = 0.d0 
            rhobar(j,k)= 0.d0 
            n2bar(j,k)= 0.d0 
            zcbar(j,k)= 0.d0 
            do it=1,ntr 
               cbar(it,j,k) = 0.d0 
            end do 
                                                                        
            pv1bar(j,k)= 0.d0 
            pv2bar(j,k)= 0.d0 
            pv3bar(j,k)= 0.d0 
            pvbar(j,k)= 0.d0 
                                                                        
            do i=1,NI 
               ubar(j,k)= ubar(j,k) + u(i,j,k,n)*UL 
               vbar(j,k)= vbar(j,k) + v(i,j,k,n)*UL 
               wbar(j,k)= wbar(j,k) + w(i,j,k,n)*WL 
               bbar(j,k)= bbar(j,k) + b(i,j,k) 
               rhobar(j,k)= rhobar(j,k) + rho(i,j,k) 
               n2bar(j,k)= n2bar(j,k) + freqN2(i,j,k) 
               zcbar(j,k) = zcbar(j,k) + zc(i,j,k) 
               do it=1,ntr 
                  cbar(it,j,k)= cbar(it,j,k) + Tr(it,i,j,k,n) 
               end do 
                                                                        
               pv1bar(j,k) = pv1bar(j,k) + pv1(i,j,k) 
               pv2bar(j,k) = pv2bar(j,k) + pv2(i,j,k) 
               pv3bar(j,k) = pv3bar(j,k) + pv3(i,j,k) 
               pvbar(j,k) = pv1bar(j,k)+pv2bar(j,k)+pv3bar(j,k) 
            end do 
            ubar(j,k)= ubar(j,k)/dble(NI) 
            vbar(j,k)= vbar(j,k)/dble(NI) 
            wbar(j,k)= wbar(j,k)/dble(NI) 
            bbar(j,k)= bbar(j,k)/dble(NI) 
            rhobar(j,k)= rhobar(j,k)/dble(NI) 
            n2bar(j,k)= n2bar(j,k)/dble(NI) 
            zcbar(j,k)= DL*zcbar(j,k)/dble(NI) 
            do it=1,ntr 
               cbar(it,j,k)= cbar(it,j,k)/dble(NI) 
            end do 
                                                                        
            pv1bar(j,k)= pv1bar(j,k)/dble(NI) 
            pv2bar(j,k)= pv2bar(j,k)/dble(NI) 
            pv3bar(j,k)= pv3bar(j,k)/dble(NI) 
            pvbar(j,k)= pvbar(j,k)/dble(NI) 
         end do 
      end do 
                                                                        
      do j=1,NJ 
         vybar(j)= 0.d0 
         do i=1,NI 
            vybar(j) = vybar(j) + vy(i,j) 
         end do 
         vybar(j)= vybar(j)/dble(NI) 
         fcor(j) = FPAR*ffc(NI/2,j) 
      end do 
                                                                        
      do k=1,NK 
         do j=1,NJ 
            vbbar(j,k)= 0.d0 
            wbbar(j,k)= 0.d0 
            wpvbar(j,k)= 0.d0 
            do it=1,ntr 
               vcbar(it,j,k)= 0.d0 
               wcbar(it,j,k)= 0.d0 
            end do 
            do i=1,NI 
               vp= UL*v(i,j,k,n) - vbar(j,k) 
               wp= WL*w(i,j,k,n) - wbar(j,k) 
               bp = b(i,j,k)- bbar(j,k) 
               do it=1,ntr 
                  cp = Tr(it,i,j,k,n)- cbar(it,j,k) 
                  vcbar(it,j,k)= vcbar(it,j,k) + vp*cp 
                  wcbar(it,j,k) = wcbar(it,j,k) + wp*cp 
               end do 
                                                                        
               vbbar(j,k)= vbbar(j,k) + vp*bp 
               wbbar(j,k) = wbbar(j,k) + wp*bp 
                                                                        
               wpvbar(j,k) = wpvbar(j,k) + w(i,j,k,n)*WL*               &
     &              ( pv1(i,j,k)+pv2(i,j,k)+pv3(i,j,k) )                
            end do 
            vbbar(j,k)= vbbar(j,k)/dble(NI) 
            wbbar(j,k)= wbbar(j,k)/dble(NI) 
            wpvbar(j,k)= wpvbar(j,k)/dble(NI) 
            do it=1,ntr 
               vcbar(it,j,k)= vcbar(it,j,k)/dble(NI) 
               wcbar(it,j,k)= wcbar(it,j,k)/dble(NI) 
            end do 
         end do 
      end do 
!     Y-derivative                                                      
!     ---------------                                                   
      do k=1,NK 
         do j=1,NJ 
            jup= j+1 
            jdn=j-1 
            if (j.eq.1) jdn=j 
            if (j.eq.NJ) jup= j 
            dyc= (yc(jup)-yc(jdn))*1.d3 
            bybar(j,k)= (bbar(jup,k)-bbar(jdn,k))/dyc 
!     b_y is in per s^2                                                 
            do it=1,ntr 
               cybar(it, j,k)= (cbar(it,jup,k)-cbar(it,jdn,k))/dyc 
            end do 
!     c_y is in per m, if c is non-dim                                  
         end do 
      end do 
                                                                        
!     Z-derivative                                                      
!     -----------------                                                 
      do k=1,NK 
         kup=k+1 
         kdn=k-1 
         if (k.eq.1) kdn=k 
         if (k.eq.NK) kup=k 
         do j=1,NJ 
                                           ! already multipl by DL      
            dz= zcbar(j,kup)-zcbar(j,kdn) 
            bzbar(j,k)= (bbar(j,kup)-bbar(j,kdn))/dz 
            do it=1,ntr 
               czbar(it,j,k)= (cbar(it,j,kup)-cbar(it,j,kdn))/dz 
            end do 
         end do 
      end do 
                                                                        
      par=1.d-3 
      do k=1,NK 
         do j=1,NJ 
            psieddy(j,k)= par*(par*vbbar(j,k)*bzbar(j,k)                &
     &           - wbbar(j,k) * bybar(j,k)/par) /                       &
     &           (bybar(j,k)*bybar(j,k) +par*par*bzbar(j,k)*bzbar(j,k)) 
!     Eddy stream function based on tracer blows up where cy,cz are zero
!            do it=1,ntr                                                
!               csfeddy(it,j,k)= par*(par*vcbar(it,j,k)*czbar(it,j,k)   
!     &              - wcbar(it,j,k) *cybar(it,j,k)/par) /              
!     &              (cybar(it,j,k)*cybar(it,j,k)                       
!     &              +par*par*czbar(it,j,k)*czbar(it,j,k))              
!            end do                                                     
         end do 
      end do 
                                                                        
      n=0 
!     WL = EPS*delta*UL                                                 
      do j=1,NJ 
         psivprev= 0.d0 
         do k=1,NK 
             psiv(j,k)= psivprev -DL*vbar(j,k)/wzbar(j,k) 
             psivprev= psiv(j,k) 
         end do 
      end do 
      do k=1,NK 
         psiwprev= 0.d0 
         do j=1,NJ 
            psiw(j,k)= psiwprev + LEN*wbar(j,k)/vybar(j) 
            psiwprev= psiw(j,k) 
         end do 
      end do 
                                                                        
!     Calculate divergence of F_reynolds and F_eddy (Residual is the sum
!     ----------------------------------------------                    
      do k=1,NK 
         kup=k+1 
         kdn=k-1 
         if (k.eq.1) kdn=k 
         if (k.eq.NK) kup=k 
         do j=1,NJ 
            jup= j+1 
            jdn=j-1 
            if (j.eq.1) jdn=j 
            if (j.eq.NJ) jup= j 
            dyc= (yc(jup)-yc(jdn))*1.d3 
                                           ! already multipl by DL      
            dz= zcbar(j,kup)-zcbar(j,kdn) 
            FReyndiv(j,k)= (vbbar(jup,k)-vbbar(jdn,k))/dyc              &
     &           + (wbbar(j,kup)-wbbar(j,kdn))/dz                       
            Feddydiv(j,k)= (psieddy(jup,k)*bzbar(jup,k)-                &
     &           psieddy(jdn,k)*bzbar(jdn,k))/dyc                       &
     &           - (psieddy(j,kup)*bybar(j,kup)-                        &
     &           psieddy(j,kdn)*bybar(j,kdn))/dz                        
         end do 
      end do 
                                                                        
      return 
      END                                           
