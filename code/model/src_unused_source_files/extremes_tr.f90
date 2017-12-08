subroutine extremes(step) 
  !     -----------------------------------                               
  USE header
  !     Modified for tracer                                               
      integer i,j,k,step 
      REAL(kind=rc_kind) :: umax,vmax,wmax,Trmax(ntr),smax,Tmax,umin,vmin,wmin,  &
          smin,Tmin,Trmin(ntr)                                               
      REAL(kind=rc_kind) :: wup(NK),wdn(NK),const 
                                                                        
!     Compute min and max quantities                                    
!     ------------------------------                                    
      umax= -1.d16 
      vmax= -1.d16 
      wmax= -1.d16 
      smax= -1.d16 
      Tmax= -1.d16 
      do it=1,ntr 
         Trmax(it)= -1.d16 
      end do 
!                                                                       
      umin= 1.d16 
      vmin= 1.d16 
      wmin= 1.d16 
      smin= 1.d16 
      Tmin= 1.d16 
      do it=1,ntr 
         Trmin(it)= 1.d16 
      end do 
!                                                                       
      do k=1,NK 
         do j=1,NJ 
            do i=1,NI 
               if (u(i,j,k,0).gt.umax) umax=u(i,j,k,0) 
               if (v(i,j,k,0).gt.vmax) vmax=v(i,j,k,0) 
               if (w(i,j,k,0).gt.wmax) wmax=w(i,j,k,0) 
               if (s(i,j,k,0).gt.smax) smax=s(i,j,k,0) 
               if (T(i,j,k,0).gt.Tmax) Tmax=T(i,j,k,0) 
               do it=1,ntr 
                  if (Tr(it,i,j,k,0).gt.Trmax(it)) Trmax(it)=           &
                       Tr(it,i,j,k,0)                                    
               end do 
               if (u(i,j,k,0).lt.umin) umin=u(i,j,k,0) 
               if (v(i,j,k,0).lt.vmin) vmin=v(i,j,k,0) 
               if (w(i,j,k,0).lt.wmin) wmin=w(i,j,k,0) 
               if (s(i,j,k,0).lt.smin) smin=s(i,j,k,0) 
               if (T(i,j,k,0).lt.Tmin) Tmin=T(i,j,k,0) 
               do it=1,ntr 
                  if (Tr(it,i,j,k,0).lt.Trmin(it)) Trmin(it)=            &
                       Tr(it,i,j,k,0)                                    
               end do 
            end do 
         end do 
      end do 
                                                                        
      open(unit=88,file='umax.out',access='append') 
      open(unit=89,file='umin.out',access='append') 
                                                                        
      write(88,11) step,smax,(Trmax(it),it=1,ntr),umax,vmax,wmax 
      write(89,11) step,smin,(Trmin(it),it=1,ntr),umin,vmin,wmin 
                                                                        
      close(88); close (89) 
   11 format(I5,9(f8.4)) 
                                                                        
!     FIND wup and wdn: the average vertical velocity in up and dn direc
                                                                        
      const= 1.d0/dble(NJ*NI) 
      do k=NK,1,-1 
         wup(k) = 0.d0 
         wdn(k) = 0.d0 
         do j=1,NJ 
            do i=1,NI 
               wup(k) = wup(k) + 0.5*(w(i,j,k,0) +abs(w(i,j,k,0))) 
               wdn(k) = wdn(k) + 0.5*(abs(w(i,j,k,0))-w(i,j,k,0)) 
            end do 
         end do 
         wup(k) = wup(k)*const 
         wdn(k) = wdn(k)*const 
      end do 
                                                                        
      open(unit=78,file='wup.out',access='append') 
      open(unit=79,file='wdn.out',access='append') 
                                                                        
      write(78,12) step,(wup(k),k=NK,1,-1) 
      write(79,12) step,(wdn(k),k=NK,1,-1) 
                                                                        
      close(78); close (79) 
   12 format(I5,64(E11.4)) 
                                                                        
                                                                        
      return 
      END                                           
