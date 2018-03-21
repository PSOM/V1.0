subroutine sigma 
!     ------------------------                                          
  USE header
!     ONLY FOR THE TOPMOST GRID CELL                                    
!     modified a for periodicew bc                                      
!     This subroutine updates all those quantities that are a function  
!     of time and  (time and depth). It is called  every time step.     
!     The surface bc  wfbct is updated here.                            
  implicit none
  integer i,j,k 
  REAL(kind=rc_kind) :: hpd,hpdinv,hu,hv,hx,hy,hxpdx,hypdy,z,temp,be2,wxk,wyk,d2,dnk,dnkm1,sig,zpd
  REAL(kind=rc_kind) :: g13(0:NI+1,0:NJ+1,0:NK+1),g23(0:NI+1,0:NJ+1,0:NK+1)          
!     Ddx and Ddy are known from init                                   

 ! ------------------- 

   dnkm1= dble(NK-1) 
   dnk= dble(NK) 
                                                                        
   be2= beta*EPS*EPS 

                                                                        
 ! ------------------- 

  do j=0,NJ+1 
    do i=0,NI+1 
      !     All these variables are functions of time                         
      hpd= h(i,j)*HDL +dztop 
      hpdinv= 1.d0/hpd 

      ! Computation of hu:
      if (i.eq.0) then 
        hu= HDL*(h(i+1,j)-h(NI-1,j)) 
       else if  (i.eq.NI+1) then 
        hu= HDL*(h(2,j)-h(i-1,j)) 
       else 
        hu= 0.5d0*HDL*(h(i+1,j)-h(i-1,j)) 
      endif 

      ! Computation of hv:
      if (j.eq.0) then 
        hv= HDL*(h(i,j+1)-h(i,j)) 
       else if (j.eq.NJ+1) then 
        hv= HDL*(h(i,j)-h(i,j-1)) 
       else 
        hv= 0.5d0*HDL*(h(i,j+1)-h(i,j-1)) 
      endif 

      ! Computation of hx and hy:
      hx= hu*ux(i,j) + hv*vx(i,j) 
      hy= hu*uy(i,j) + hv*vy(i,j) 

      !     wz is not a function of depth when the sigma lines are equally spa
      !     Hence wz is wz(i,j,time). For a stretched grid wz would be w(i,j,k
      !     Then Jac which is now Jac(i,j,time) would  become  Jac(i,j,k,time)

      ! Computation of wx, wy, g13, g23:                                                                        
      do k=NK,NK+1 
        wz(i,j,k)= hpdinv 
        Jac(i,j,k)= J2d(i,j)/wz(i,j,k) 
        !     All these variables are functions of time and depth               
        !     now hdt computed in vhydro already contains HDL                   
        sig= dble(k)-0.5d0 
        z= (sig -dnkm1)*hpd -dztop                        ! z is dimensionless and is equal to  z*/DL                       
        wx(i,j,k)= (dnkm1-sig)*hx*hpdinv 
        wy(i,j,k)= (dnkm1-sig)*hy*hpdinv 
        g13(i,j,k)= ux(i,j)*wx(i,j,k) +uy(i,j)*wy(i,j,k) 
        g23(i,j,k)= vx(i,j)*wx(i,j,k) +vy(i,j)*wy(i,j,k) 
        ! g33(i,j,k)= wx(i,j,k)*wx(i,j,k) +wy(i,j,k)*wy(i,j,k) +  be2*wz(i,j)*wz(i,j)                                
      enddo

      ! Computation of wt (useless...):
      do k=NK-1,NK 
         sig= dble(k) 
         wt(i,j,k)= (dnkm1-sig)*hdt(i,j)*hpdinv  !     wt is defined at cell faces                                       
      end do 

      wzk(i,j,NK)= hpdinv 
      wzk(i,j,NK-1)= 0.5*(wz(i,j,NK) + wz(i,j,NK-1)) 
                                                                  
      !     We evaluate gk at i=0,NI+1,j=0,NJ+1 only because these are used   
      !     at k=0,NK to fill the horizontal edges in mgpfill or hnbc.        
      do k=NK-1,NK 
        sig= dble(k) 
        wxk= (dnkm1-sig)*hx*hpdinv 
        wyk= (dnkm1-sig)*hy*hpdinv 
        gqk(i,j,k,1)= qpr*Jac(i,j,k)*(ux(i,j)*wxk +uy(i,j)*wyk) 
        gqk(i,j,k,2)= qpr*Jac(i,j,k)*(vx(i,j)*wxk +vy(i,j)*wyk) 
        gqk(i,j,k,3)= Jac(i,j,k)*(qpr*(wxk*wxk +wyk*wyk) + be2*wz(i,j,k)*wz(i,j,k))                            
      enddo

    enddo ! <- i
  enddo ! <- j  

                                                                        
  !- ======                                                            
     k= NK 
  !- ======                                                            

  do i=0,NI 
    do j=1,NJ 
      Jifc(i,j,k)= 0.5d0*(Jac(i,j,k)+ Jac(i+1,j,k)) 
      gi(i,j,k,1)= 0.5d0*(g11(i,j) +g11(i+1,j))*Jifc(i,j,k) 
      gi(i,j,k,2)= 0.5d0*(g12(i,j) +g12(i+1,j))*Jifc(i,j,k) 
      gqi(i,j,k,1)= qpr*gi(i,j,k,1) 
      gqi(i,j,k,2)= qpr*gi(i,j,k,2) 
    enddo
  enddo
                                                                    
  do i=1,NI 
    do j=0,NJ 
      Jjfc(i,j,k)= 0.5d0*(Jac(i,j,k)+ Jac(i,j+1,k)) 
      gj(i,j,k,1)= 0.5d0*(g12(i,j) +g12(i,j+1))*Jjfc(i,j,k) 
      gj(i,j,k,2)= 0.5d0*(g22(i,j) +g22(i,j+1))*Jjfc(i,j,k) 
      gqj(i,j,k,1)= qpr*gj(i,j,k,1) 
      gqj(i,j,k,2)= qpr*gj(i,j,k,2) 
    enddo
  enddo
                                                                    
  do j=1,NJ 
    do i=0,NI 
      gi3(i,j,k)= 0.5d0*(g13(i,j,k) +g13(i+1,j,k))*Jifc(i,j,k) 
      gqi3(i,j,k)= qpr*gi3(i,j,k) 
    enddo
  enddo

  do j=0,NJ 
    do i=1,NI 
      gj3(i,j,k)= 0.5d0*(g23(i,j,k) +g23(i,j+1,k))*Jjfc(i,j,k) 
      gqj3(i,j,k)= qpr*gj3(i,j,k) 
    enddo
  enddo

 do i=0,NI
   do j=0,NJ
     do k=NK,NK+1
       Jacinv(i,j,k)=1./Jac(i,j,k)
     enddo
   enddo
 enddo


return 
END                                           
