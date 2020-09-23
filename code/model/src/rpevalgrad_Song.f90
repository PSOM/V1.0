subroutine rpevalgrad_Song(n)

      USE header
      implicit none
      integer i,j,k,n
      REAL(kind=rc_kind) :: z,hpd,zt,zb,temp,rg,dep,ht
      real :: const,dz,Jsum
      REAL(kind=rc_kind) :: Jx(0:NI,0:NK),Jx2(NI,0:NK),Jy(0:NJ,0:NK),Jy2(NJ,0:NK)
      REAL(kind=rc_kind) :: dzxi(0:NI,0:NK),dzeta(NJ,0:NK),dzsig(NJ,0:NK),dzsig_x(0:NI,0:NK)

      if (rect.eqv.(.false.)) then 
!        need to modify grpifc,grpjfc to contain cross diff terms       
         write(6,*) 'modify grpifc,grpjfc' 
         stop 
      end if 
! c     Evaluate Jacobian at the edge-centers
! c     multiply press and press gradient by const
      const= 10.d0*gpr/P1
      call findzall 
      call evalrho(rho,n)

! ===================================================================================
! 
      do 10 i=0,NI+1 
         do 20 j=0,NJ+1 
            hpd= h(i,j)*HDL +D(i,j) 
            dep= D(i,j) 
            ht= h(i,j)*HDL 
!     for k=NK (firt for the toppmost layer of cells)                   
!     z is the dimensional z and is therefore multiplied by DL          
            k= NK 
            z= DL*zc(i,j,k) 
            zt= DL*zf(i,j,k) 
            zb= DL*zf(i,j,k-1) 
            rg= (rho(i,j,k) - R0)*const 
!            rg= (rhreal - R0)*const *R0/rhreal -non-Boussinesq for this
            rp(i,j,k)= rg*(zt-z) 
            temp= rg*(z-zb) 
!                                                                       
!****************                                                       
!            rp(i,j,k)= 0.d0                                            
!*******************                                                    
!     for the rest of the column                                        
!     remember that findz will not work for the k=0 cell where sig is -v
!     We will simply fill rp at k=0 using extrapolation - since this    
!     is what we use to fill s and T in any case.                       
!            do 30 k=NK-1,0,-1                                          
            do 30 k=NK-1,1,-1 
!     z is the dimensional z and is therefore multiplied by DL          
!               z= DL*findz(dfloat(k)-0.5,dztop,dep,ht,NK)              
               z= DL*zc(i,j,k) 
               zt= zb 
               zb= DL*zf(i,j,k-1) 
               rg= (rho(i,j,k) - R0)*const 

!     rg= (rhreal - R0)*const *R0/rhreal-non-Boussinesq for this term   
               rp(i,j,k)= rg*(zt-z) +temp +rp(i,j,k+1) 
               temp= rg*(z-zb) 
   30       continue 
!-            rp(i,j,0)= 2.d0*rp(i,j,1) -rp(i,j,2)                      
   20    continue 
   10 continue 
!                                                                       

!     Boundary conditions                                               
      do k=1,NK 
         do i=1,NI 
            rp(i,0,k)= rp(i,1,k) 
            rp(i,NJ+1,k)= rp(i,NJ,k) 
         end do 
!     periodic-ew boundaries                                            
         do j=0,NJ+1 
            rp(NI+1,j,k)= rp(1,j,k) 
            rp(0,j,k)= rp(NI,j,k) 
         end do 
      end do 
      
! ===================================================================================
! c     y-direction
! c     -----------
! c     sloping sigma surfaces
! c     solid boundaries

      do i=1,NI
         do j=1,NJ-1
            k=0
!            Jy(j,k)= rho(i,j+1,k+1) - rho(i,j,k+1) 
!            dzeta(j,k)= (zc(i,j+1,k+1) -zc(i,j,k+1) )*DL
            Jy(j,k)= 0.5*(rho(i,j+1,k+1)+rho(i,j+1,k+2) - (rho(i,j,k+1) +rho(i,j,k+2)) )
            dzeta(j,k)= 0.5*(zc(i,j+1,k+1) +zc(i,j+1,k+2) - (zc(i,j,k+1) + zc(i,j,k+2)) )*DL
            do k=1,NK-1
               Jy(j,k)= 0.5*(rho(i,j+1,k)+rho(i,j+1,k+1) - (rho(i,j,k) +rho(i,j,k+1) ))
               dzeta(j,k)= 0.5*(zc(i,j+1,k) +zc(i,j+1,k+1) - (zc(i,j,k) +zc(i,j,k+1) ) )*DL
            end do
            k= NK
            Jy(j,k)= 0.5*( rho(i,j+1,k)+rho(i,j+1,k-1) - (rho(i,j,k) +rho(i,j,k-1)))
            dzeta(j,k)= 0.5*(zc(i,j+1,k)+zc(i,j+1,k-1) - (zc(i,j,k)+zc(i,j,k-1)) )*DL
!            Jy(j,k)= rho(i,j+1,k) -rho(i,j,k)
!            dzeta(j,k)= (zc(i,j+1,k) -zc(i,j,k))*DL
         end do

         do j=1,NJ-1
            k=0
            Jy2(j,0)= 0.5*(rho(i,j,k+2) +rho(i,j+1,k+2) - (rho(i,j,k+1) + rho(i,j+1,k+1)) )
            dzsig(j,1)= 0.5*(zc(i,j,2) +zc(i,j+1,2) - (zc(i,j,1) + zc(i,j+1,1) ))*DL
            do k=1,NK-1
               Jy2(j,k)= 0.5*(rho(i,j,k+1) +rho(i,j+1,k+1) - (rho(i,j,k) + rho(i,j+1,k)) )
               dzsig(j,k)= 0.5d0*(zc(i,j,k+1) +zc(i,j+1,k+1) - (zc(i,j,k) +zc(i,j+1,k)) )*DL
            end do
            Jy2(j,NK)= 0.5*(rho(i,j,NK) +rho(i,j+1,NK) - (rho(i,j,NK-1) + rho(i,j+1,NK-1) ) )
     !     use same stencil for z as for rho
            dzsig(j,NK)= 0.5*(zc(i,j,NK) +zc(i,j+1,NK) - (zc(i,j,NK-1) + zc(i,j+1,NK-1) ) )*DL
         end do

         do k=0,NK
            do j=1,NJ-1
               Jy(j,k)= Jy(j,k)*dzsig(j,k)
               Jy2(j,k)= Jy2(j,k)*dzeta(j,k)
               Jy(j,k)= Jy(j,k) -Jy2(j,k)
            end do
         end do
         ! Solid boundaries
         do k=0,NK
            Jy(0,k)= 0.0
            Jy(NJ,k)= 0.0
         end do
         
         do j=1,NJ-1
            k= NK
            Jsum= Jy(j,k)*0.5
            grpjfc(i,j,NK)= Jsum*const*gj(i,j,k,2)
            do k=NK-1,1,-1
               Jsum= Jsum + Jy(j,k)
               grpjfc(i,j,k)= Jsum*const*gj(i,j,k,2)
            end do
         end do
         do k=1,NK
            grpjfc(i,0,k)= 0.d0
            grpjfc(i,NJ,k)= 0.d0
         end do
         do j=1,NJ
            k= NK
            drpy(i,j,k)= const*vy(i,j)*0.5 *( Jy(j,NK) + Jy(j-1,NK))*0.5
            do k=NK-1,1,-1
               drpy(i,j,k)= drpy(i,j,k+1) +const*vy(i,j)*0.5 *( Jy(j,k) + Jy(j-1,k))
            end do
         end do
      end do
      
! ===================================================================================
! c     x- direction
! c     ------------------------
! c     sloping sigma surfaces
! c     periodic boundaries

      do j=1,NJ
         do i=1,NI-1
            k=0
            Jx(i,k)= rho(i+1,j,k+1) - rho(i,j,k+1) 
            dzxi(i,k)= (zc(i+1,j,k+1) -zc(i,j,k+1) )*DL
            do k=1,NK-1
               Jx(i,k)= 0.5*(rho(i+1,j,k)+rho(i+1,j,k+1) - (rho(i,j,k) +rho(i,j,k+1) ))
               dzxi(i,k)= 0.5*(zc(i+1,j,k) +zc(i+1,j,k+1) - (zc(i,j,k) +zc(i,j,k+1) ) )*DL
            end do
            k= NK
            Jx(i,k)= rho(i+1,j,k) -rho(i,j,k)
            dzxi(i,k)= (zc(i+1,j,k) -zc(i,j,k))*DL
         end do
         i= NI
         k=0
         Jx(i,k)= rho(1,j,k+1) - rho(i,j,k+1) 
         dzxi(i,k)= (zc(1,j,k+1) -zc(i,j,k+1) )*DL
         do k=1,NK-1
            Jx(i,k)= 0.5*(rho(1,j,k)+rho(1,j,k+1) - (rho(i,j,k) +rho(i,j,k+1) ))
            dzxi(i,k)= 0.5*(zc(1,j,k) +zc(1,j,k+1) - (zc(i,j,k) +zc(i,j,k+1) ) )*DL
         end do
         k= NK
         Jx(i,k)= rho(1,j,k) -rho(i,j,k)
         dzxi(i,k)= (zc(1,j,k) -zc(i,j,k))*DL

         do i=1,NI-1
            k=0
            Jx2(j,0)= 0.5*(rho(i,j,k+2) +rho(i+1,j,k+2) - (rho(i,j,k+1) + rho(i+1,j,k+1)) )
            dzsig_x(i,1)= 0.5*(zc(i,j,2) +zc(i+1,j,2) - (zc(i,j,1) + zc(i+1,j,1) ))*DL
            do k=1,NK-1
               Jx2(i,k)= 0.5*(rho(i,j,k+1) +rho(i+1,j,k+1) - (rho(i,j,k) + rho(i+1,j,k)) )
               dzsig_x(i,k)= 0.5d0*(zc(i,j,k+1) +zc(i+1,j,k+1) - (zc(i,j,k) +zc(i+1,j,k)) )*DL
            end do
            Jx2(i,NK)= 0.5*(rho(i,j,NK) +rho(i+1,j,NK) - (rho(i,j,NK-1) + rho(i+1,j,NK-1) ) )
            dzsig_x(i,NK)= 0.5*(zc(i,j,NK) +zc(i+1,j,NK) - (zc(i,j,NK-1) + zc(i+1,j,NK-1) ) )*DL
         end do
         i=NI
         k=0
         Jx2(j,0)= 0.5*(rho(i,j,k+2) +rho(1,j,k+2) - (rho(i,j,k+1) + rho(1,j,k+1)) )
         dzsig_x(i,1)= 0.5*(zc(i,j,2) +zc(1,j,2) - (zc(i,j,1) + zc(1,j,1) ))*DL
         do k=1,NK-1
            Jx2(i,k)= 0.5*(rho(i,j,k+1) +rho(1,j,k+1) - (rho(i,j,k) + rho(1,j,k)) )
            dzsig_x(i,k)= 0.5d0*(zc(i,j,k+1) +zc(1,j,k+1) - (zc(i,j,k) +zc(1,j,k)) )*DL
         end do
         Jx2(i,NK)= 0.5*(rho(i,j,NK) +rho(1,j,NK) - (rho(i,j,NK-1) + rho(1,j,NK-1) ) )
         dzsig_x(i,NK)= 0.5*(zc(i,j,NK) +zc(1,j,NK) - (zc(i,j,NK-1) + zc(1,j,NK-1) ) )*DL

         do k=0,NK
            do i=1,NI
               Jx(i,k)= Jx(i,k)*dzsig_x(i,k)
               Jx2(i,k)= Jx2(i,k)*dzxi(i,k)
               Jx(i,k)= Jx(i,k) -Jx2(i,k)
            end do
         end do
         ! Periodic ew boundaries
         do k=0,NK
            Jx(0,k)= Jx(NI,k)
         end do
         
         do i=1,NI
            k= NK
            Jsum= Jx(i,k)*0.5
            grpifc(i,j,NK)= Jsum*const*gi(i,j,k,1)
            do k=NK-1,1,-1
               Jsum= Jsum + Jx(i,k)
               grpifc(i,j,k)= Jsum*const*gi(i,j,k,1)
            end do
         end do
         !     periodic ew boundary
         do k=1,NK
            grpifc(0,j,k)= grpifc(NI,j,k)
         end do
         do i=1,NI
            k= NK
            drpx(i,j,k)= const*ux(i,j)*0.5 *( Jx(i,NK) + Jx(i-1,NK))*0.5
            do k=NK-1,1,-1
               drpx(i,j,k)= drpx(i,j,k+1) +const*ux(i,j)*0.5 *( Jx(i,k) + Jx(i-1,k))
            end do
         end do
      end do


      return
      end

