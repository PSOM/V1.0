subroutine mixing_vertical(var,vardif,m,step)
  !     ---------------------------------------------                     
#include "cppdefs.h"

  USE header
  !     Wind Stress can be specified here                                 
  !     Kz*du/dz = -Tx/rho,   Kz*dv/dz = -Ty/rho                          
  !     use level m                                                       
  !     computes d2s/dz2 at the cell centers.                             
  !     fricu and fricv are u,v diffusion terms (in header)               
  !     that need to be saved (if m.eq.0) for n2budget)                   

  implicit none 
  !REAL(kind=rc_kind) :: dSin,posneg,facreverse
  integer i,j,k,m,step,nday
  REAL(kind=rc_kind) :: zdep,yw,yset,ymid,ustar,Ekdepth,ff,zz     
                                                                   
  REAL(kind=rc_kind), dimension(0:NI+1,0:NJ+1, 0:NK+1) :: var 
  REAL(kind=rc_kind) :: vardif(NI,NJ,NK)

  REAL(kind=rc_kind) :: dvardzfc(0:NK), sw_source(NK)
  REAL(kind=rc_kind) :: Kdudzb,Kdvdzb,Kdudzt,Kdvdzt,wzkth,fact,fac_temp,heatfac,J_A1, J_A2
  REAL(kind=rc_kind) :: fac,facb,dfactor,facy,rhoinv,diff
  REAL(kind=rc_kind) :: wgt,day,ztransit,zextent,thy                                         
  REAL(kind=rc_kind), parameter :: Kzmax= 1.d-3, KzmaxTr=1.d-3 !  The viscosity is Kz*Kzmax 
  REAL(kind=rc_kind), parameter :: Cp= 4187.d0 !  Specific heat capacity of water

  INTEGER OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS,CHUNK,NPROC

  facb = RR*DL      ! Linear Drag                                                       
  ! facb = RR*DL*UL ! Quadratic drag                                                    

  fac= 1.d0/(UL*DL*delta) 
  fact = DL/UL 
  fac_temp= DL

  J_A1= J_A/J_lambda1             ! For shortwave penetration
  J_A2= (1 - J_A)/J_lambda2


  !Calculate the divergence of the diffusive flux with BOUNDARY CONDITIONS AND SOURCE TERMS
  !-------------------------------------
     if(selvar==1) then              ! TEMPERATURE 
        ! fluxes at top and bottom, source term for heating, and step
!$OMP PARALLEL DO PRIVATE(j,i,k,dvardzfc,dfactor) SHARED(vardif) COLLAPSE(2)

        do j=1,NJ 
           do i=1,NI 

       !CALCULATE Diffusive Fluxes across the k face of each cell 

              do k=1,NK-1 
                 dvardzfc(k)= wz(i,j,k)*KzTr* (var(i,j,k+1)- var(i,j,k)) 
              end do
              dvardzfc(NK)  = fac_temp *( - qloss(j) )/(R0*Cp)    ! Cp is the specific heat capacity of water
              dvardzfc(0)= 0.d0   ! no flux at bottom 

              sw_source= 0.d0

              ! Will prescribe the heat loss qloss as a flux at the surface (BC at top surface)
              ! Will prescribe the SW penetrative flux as a source by taking d/dz (of the penetrative heat flux).
              !      Kdfluxdzt = DL*( swr(j) - qloss(j) )/(R0*4187.d0)
              !Divergence of the penetrative heat flux is source
              heatfac= TL*swr(j)/(R0*Cp)     ! multiply by TL to make it  per non-dim time
              do k=1,NK
                 zz= zc(i,j,k)*DL
!               swrd(k) = swr(j)*( J_A*exp(zc(NI/2,NJ/2,k)*DL/J_lambda1) + (1 - J_A)*exp(zc(NI/2,NJ/2,k)*DL/J_lambda2)   )
                 sw_source(k) = heatfac *( J_A1 *exp(zz/J_lambda1) + J_A2 *exp(zz/J_lambda2)   )
               ! This is the heating source term (divergence of SW penetrative flux) prescribed at cell centers
               !        Kdfluxdzz(k) = DL*( swrd(k) )/(R0*4187.d0)
              end do

              do k=1,NK
                 vardif(i,j,k)= fac*Jac(i,j, k)* ( wz(i,j, k)*( dvardzfc(k) -dvardzfc(k-1))  + sw_source(k) )
              end do

           enddo ! i
        enddo ! j
!$OMP END PARALLEL DO                                                                   


     else if (selvar==2) then              ! SALINITY
!--------------------------------------------------------------------        
!$OMP PARALLEL DO PRIVATE(j,i,k,dvardzfc,dfactor) SHARED(vardif) COLLAPSE(2)

        do j=1,NJ
           do i=1,NI

              do k=1,NK-1 
                 dvardzfc(k)= wz(i,j,k)*KzTr* (var(i,j,k+1)- var(i,j,k)) 
              end do
              dvardzfc(0)= 0.d0   ! no flux at top and bottom is the default
              dvardzfc(NK)= 0.0 
              sw_source= 0.d0
              do k=1,NK
                 vardif(i,j,k)= fac*Jac(i,j, k)*wz(i,j, k)*( dvardzfc(k) -dvardzfc(k-1))  
              end do
           enddo ! i
        enddo ! j
!$OMP END PARALLEL DO                                                                   


     else if (selvar==3) then              ! U Momentum
!--------------------------------------------------------------------        
!        fac= 1.d0/(UL*DL*delta) 
!        fact = DL/UL 
!$OMP PARALLEL DO PRIVATE(j,i,k,dvardzfc,dfactor) SHARED(vardif) COLLAPSE(2)
        do j=1,NJ
           do i=1,NI

              do k=1,NK-1 
                 dvardzfc(k)= wz(i,j,k)*KzMom* (var(i,j,k+1)- var(i,j,k)) 
              end do

              rhoinv = 1.d0/rho(i,j,NK) 
              dvardzfc(NK)= stress_top_x(i,j)*rhoinv*fact 
              !     Quadratic drag at bottom
              !=            Kdudzb= facb*u(i,j,1,m)*abs(u(i,j,1,m))                  
              !     Linear Drag at bottom  
              dvardzfc(0)= facb*u(i,j,1,m) 
              do k=1,NK
                 vardif(i,j,k)= fac*Jac(i,j,k)*wz(i,j,k)*( dvardzfc(k) -dvardzfc(k-1)) 
              end do
           enddo ! i
        enddo ! j
!$OMP END PARALLEL DO   

     else if (selvar==4) then              ! V Momentum
!--------------------------------------------------------------------        
!$OMP PARALLEL DO PRIVATE(j,i,k,dvardzfc,dfactor) SHARED(vardif) COLLAPSE(2)
        do j=1,NJ
           do i=1,NI
              rhoinv = 1.d0/rho(i,j,NK) 
              dvardzfc(NK) = stress_top_y(i,j)*rhoinv*fact 
              do k=1,NK-1 
                 dvardzfc(k)= wz(i,j,k)*KzMom* (var(i,j,k+1)- var(i,j,k)) 
              end do
            !     Linear Drag at bottom  
              dvardzfc(0)= facb*v(i,j,1,m) 
              do k=1,NK
                 vardif(i,j,k)= fac*Jac(i,j, k)*wz(i,j,k)*( dvardzfc(k) -dvardzfc(k-1)) 
              end do
           enddo ! i
        enddo ! j
!$OMP END PARALLEL DO                                                                   
!--------------------------------------------------------------------        
     else if (selvar==5) then              ! W Momentum
!$OMP PARALLEL DO PRIVATE(j,i,k,dvardzfc,dfactor) SHARED(vardif) COLLAPSE(2)
        do j=1,NJ
           do i=1,NI
              dvardzfc(NK) = 0.d0
              dvardzfc(0)= 0.d0
              do k=1,NK-1 
                 dvardzfc(k)= wz(i,j,k)*KzMom* (var(i,j,k+1)- var(i,j,k)) 
              end do
            
              do k=1,NK
                 vardif(i,j,k)= fac*Jac(i,j, k)* wz(i,j, k)*( dvardzfc(k) -dvardzfc(k-1))
              end do
              if ((j.le.10).or.(j.gt.(NJ-10))) then                        ! DAMPING at boundaries
                 do k=1,NK
                    vardif(i,j,k)= vardif(i,j,k)  - fac*Jac(i,j, k)* var(i,j,k)* TL/43200.   ! 1/2 day restoring to zero w
                 end do
              end if

           enddo ! i
        enddo ! j
!$OMP END PARALLEL DO                                                                   
!--------------------------------------------------------------------        
     else if (selvar> 5.5) then              ! Tracers
!$OMP PARALLEL DO PRIVATE(j,i,k,dvardzfc,dfactor) SHARED(vardif) COLLAPSE(2)
        do j=1,NJ
           do i=1,NI
              dvardzfc(NK) = 0.d0
              dvardzfc(0)= 0.d0
              do k=1,NK-1 
                 dvardzfc(k)= wz(i,j,k)*KzTr* (var(i,j,k+1)- var(i,j,k)) 
              end do
              
              do k=1,NK
                 vardif(i,j,k)= fac*Jac(i,j, k)*wz(i,j, k)*( dvardzfc(k) -dvardzfc(k-1))  
              end do
           enddo ! i
        enddo ! j
!$OMP END PARALLEL DO                                                                   

     end if
     !--------Done with BC and source terms
     !-------------------
     ! AND COMPUTATION OF THE FLUX DIVERGENCE
     !-------------------

return 
END                                           
                                                                        
                                                                        
!       subroutine viscosity(dudz,dvdz,drdz,i,j)
! !     ---------------------------------------------
!         USE header
! !     Compute a Richardson number based vertical viscosity
! !     coefficient Kz using the final velocities u,v and density
! !     with the n=0 subscript
! !     Kz is situated at k cell faces
! !     The algorithm is from Rev Geophys., p 373, 1994. Large et al.
! !
! !     Assumes that rho is evaluated just prior to this call
!
!       integer i,j,k
!
!       integer n1
!       parameter (n1=3)
!       REAL(kind=rc_kind) :: fac,bvfreq,grho,RiCr,Ri,vshear
!       REAL(kind=rc_kind) :: dudz(0:NK),dvdz(0:NK),drdz(0:NK)
!       parameter (RiCr= 0.7d0)
!
!       grho= 9.81/R0
! !     fac= DL/(UL*UL)
!       DLinv = 1.0/DL
!       fac = UL*UL/(DL*DL)
! !
! !     Set kz(k=0)= 1. It is required only for the s,T equations, since
! !     for the momentum equations,  K*du/dz= ru.
! !     Kz is not needed at k=0 and NK since stress bcs are used.
!       Kz(i,j,0)= 0.d0
!       Kz(i,j,NK)= 0.d0
! !      if  (i.eq.16) write(100,*) 'j = ', j
!       do k=1,NK-1
!          bvfreq=  -grho*drdz(k)*DLinv
! !     BVfreq is N**2 and is in s^-2 if re-dim by DLinv
!          vshear= (dudz(k)*dudz(k) + dvdz(k)*dvdz(k))*fac
!          if (vshear.eq.0.d0) then
! !         if (vshear.le.1.d-12) then
!             Ri= 100.d0
!          else
!             Ri= bvfreq/vshear
!          end if
! !     unstable density profile => Ri < 0
!          if (Ri.le.0.d0) then
!             Kz(i,j,k)= 1.d0
!          else if (Ri.lt.RiCr) then
!             Kz(i,j,k)= (1.d0 - (Ri*Ri)/(RiCr*RiCr))**n1
!          else
!             Kz(i,j,k)= 0.d0
!          endif
! !         if ((i.eq.NI/2).and.(j.eq.NJ/2)) write(6,*) Ri ,bvfreq,vshear,
! !     &        Kz(i,j,k)
!       end do
!
! !     The value of Kz to be used is  Kz(k) *Kzmax
!
!       return
!       END
