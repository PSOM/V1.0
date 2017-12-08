MODULE interface

  USE mod_param
  USE mod_vel
  USE mod_grid
  USE common2tracmass

  
  IMPLICIT NONE

CONTAINS
    subroutine initialize_tracmass
      
      IMT = NI; JMT = NJ ; KM = NK;

      allocate ( uflux(0:imt,jmt,km,nst), vflux(imt,0:jmt,km,nst) )
      allocate ( wflux(imt,jmt,0:km,nst) )
      allocate ( hs(imt,jmt,nst) )
      allocate ( kmt(imt,jmt))
      allocate ( dxv(imt,jmt),dyu(imt,jmt))
      allocate ( dz(km) )
      allocate ( dxdy(imt,jmt) )
      allocate ( dvol(imt,jmt,km) )
      allocate ( dzt(imt,jmt,km) )
      allocate ( trj(NTRACmax,7), nrj(NTRACmax,6) )

      hs = 0
      kmt = km

      dzt= DL/wz(1:NI,1:NJ,1:NK)
      dxv= LEN/ux(1:NI,1:NJ)
      dyu= LEN/vy(1:NI,1:NJ)
      dxdy= LEN*LEN*J2d(1:NI,1:NJ)
      dvol =LEN*LEN*DL*Jac(1:NI,1:NJ,1:NK)

    end subroutine initialize_tracmass


    subroutine set_tracmass_fluxes(step)
  
      INTEGER                                :: step
      
      uflux(:,:,:,step) = uf*LEN*DL
      vflux(:,:,:,step) = vf*LEN*DL
      wflux(:,:,:,step) = wf*LEN*LEN
  
    end subroutine set_tracmass_fluxes
  
end MODULE interface
