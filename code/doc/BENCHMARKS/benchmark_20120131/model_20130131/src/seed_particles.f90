subroutine seed_particles
  !---------------------
  USE interface

  
  trj(1,1) = IMT/2          ! Current grid position  x
  trj(1,2) = JMT/2          ! Current grid position  y
  trj(1,3) = KM             ! Current grid position  z
  trj(1,4) = 0.             ! tt
  trj(1,5) = 1.             ! subvol
  trj(1,6) = 1.             ! arct
  trj(1,7) = 1.             ! t0
  nrj(1,1) = IMT/2          ! grid index i
  nrj(1,2) = JMT/2          ! grid index j
  nrj(1,3) = KM             ! grid index k
  nrj(1,4) = 1.             ! niter
  nrj(1,5) = 1              ! ts
end subroutine seed_particles
