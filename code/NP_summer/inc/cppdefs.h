#undef  runtracmass
#define periodic_ew
#undef  periodic_ns
#undef  allow_particle
#undef  rhoonly
#define relaxationflag
#undef  fixed_bottom_thickness
#define file_output
#define file_output_cdf
#define file_output_bin
#undef  gotm_call
#undef  implicit
#undef  parallel
! Stops the run after grid is calculated and generates output files containing the grid points.
#undef prerun
