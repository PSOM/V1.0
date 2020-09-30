#undef   runtracmass
#define  periodic_ew
#undef   allow_particle
!#undef   rhoonly
#define   rhoonly

#undef   relaxation
#define  fixed_bottom_thickness
#undef   gotm_call
#define  file_output
#define  file_output_cdf
#undef   file_output_bin
#undef   implicit
   
   !#define  biharmonic_horizontal
#undef   biharmonic_horizontal
   
#undef  sigma_stretch
!#define  sigma_stretch