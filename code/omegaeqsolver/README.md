# PSOM_omegaeqsolver
Compute the vertical velocities from the PSOM model outputs through the omega-equation.

This is an offline code that works with the strain output files of any PSOM experiment.

You can compile this code as an experiment with the tools/compile script but you should
replace the size.h file with the one from your experiment and eventually also the
routines that can affects the grid.

To run this code you should pass the namelist of the experiment in which you want to computed
the vertical velocities. The code will read the strain files from your output dir and it wil save files
named omega_XXXXX.cdf with the variable wqgeo.
