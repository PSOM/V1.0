
##########
# Header
#---------

# Compiler
FC=ifort

# Libraries
LIB_DIR=-L/vortexfs1/apps/netcdf-intel/lib  
INCLUDE=-I/vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc -I/vortexfs1/home/mdever/.conda/envs/mat-env/include 
LNK_FLAGS=-lnetcdff -L/vortexfs1/apps/netcdf-intel/lib -L/vortexfs1/apps/zlib-1.2.11-intel/lib -L/vortexfs1/apps/szip-2.1.1-intel/lib -L/vortexfs1/apps/hdf-4.2.13-intel/lib -L/vortexfs1/apps/hdf5-1.10.2-intel/lib -L/vortexfs1/apps/jpeg-9c/lib -lnetcdf -lnetcdf

# Exectuable name
EXE=nh

###########
# The executable
#----------

$(EXE) : header.o particles.o relaxation.o main.o ini_setup.o ini_h.o ini_tracer.o write_bin_mod.o momentum.o potdens.o facediv.o intpol_3rdorder.o srcface_3rdorder.o vcenter.o vface.o advect.o advection_and_mixing.o cdiv.o cpcors.o cpfine.o efill.o mgrid.o prolong.o resid.o restrict.o mgpfill.o pbc.o sor.o spline.o seval.o sigma_toplayer.o staticsigma.o hsolve.o hfill.o chfine.o hbc.o vhydro.o evalrho_sT.o evalrho_rho.o evalrho.o psolve.o uvchy.o newcor.o newsrc.o coriolis.o hsave.o smooth.o setbc.o linerelax_periodicew.o dgtsl.o diag_energy.o findz_topmoves.o meanh.o rpevalgrad_Song.o correctbc.o conadjust_sT.o mixing_vertical.o checks.o ini_uv.o ran3.o mprove.o dgtsv.o pcorrect.o diag_vort.o calcskfc.o cfdiv.o mixing_horizontal.o mixing_isopycnal.o mixing_isopycnal_biharmonic.o srcface_nopy.o tracersource_nut.o diag_streamfunction.o diag_n2budget.o diag_n2.o diag_pv.o sigma2z.o ini_st.o utils.o write_bin.o ini_topog.o rpevalgrad_Sche.o findsigma.o wind_stress.o heat_flux.o uv_geostrophic.o write_cdf_2D_geopotential.o write_cdf.o write_cdf_1D_mooring.o write_cdf_2D_x.o write_cdf_2D_y.o write_cdf_2D_sigma.o write_cdf_2D_isopycnal.o write_cdf_3D_strain.o write_cdf_3D.o 
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -o $(EXE) header.o particles.o relaxation.o main.o ini_setup.o ini_h.o ini_tracer.o write_bin_mod.o momentum.o potdens.o facediv.o intpol_3rdorder.o srcface_3rdorder.o vcenter.o vface.o advect.o advection_and_mixing.o cdiv.o cpcors.o cpfine.o efill.o mgrid.o prolong.o resid.o restrict.o mgpfill.o pbc.o sor.o spline.o seval.o sigma_toplayer.o staticsigma.o hsolve.o hfill.o chfine.o hbc.o vhydro.o evalrho_sT.o evalrho_rho.o evalrho.o psolve.o uvchy.o newcor.o newsrc.o coriolis.o hsave.o smooth.o setbc.o linerelax_periodicew.o dgtsl.o diag_energy.o findz_topmoves.o meanh.o rpevalgrad_Song.o correctbc.o conadjust_sT.o mixing_vertical.o checks.o ini_uv.o ran3.o mprove.o dgtsv.o pcorrect.o diag_vort.o calcskfc.o cfdiv.o mixing_horizontal.o mixing_isopycnal.o mixing_isopycnal_biharmonic.o srcface_nopy.o tracersource_nut.o diag_streamfunction.o diag_n2budget.o diag_n2.o diag_pv.o sigma2z.o ini_st.o utils.o write_bin.o ini_topog.o rpevalgrad_Sche.o findsigma.o wind_stress.o heat_flux.o uv_geostrophic.o write_cdf_2D_geopotential.o write_cdf.o write_cdf_1D_mooring.o write_cdf_2D_x.o write_cdf_2D_y.o write_cdf_2D_sigma.o write_cdf_2D_isopycnal.o write_cdf_3D_strain.o write_cdf_3D.o $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


###########
# The modules
#----------

header.mod :  header.o /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/header.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/cppdefs.h /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/size.h /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/header_expe.h  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/header.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS) 

header.o :  /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/header.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/cppdefs.h /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/size.h /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/header_expe.h 
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/header.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)



particles.mod : header.mod particles.o /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/particles.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/cppdefs.h  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/particles.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS) 

particles.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/particles.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/cppdefs.h 
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/particles.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


###########
# The objects
#----------

relaxation.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/relaxation.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/relaxation.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


main.o : header.mod particles.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/main.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/cppdefs.h /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/main_declarations.h /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/ini_param.h 
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/main.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


ini_setup.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/ini_setup.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/cppdefs.h 
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/ini_setup.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


ini_h.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/ini_h.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/ini_h.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


ini_tracer.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/ini_tracer.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/ini_tracer.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


write_bin_mod.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/write_bin_mod.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/write_bin_mod.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


momentum.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/momentum.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/momentum.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


potdens.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/potdens.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/potdens.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


facediv.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/facediv.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/facediv.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


intpol_3rdorder.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/intpol_3rdorder.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/intpol_3rdorder.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


srcface_3rdorder.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/srcface_3rdorder.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/srcface_3rdorder.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


vcenter.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/vcenter.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/vcenter.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


vface.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/vface.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/vface.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


advect.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/advect.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/advect.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


advection_and_mixing.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/advection_and_mixing.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/cppdefs.h 
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/advection_and_mixing.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


cdiv.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/cdiv.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/cdiv.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


cpcors.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/cpcors.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/cpcors.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


cpfine.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/cpfine.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/cpfine.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


efill.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/efill.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/efill.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


mgrid.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/mgrid.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/mgrid.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


prolong.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/prolong.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/prolong.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


resid.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/resid.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/resid.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


restrict.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/restrict.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/restrict.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


mgpfill.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/mgpfill.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/mgpfill.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


pbc.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/pbc.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/pbc.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


sor.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/sor.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/sor.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


spline.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/spline.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/spline.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


seval.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/seval.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/seval.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


sigma_toplayer.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/sigma_toplayer.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/sigma_toplayer.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


staticsigma.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/staticsigma.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/cppdefs.h 
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/staticsigma.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


hsolve.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/hsolve.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/hsolve.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


hfill.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/hfill.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/hfill.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


chfine.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/chfine.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/chfine.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


hbc.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/hbc.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/hbc.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


vhydro.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/vhydro.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/vhydro.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


evalrho_sT.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/evalrho_sT.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/evalrho_sT.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


evalrho_rho.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/evalrho_rho.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/evalrho_rho.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


evalrho.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/evalrho.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/cppdefs.h 
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/evalrho.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


psolve.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/psolve.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/psolve.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


uvchy.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/uvchy.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/uvchy.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


newcor.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/newcor.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/newcor.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


newsrc.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/newsrc.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/newsrc.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


coriolis.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/coriolis.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/coriolis.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


hsave.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/hsave.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/hsave.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


smooth.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/smooth.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/cppdefs.h 
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/smooth.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


setbc.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/setbc.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/setbc.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


linerelax_periodicew.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/linerelax_periodicew.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/linerelax_periodicew.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


dgtsl.o :  /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/dgtsl.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/dgtsl.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


diag_energy.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/diag_energy.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/diag_energy.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


findz_topmoves.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/findz_topmoves.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/cppdefs.h 
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/findz_topmoves.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


meanh.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/meanh.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/meanh.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


rpevalgrad_Song.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/rpevalgrad_Song.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/rpevalgrad_Song.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


correctbc.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/correctbc.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/correctbc.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


conadjust_sT.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/conadjust_sT.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/conadjust_sT.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


mixing_vertical.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/mixing_vertical.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/cppdefs.h 
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/mixing_vertical.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


checks.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/checks.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/cppdefs.h 
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/checks.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


ini_uv.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/ini_uv.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/ini_uv.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


ran3.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/ran3.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/ran3.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


mprove.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/mprove.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/mprove.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


dgtsv.o :  /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/dgtsv.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/dgtsv.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


pcorrect.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/pcorrect.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/pcorrect.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


diag_vort.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/diag_vort.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/diag_vort.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


calcskfc.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/calcskfc.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/calcskfc.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


cfdiv.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/cfdiv.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/cfdiv.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


mixing_horizontal.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/mixing_horizontal.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/mixing_horizontal.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


mixing_isopycnal.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/mixing_isopycnal.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/mixing_isopycnal.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


mixing_isopycnal_biharmonic.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/mixing_isopycnal_biharmonic.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/mixing_isopycnal_biharmonic.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


srcface_nopy.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/srcface_nopy.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/srcface_nopy.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


tracersource_nut.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/tracersource_nut.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/tracersource_nut.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


diag_streamfunction.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/diag_streamfunction.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/diag_streamfunction.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


diag_n2budget.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/diag_n2budget.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/diag_n2budget.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


diag_n2.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/diag_n2.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/diag_n2.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


diag_pv.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/diag_pv.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/diag_pv.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


sigma2z.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/sigma2z.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/sigma2z.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


ini_st.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/ini_st.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/ini_st.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


utils.o :  /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/utils.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/utils.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


write_bin.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/write_bin.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/cppdefs.h 
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/write_bin.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


ini_topog.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/ini_topog.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/cppdefs.h 
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/ini_topog.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


rpevalgrad_Sche.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/rpevalgrad_Sche.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/rpevalgrad_Sche.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


findsigma.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/findsigma.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/findsigma.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


wind_stress.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/wind_stress.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/wind_stress.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


heat_flux.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/heat_flux.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/heat_flux.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


uv_geostrophic.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/uv_geostrophic.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/uv_geostrophic.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


write_cdf_2D_geopotential.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/write_cdf_2D_geopotential.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/netcdf.inc 
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/write_cdf_2D_geopotential.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


write_cdf.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/write_cdf.f90  
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/write_cdf.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


write_cdf_1D_mooring.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/write_cdf_1D_mooring.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/netcdf.inc 
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/write_cdf_1D_mooring.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


write_cdf_2D_x.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/write_cdf_2D_x.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/netcdf.inc 
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/write_cdf_2D_x.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


write_cdf_2D_y.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/write_cdf_2D_y.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/netcdf.inc 
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/write_cdf_2D_y.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


write_cdf_2D_sigma.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/write_cdf_2D_sigma.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/netcdf.inc 
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/write_cdf_2D_sigma.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


write_cdf_2D_isopycnal.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/write_cdf_2D_isopycnal.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/netcdf.inc 
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/write_cdf_2D_isopycnal.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


write_cdf_3D_strain.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/write_cdf_3D_strain.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/netcdf.inc 
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/write_cdf_3D_strain.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


write_cdf_3D.o : header.mod /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/write_cdf_3D.f90 /vortexfs1/share/mahadevanlab/PSOM/code/model_build/inc/netcdf.inc 
	$(FC)  -profile-functions -fpp -real-size 64 -mcmodel=large -qopenmp -stand f03 -O3 -u -c /vortexfs1/share/mahadevanlab/PSOM/code/model_build/src/write_cdf_3D.f90 $(INCLUDE) $(LIB_DIR) $(LNK_FLAGS)


##########
# Ending
#---------

clean :
	rm *.o *.mod

cleanmore :
	rm *.o *.mod $(EXE)

