TEST 	 = src/main
NAME2D = src/gturb_f2d2c
NAME3D = src/gturb_f3d3c
# FC       = mpif90
FC       = nvfortran
EXE	 = exe

# FCFLAGS  = -fast -gpu=managed,cc75,cuda11.5.0 -acc=gpu -cuda -cudalib=cufft -Minfo=accel -Mpreprocess -D_USE_NVTX -lnvToolsExt -mcmodel=medium
# FCFLAGS  = -fast -acc=gpu -cuda -cudalib=cufft -Minfo=all -Mpreprocess -mcmodel=medium#-fPIC
FCFLAGS  = -fast -acc=gpu -cuda -cudalib=cufft -mcmodel=medium -Minfo=all

# -ta=tesla:managed
# FCFLAGS  = -cudalib=cufft

SRCFILES = src/modules.f90 src/defk.f90 src/trunc.f90 \
					 src/nlt.f90 src/produ.f90 src/diagnostics.f90 src/rkstep.f90 \
					 src/read.f90 src/write.f90 src/inieuler.f90 \
					 src/fft_inv.f90 src/fft_dir.f90 \
					 
SRC2D = src/iniforcing_2d2c.f90 src/forcing_2d2c.f90
SRC3D = src/iniforcing_3d3c.f90 src/forcing_3d3c.f90

all: 3d move

2d: $(TEST).f90
	$(FC) $(FCFLAGS) -o $(NAME2D) $(SRCFILES) $(SRC2D) $<
	@echo 'Cleaning up...'
	@rm -rf *.o *.mod
	
3d: $(TEST).f90
	$(FC) $(FCFLAGS) -o $(NAME3D) $(SRCFILES) $(SRC3D) $<
	@echo 'Cleaning up...'
	@rm -rf *.o *.mod

move:
	@mkdir -p gturb gturb/fields gturb/files
# 	@mkdir -p gturb/Diag/fields gturb/Diag/spectra gturb/Diag/fluxes
# 	@mkdir -p gturb/history-slurm
	@scp src/jobscript gturb/jobscript
	@scp src/gturb* gturb/
	@scp src/seed.0 gturb/files/seed.000
	@scp src/startframe.dat gturb/curframe.dat
	@scp src/params.dat gturb/params.dat
	@scp src/reset.sh gturb/reset.sh
	@rm src/gturb*

verify:

clean:
	@echo 'Cleaning up...'
	@rm -rf *.o *.dwf *.pdb *.mod prof gturb/

# TEST 	   = src/main
# FC       = nvfortran
# FCFLAGS  = -fast -acc=gpu -cuda -cudalib=cufft -mcmodel=medium -Minfo=all
# 
# SRCFIL2D = src/modules.f90 src/defk.f90 src/iniforcing_2d2c.f90 src/trunc.f90 \
# 					 src/nlt.f90 src/produ.f90 src/diagnostics.f90 src/rkstep.f90 \
# 					 src/read.f90 src/write.f90 src/forcing_2d2c.f90 src/inieuler.f90 \
# 					 src/fft_inv.f90 src/fft_dir.f90 \
# 					 
# SRCFIL3D = src/modules.f90 src/defk.f90 src/iniforcing_3d3c.f90 src/trunc.f90 \
# 					 src/nlt.f90 src/produ.f90 src/diagnostics.f90 src/rkstep.f90 \
# 					 src/read.f90 src/write.f90 src/forcing_3d3c.f90 src/inieuler.f90 \
# 					 src/fft_inv.f90 src/fft_dir.f90 \
# 
# 3d: $(TEST).f90
# 	$(FC) $(FCFLAGS) -o $(gturb)$(3d3c) $(SRCFIL3D) $<
# 	@echo 'Cleaning up...'
# 	@rm -rf *.o *.mod
# 	
# 2d: $(TEST).f90
# 	$(FC) $(FCFLAGS) -o $(gturb)$(2d2c) $(SRCFIL2D) $<	
# 	@echo 'Cleaning up...'
# 	@rm -rf *.o *.mod
# 	
# mv:
# 	@mkdir -p gturb gturb/Diag gturb/Seed
# 	@mkdir -p gturb/Diag/fields gturb/Diag/spectra gturb/Diag/fluxes
# 	@mkdir -p gturb/history-slurm
# 	@scp src/jobscript gturb/jobscript
# 	@scp src/gturb3d3c gturb/gturb3d3c
# 	@scp src/gturb2d2c gturb/gturb2d2c
# 	@scp src/seed.0 gturb/Seed/seed.00000000
# 	@scp src/startframe.dat gturb/curframe.dat
# 	@scp src/params.dat gturb/params.dat
# 	@scp src/reset.sh gturb/reset.sh
# 	
# run:
# 	# $(TEST).$(EXE)
# # 	$(RUN) ./$(TEST).$(EXE)
# 
# verify:
# 
# 
# clean:
# 	@echo 'Cleaning up...'
# 	@rm -rf *.$(EXE) *.o *.dwf *.pdb *.mod prof