#!/bin/bash
FC	= mpif90
FFLAGS  = -O3
#FFLAGS  = -O3 -fallow-argument-mismatch
LFLAGS  = $(FFLAGS)
EXE     = fwi_run
SRC	= modules.f90 \
	00-main_Inversion.f90  \
	get_was_data.f90	\
	get_mcs_data.f90	\
	auxiliar.f90 \
	FWI.f90	\
	calculate_misfit_gradient.f90	\
	Misfit.f90	Misfit_L2.f90		Misfit_phase_env_Sn.f90		Misfit_CCTT.f90	\
	solver_acoustic_SS.f90  \
	solver_PML_matrices.f90 	solver_space_derivs.f90 	solver_PML_coeffs.f90  \
	solver_FD_coeffs.f90   	time_filter.f90 	Beta_calculous.f90  	filters.f90 \
	interpolation_akima.f90 	derivative.f90		polyfit.f90   	akima433.f90	fftpack.f90     \
	FA.f90  \
	tau_search.f90 CC_function.f90 fun_ceros.f90  interpola_lineal1D.f90	f1f2correlation.f90	\
	first_arrivals.f90	 \
	sub1_cycle.f90		sub2_cycle.f90		fun_DQ.f90	fourier.f90	\
	solver_forward_simple.f90	source_estimation.f90


# No need to edit below this line
.SUFFIXES: .f90 .o

OBJ = $(SRC:.f90=.o)

.f90.o:
	$(FC) $(FFLAGS) -c $<

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(LFLAGS) -o $@ $(OBJ)

$(OBJ): $(MF)

tar:
	tar cvf $(EXE).tar $(MF) $(SRC)

clean:
	rm -f $(OBJ) $(EXE) core

