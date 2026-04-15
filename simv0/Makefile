#==============================================================================
#  Makefile  --  simv0 Integration Platform
#
#  Build order (module dependencies must be compiled first):
#    1. simv0_config.f90         (no Fortran module dependencies)
#    2. simv0_method_selector.f90 (uses simv0_config)
#    3. simv0_ode_suite.f90       (uses simv0_config)
#    4. simv0_bte_coupling.f90    (uses simv0_config)
#    5. simv0_spatial_1d.f90      (uses simv0_config, simv0_bte_coupling)
#    6. simv0_driver.f90          (uses all above)
#
#  Usage:
#    make           - build executable  ./simv0
#    make run       - build and run
#    make clean     - remove objects, .mod files, executable, data output
#    make plots     - generate Gnuplot visualisations (requires gnuplot)
#==============================================================================

FC       := gfortran
FFLAGS   := -O2 -Wall -Wno-unused-variable
TARGET   := simv0

SRCS :=	simv0_config.f90          \
	simv0_method_selector.f90 \
	simv0_ode_suite.f90       \
	simv0_bte_coupling.f90    \
	simv0_spatial_1d.f90      \
	simv0_driver.f90

OBJS := $(SRCS:.f90=.o)

#  Default target
.PHONY: all run clean plots

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

#  Per-module compilation rules (explicit order for .mod dependency chain)
simv0_config.o: simv0_config.f90
	$(FC) $(FFLAGS) -c $<

simv0_method_selector.o: simv0_method_selector.f90 simv0_config.o
	$(FC) $(FFLAGS) -c $<

simv0_ode_suite.o: simv0_ode_suite.f90 simv0_config.o
	$(FC) $(FFLAGS) -c $<

simv0_bte_coupling.o: simv0_bte_coupling.f90 simv0_config.o
	$(FC) $(FFLAGS) -c $<

simv0_spatial_1d.o: simv0_spatial_1d.f90 simv0_config.o simv0_bte_coupling.o
	$(FC) $(FFLAGS) -c $<

simv0_driver.o: simv0_driver.f90 simv0_config.o simv0_method_selector.o \
                simv0_ode_suite.o simv0_bte_coupling.o simv0_spatial_1d.o
	$(FC) $(FFLAGS) -c $<

run: $(TARGET)
	./$(TARGET)

plots: run
	@if command -v gnuplot >/dev/null 2>&1; then \
	    gnuplot plots/plot_simv0.gp; \
	    echo "Plots generated in plots/"; \
	else \
	    echo "gnuplot not found; skipping plot generation"; \
	fi

clean:
	rm -f $(OBJS) *.mod $(TARGET)
	rm -f simv0_decoupled.dat simv0_coupled.dat simv0_spatial.dat
	rm -f simv0_comparison.dat simv0_diagnostic.dat
