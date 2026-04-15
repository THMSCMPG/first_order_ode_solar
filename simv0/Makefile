#==============================================================================
#  Makefile  --  simv0 Integration Platform
#
#  Directory structure:
#    bin/        -- compiled executable  (bin/simv0)
#    obj/        -- object files (.o) and Fortran modules (.mod)
#    data/       -- output .dat files (written by the driver)
#    plots/      -- gnuplot scripts and generated .png files
#    configs/    -- example configuration files (.cfg)
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
#    make            - build executable to bin/simv0
#    make run        - build and run (outputs .dat files to data/)
#    make plots      - generate Gnuplot visualisations in plots/
#    make clean      - remove bin/ and obj/ directories
#    make clean-data - remove data/*.dat files
#    make clean-all  - remove all generated files
#==============================================================================

FC       := gfortran
FFLAGS   := -O2 -Wall -Wno-unused-variable
BIN_DIR  := bin
OBJ_DIR  := obj
DATA_DIR := data
PLOT_DIR := plots

TARGET   := $(BIN_DIR)/simv0

SRCS :=	simv0_config.f90          \
	simv0_method_selector.f90 \
	simv0_ode_suite.f90       \
	simv0_bte_coupling.f90    \
	simv0_spatial_1d.f90      \
	simv0_driver.f90

OBJS := $(SRCS:%.f90=$(OBJ_DIR)/%.o)

#  Default target
.PHONY: all run clean clean-data clean-all plots

all: $(TARGET)

# Create output directories as needed
$(BIN_DIR) $(OBJ_DIR) $(DATA_DIR):
	mkdir -p $@

$(TARGET): $(BIN_DIR) $(OBJ_DIR) $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(filter %.o,$^)

#  Per-module compilation rules (explicit order for .mod dependency chain)
$(OBJ_DIR)/simv0_config.o: simv0_config.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c -o $@ $< -J$(OBJ_DIR)

$(OBJ_DIR)/simv0_method_selector.o: simv0_method_selector.f90 $(OBJ_DIR)/simv0_config.o | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c -o $@ $< -J$(OBJ_DIR)

$(OBJ_DIR)/simv0_ode_suite.o: simv0_ode_suite.f90 $(OBJ_DIR)/simv0_config.o | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c -o $@ $< -J$(OBJ_DIR)

$(OBJ_DIR)/simv0_bte_coupling.o: simv0_bte_coupling.f90 $(OBJ_DIR)/simv0_config.o | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c -o $@ $< -J$(OBJ_DIR)

$(OBJ_DIR)/simv0_spatial_1d.o: simv0_spatial_1d.f90 $(OBJ_DIR)/simv0_config.o $(OBJ_DIR)/simv0_bte_coupling.o | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c -o $@ $< -J$(OBJ_DIR)

$(OBJ_DIR)/simv0_driver.o: simv0_driver.f90 $(OBJ_DIR)/simv0_config.o $(OBJ_DIR)/simv0_method_selector.o \
                            $(OBJ_DIR)/simv0_ode_suite.o $(OBJ_DIR)/simv0_bte_coupling.o $(OBJ_DIR)/simv0_spatial_1d.o | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c -o $@ $< -J$(OBJ_DIR)

# Run from data/ so .dat output files land in data/
run: $(TARGET) $(DATA_DIR)
	@echo "Running $(TARGET)..."
	cd $(DATA_DIR) && ../$(TARGET)
	@echo "✓ Output written to $(DATA_DIR)/"

# Run gnuplot from plots/ so ../data/ paths resolve correctly
plots:
	mkdir -p $(PLOT_DIR)
	@if command -v gnuplot >/dev/null 2>&1; then \
	    cd $(PLOT_DIR) && gnuplot plot_simv0.gp; \
	    echo "✓ Plots generated in $(PLOT_DIR)/"; \
	else \
	    echo "gnuplot not found; skipping plot generation"; \
	fi

# Remove build artifacts only (preserves data/)
clean:
	rm -rf $(BIN_DIR) $(OBJ_DIR)

# Remove data output files only
clean-data:
	rm -f $(DATA_DIR)/*.dat

# Remove all generated files
clean-all: clean clean-data
	rm -f $(PLOT_DIR)/*.png
