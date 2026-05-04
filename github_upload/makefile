# ==========================
# Fortran Makefile
# ==========================

# Compiler
FC = gfortran

# Compiler flags for standard build (debug + optimization)
FFLAGS = -c -g -O2 -Wall -fcheck=all

# Target executable
TARGET = disk_model

# Source files
SRCS = disk_global_constants.f90 \
       disk_params.f90 \
       disk_grid.f90 \
       disk_radmc_subroutines.f90 \
       disk_radmc_inpfiles.f90 \
       disk_density.f90 \
       disk_opacity.f90 \
       disk_mixopacity.f90 \
       disk_structure.f90 \
       disk_main_program.f90

# Object files
OBJS = $(SRCS:.f90=.o)

# ==========================
# Phony targets
# ==========================
.PHONY: all clean debug

# Default target: build everything
all: $(TARGET)

# Debug build: override FFLAGS
debug: FFLAGS = -c -g -O0 -Wall -fcheck=all -fbacktrace
debug: clean all

# Link the final executable
$(TARGET): $(OBJS)
	$(FC) $(OBJS) -o $(TARGET)

# ==========================
# Compilation rules
# ==========================
disk_global_constants.o: disk_global_constants.f90
	$(FC) $(FFLAGS) -o $@ $<

disk_params.o: disk_params.f90 disk_global_constants.o
	$(FC) $(FFLAGS) -o $@ $<

disk_grid.o: disk_grid.f90 disk_params.o
	$(FC) $(FFLAGS) -o $@ $<

disk_radmc_subroutines.o: disk_radmc_subroutines.f90
	$(FC) $(FFLAGS) -o $@ $<

disk_radmc_inpfiles.o: disk_radmc_inpfiles.f90 disk_radmc_subroutines.o
	$(FC) $(FFLAGS) -o $@ $<

disk_density.o: disk_density.f90 disk_global_constants.o
	$(FC) $(FFLAGS) -o $@ $<

disk_opacity.o: disk_opacity.f90 disk_global_constants.o disk_radmc_subroutines.o
	$(FC) $(FFLAGS) -o $@ $<

disk_mixopacity.o: disk_mixopacity.f90 disk_global_constants.o
	$(FC) $(FFLAGS) -o $@ $<

disk_structure.o: disk_structure.f90 disk_global_constants.o disk_params.o disk_grid.o disk_density.o disk_radmc_subroutines.o disk_radmc_inpfiles.o disk_opacity.o disk_mixopacity.o
	$(FC) $(FFLAGS) -o $@ $<

disk_main_program.o: disk_main_program.f90 disk_global_constants.o disk_params.o disk_grid.o disk_opacity.o disk_structure.o disk_radmc_inpfiles.o disk_mixopacity.o
	$(FC) $(FFLAGS) -o $@ $<

# ==========================
# Clean target
# ==========================
clean:
	-rm -f *.o *.mod $(TARGET)

