# Compiler and flags
FC=gfortran
FFLAGS=-O2 -g -Wall -fcheck=all -cpp
OBJDIR=obj

# Program name
PROGRAM = radiative_transfer

# List of module files
MODULES = configure.f90 \
          common_grid.f90 \
          common_dust.f90 \
          common_boundary.f90 \
          common_montecarlo.f90 \
          common_montecarloarrays.f90 \
          common_diffusion.f90 \
          common_vstruct.f90 \
          numerical_receipe.f90

# List of source files
SOURCES = $(MODULES) \
          dust_main.f90 \
          diffusion.f90 \
          vertical_structure.f90 \
          montecarlo.f90 \
          main.f90 \
          main_program.f90

# Object files (OBJDIR prefixed)
OBJECTS = $(SOURCES:.f90=.o)

# Default rule
all: $(PROGRAM)

# Compile the program
$(PROGRAM): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $(OBJECTS)

# Compile .f90 to .o
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# Clean rule
clean:
	rm -f *.o *.mod $(PROGRAM)

