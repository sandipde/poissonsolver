#__________________________________________________________________________
#   PROGRAM BY SANDIP DE      26.01.09
#__________________________________________________________________________
# Start by clearing the list of suffixes
.SUFFIXES:

# Default suffix list
.SUFFIXES: .f90 .f .o


# List of object files
OBJS =  modules.o              \
 	poisson_allocate_deallocate.o \
        poisson.o			\
        read_input.o            \
        gen_mesh.o       \
	charge_density.o \
	double_grid_trans.o \
	fft_grid_trans.o      \
	fft.o             \
	gen_hartree_pot.o \
	write_in_vmd_format.o \
 	del2.o \
	check_result.o 
EXECNAME =poisson.x
FC = gfortran 

#FCFLAGS = -r8  -O3 -i8 -tpp7 -axW -xW -ftz -lowercase -mp1 -complex_limited_range -cm -w -ipo
FCFLAGS = -O3 
LFLAGS= -lfftw3 -lm  -lsvml
# Default rule for making a .o file from a .f90 file
.f90.o:
	$(FC) -c $(FCFLAGS) $<

# Perform a complete compilation (**default action**)
all: $(EXECNAME)

$(EXECNAME): $(OBJS)
	$(FC) -o $(EXECNAME) $(OBJS) $(LFLAGS)
clean:
	rm -f **.o **.mod
distclean:
	rm -f **.o **.mod *.dat *.x
