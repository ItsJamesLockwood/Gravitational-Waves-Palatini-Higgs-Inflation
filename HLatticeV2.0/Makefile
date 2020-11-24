############# gfortran ###############################
FC = gfortran
#FFLAGS = -cpp -fopenmp --free-line-length-none -fmax-identifier-length=63 -O2 -mcmodel=medium
FFLAGS = -cpp -fopenmp --free-line-length-none -fmax-identifier-length=63 -O2
####################################################
########### if you use ifort #########################
#FC = ifort
#FFLAGS = -fpp2 -openmp -O1 
#FFLAGS = -fpp2 -openmp -O1 -mcmodel=medium -i-dynamic -i8
########################################################

UTILOBJ = parameters.o  utils.o fft.o
MAINOBJ = define_fields.o model.o mutils.o io.o init.o evolve.o 
OBJ = $(UTILOBJ) $(MAINOBJ) 



################################################################
default: HLattice

%.o:%.f90
	$(FC) $(FFLAGS) -c $*.f90

HLattice: $(OBJ) driver.o
	$(FC) $(FFLAGS) $^ -o $@  	

*.o : configure.h

clean:
	rm -f *.o *.mod *.*~ \#*.*\# HLattice 2>/dev/null

wash:
	rm -f data/test* 

clear:
	rm -f data/*?.*
