# Makefile
# #
FC = gfortran     # Fortran 90 compiler


OBJ = make_mesh.o

#FCFLAGS= -O -g -fbounds-check -fcheck=all fbacktrace
FCFLAGS= -O2


#############################
.SUFFIXES: .f90

all: make_mesh.x 

############################
# Objects
############################

make_mesh.o : make_mesh.f90
	$(FC) $(FCFLAGS) -c $< 

make_mesh.x :$(OBJ)
	$(FC)  $(FFLAGS) -o $@ $(OBJ)


clean:
	rm *.o





