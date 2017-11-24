# Makefile
# #
#FC = ifort     # Fortran 90 compiler
FC = gfortran     # Fortran 90 compiler


OBJ = utils.o lateration.o typedef.o makepdf.o k223d.o

FCFLAGS= -O2 
#FCFLAGS= -O -g -fbounds-check -fcheck=all -fbacktrace


#############################
.SUFFIXES: .f90

all: k223d.x

############################
# Objects
############################

utils.o: utils.f90 $(DEC)
	$(FC) $(FCFLAGS) -c $<

typedef.o: typedef.f90 $(DEC)
	$(FC) $(FCFLAGS) -c $<

lateration.o: lateration.f90 $(DEC)
	$(FC) $(FCFLAGS) -c $<

makepdf.o: makepdf.f90 $(DEC)
	$(FC) $(FCFLAGS) -c $<

k223d.o : k223d.f90
	$(FC) $(FCFLAGS) -c $<

k223d.x :$(OBJ)
	$(FC)  $(FFLAGS) -o $@ $(OBJ)


clean:
	rm *.o *.mod
