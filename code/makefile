# Makefile
# #
# FC = ifort     # Fortran 90 compiler
FC = gfortran     # Fortran 90 compiler


# OBJ = utils.o lateration.o typedef.o makepdf.o k223d.o
# OBJ =  generic.o LAT_mesh.o LAT_source.o lists.o LAT_mesh_utils.o distance.o time.o typedef.o makepdf.o k223d.o
OBJ =  generic.o LAT_mesh.o LAT_time.o LAT_source.o lists.o LAT_mesh_utils.o distance.o time.o typedef.o makepdf.o k223d.o

#debug ifort:
# FCFLAGS= -g -check all -assume byterecl -traceback -debug all -fpe0
#run ifort
# FCFLAGS= -O3 -fast
#debug gfortran:
# FCFLAGS= -O -g -fbounds-check -ffpe-trap=invalid -fcheck=all -fbacktrace
#run gfortran
FCFLAGS= -Ofast

#############################
.SUFFIXES: .f90

all: k223d.x

############################
# Objects
############################

generic.o: generic.f90
						$(FC) $(FCFLAGS) -c $<

LAT_mesh.o: LAT_mesh.f90
						$(FC) $(FCFLAGS) -c $<

LAT_time.o: LAT_time.f90
						$(FC) $(FCFLAGS) -c $<

LAT_source.o: LAT_source.f90
						$(FC) $(FCFLAGS) -c $<

lists.o: lists.f90
						$(FC) $(FCFLAGS) -c $<

LAT_mesh_utils.o: LAT_mesh_utils.f90
						$(FC) $(FCFLAGS) -c $<

distance.o: distance.f90
						$(FC) $(FCFLAGS) -c $<

time.o: time.f90
						$(FC) $(FCFLAGS) -c $<

# utils.o: utils.f90 $(DEC)
# 	$(FC) $(FCFLAGS) -c $<

typedef.o: typedef.f90 $(DEC)
	$(FC) $(FCFLAGS) -c $<

# lateration.o: lateration.f90 $(DEC)
# 	$(FC) $(FCFLAGS) -c $<

makepdf.o: makepdf.f90 $(DEC)
	$(FC) $(FCFLAGS) -c $<

k223d.o : k223d.f90
	$(FC) $(FCFLAGS) -c $<

k223d.x :$(OBJ)
	$(FC)  $(FFLAGS) -o $@ $(OBJ)


clean:
	rm *.o *.mod *.i90 
