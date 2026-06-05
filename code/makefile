# Makefile for k223d
#
# Usage:
#   make                -> standalone executable k223d.x (default)
#   make libk223d.so    -> shared library for Python wrapper
#   make all            -> both executable and shared library
#   make install        -> copy libk223d.so to the Python wrapper directory
#   make debug          -> debug build with full checks
#   make clean          -> remove compiled files
#
# Override the install destination:
#   make install PYTHON_LIB_DIR=/your/path
#
# Switch compiler by changing FC:
#   FC = gfortran  (default)
#   FC = ifort

FC = gfortran

PYTHON_LIB_DIR = ../examples/local_py_scripts

#------------------------------------------------------------------------------
# Debug flags: full runtime checks, no optimisation, stack trace on error
#------------------------------------------------------------------------------
FCFLAGS_DEBUG_GFORTRAN = -g -O0 \
                         -fcheck=all \
                         -ffpe-trap=invalid,zero,overflow \
                         -fbacktrace \
                         -Wall -Wextra \
                         -Wimplicit-interface

FCFLAGS_DEBUG_IFORT    = -g -O0 \
                         -check all \
                         -fpe0 \
                         -traceback \
                         -warn all

#------------------------------------------------------------------------------
# Production flags: high optimisation, current CPU, no fast-math
# -fPIC included always: negligible cost on x86-64 and required for
# libk223d.so without needing separate object files for each target.
#
# NOTE: -Ofast/-fast are intentionally avoided. Both enable -ffast-math which
#       breaks IEEE 754 compliance and can corrupt the infinity sentinel values
#       used in the distance and time calculations. -O3 + -march=native gives
#       nearly identical speed without these risks.
#------------------------------------------------------------------------------
FCFLAGS_PROD_GFORTRAN  = -O3 -march=native -funroll-loops -fPIC

FCFLAGS_PROD_IFORT     = -O3 -xHost -ipo -fPIC

#------------------------------------------------------------------------------
# Active flags: change FCFLAGS here or override on command line
#   e.g.  make FCFLAGS="$(FCFLAGS_DEBUG_GFORTRAN)"
#------------------------------------------------------------------------------
FCFLAGS = $(FCFLAGS_PROD_GFORTRAN)

#------------------------------------------------------------------------------
# Object list — order sets compilation sequence for module dependencies
#------------------------------------------------------------------------------
OBJ = generic.o      \
      lists.o        \
      LAT_mesh.o     \
      LAT_time.o     \
      LAT_source.o   \
      forparse.o     \
      mesh_geom.o    \
      distance.o     \
      time.o         \
      mesh_io.o       \
      source_model.o \
      k223d_core.o   \
      k223d.o

# Shared library uses same objects minus k223d.o, plus k223d_capi.o
CAPI_OBJ = generic.o      \
           lists.o        \
           LAT_mesh.o     \
           LAT_time.o     \
           LAT_source.o   \
           mesh_geom.o    \
           distance.o     \
           time.o         \
           mesh_io.o       \
           source_model.o \
           k223d_core.o   \
           k223d_capi.o

#------------------------------------------------------------------------------
# Default target: standalone executable only
#------------------------------------------------------------------------------
k223d.x: $(OBJ)
	$(FC) $(FCFLAGS) -o $@ $(OBJ)

# Build both executable and shared library
all: k223d.x libk223d.so

#------------------------------------------------------------------------------
# Compilation rules with explicit module dependencies
# Each target lists the .o files whose modules it USEs so make rebuilds
# correctly when an upstream module changes.
#------------------------------------------------------------------------------

generic.o: generic.f90
	$(FC) $(FCFLAGS) -c $<

lists.o: lists.f90 generic.o
	$(FC) $(FCFLAGS) -c $<

LAT_mesh.o: LAT_mesh.f90 generic.o lists.o
	$(FC) $(FCFLAGS) -c $<

LAT_time.o: LAT_time.f90 generic.o LAT_mesh.o
	$(FC) $(FCFLAGS) -c $<

LAT_source.o: LAT_source.f90 generic.o LAT_mesh.o
	$(FC) $(FCFLAGS) -c $<

forparse.o: forparse.f90 generic.o
	$(FC) $(FCFLAGS) -c $<

mesh_geom.o: mesh_geom.f90 generic.o LAT_mesh.o LAT_source.o distance.o LAT_time.o
	$(FC) $(FCFLAGS) -c $<

distance.o: distance.f90 generic.o LAT_mesh.o lists.o
	$(FC) $(FCFLAGS) -c $<

time.o: time.f90 generic.o LAT_mesh.o LAT_time.o
	$(FC) $(FCFLAGS) -c $<

mesh_io.o: mesh_io.f90 generic.o LAT_mesh.o LAT_source.o
	$(FC) $(FCFLAGS) -c $<

model_setup.o: model_setup.f90 generic.o LAT_mesh.o LAT_source.o
	$(FC) $(FCFLAGS) -c $<
      
source_model.o: source_model.f90 generic.o LAT_mesh.o LAT_source.o \
                mesh_geom.o distance.o time.o model_setup.o
	$(FC) $(FCFLAGS) -c $<

k223d_core.o: k223d_core.f90 generic.o LAT_mesh.o LAT_source.o \
              mesh_geom.o source_model.o model_setup.o
	$(FC) $(FCFLAGS) -c $<

k223d.o: k223d.f90 k223d_core.o LAT_source.o mesh_io.o forparse.o
	$(FC) $(FCFLAGS) -c $<

k223d_capi.o: k223d_capi.f90 k223d_core.o LAT_source.o LAT_mesh.o
	$(FC) $(FCFLAGS) -c $<

#------------------------------------------------------------------------------
# Link: standalone executable
#------------------------------------------------------------------------------
k223d.x: $(OBJ)
	$(FC) $(FCFLAGS) -o $@ $(OBJ)

#------------------------------------------------------------------------------
# Link: shared library for Python ctypes wrapper
#------------------------------------------------------------------------------
libk223d.so: $(CAPI_OBJ)
	$(FC) -shared $(FCFLAGS) -o $@ $(CAPI_OBJ)

#------------------------------------------------------------------------------
# Convenience targets
#------------------------------------------------------------------------------
debug:
	$(MAKE) FCFLAGS="$(FCFLAGS_DEBUG_GFORTRAN)"

production:
	$(MAKE) FCFLAGS="$(FCFLAGS_PROD_GFORTRAN)"

install: libk223d.so
	cp libk223d.so $(PYTHON_LIB_DIR)/libk223d.so

clean:
	rm -f *.o *.mod *.i90

.PHONY: all debug production install clean
