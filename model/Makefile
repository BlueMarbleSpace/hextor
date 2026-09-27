#
# Makefile for the HEXTOR driver
#
# FC must be an HDF5-aware Fortran wrapper (h5fc/h5pfc) around the compiler
# the flags are written for.  On merlin that is /opt/hdf5/bin/h5fc, an HDF5
# built against Intel ifx; source /opt/intel/oneapi/setvars.sh before building
# or running, or the driver dies at start-up looking for libimf.so.  On another
# machine override it on the command line (`make FC=h5pfc`) rather than editing
# this file.  Paths are relative to this directory, so nothing else is
# machine-specific.
FC = /opt/hdf5/bin/h5fc

# Flags for Intel ifx.  The published binaries were built with -parallel in
# this list, which ifx does not support and ignored with a warning (#10148),
# so it is dropped here; auto-parallelisation was never in effect.
FLAGS = -r8 -align all -warn nogeneral -O3 -g
# Flags for gfortran
#FLAGS = -fdefault-real-8 -Wall -frecursive -O3 -g

TARGET  = driver
OBJECTS = radiation.o driver.o

$(TARGET) : $(OBJECTS)
	$(FC) $(FLAGS) -o $(TARGET) $(OBJECTS)

radiation.o : radiation/radiation.f90
	$(FC) $(FLAGS) -c radiation/radiation.f90
# driver.f USEs the radiation module, so its .mod file must exist first.
driver.o : driver.f radiation.o
	$(FC) $(FLAGS) -c driver.f

clean:
	rm -f $(OBJECTS) $(TARGET) radiation_mod.mod
new:
	$(MAKE) clean
	$(MAKE)
