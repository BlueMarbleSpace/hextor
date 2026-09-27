#
# Makefile for the standalone radiation-module probes: radtest (main.f90),
# compare (compare.f90) and tabletest (tabletest.f90, driven by
# tools/check_table_reader.py).
#
# FC and the flag rules are as in ../Makefile: an HDF5-aware wrapper around
# Intel ifx on merlin, overridden on the command line elsewhere.  Paths are
# relative to this directory.  `make` only relinks a probe when radiation.f90
# is newer than its object, so run `make clean` first after touching the
# module, or a stale probe validates a stale check (see CLAUDE.md).
FC = /opt/hdf5/bin/h5fc

# Flags for Intel ifx
FLAGS = -r8 -align all -warn nogeneral
# Flags for gfortran
#FLAGS = -fdefault-real-8 -Wall -frecursive -O3 -g

TARGET  = radtest
OBJECTS = radiation.o main.o

$(TARGET) : $(OBJECTS)
	$(FC) $(FLAGS) -o $(TARGET) $(OBJECTS)

compare : radiation.o compare.o
	$(FC) $(FLAGS) -o compare radiation.o compare.o

tabletest : radiation.o tabletest.o
	$(FC) $(FLAGS) -o tabletest radiation.o tabletest.o

radiation.o : radiation.f90
	$(FC) $(FLAGS) -c radiation.f90
main.o : main.f90 radiation.o
	$(FC) $(FLAGS) -c main.f90
compare.o : compare.f90 radiation.o
	$(FC) $(FLAGS) -c compare.f90
tabletest.o : tabletest.f90 radiation.o
	$(FC) $(FLAGS) -c tabletest.f90

clean:
	rm -f $(OBJECTS) $(TARGET) compare compare.o tabletest tabletest.o radiation_mod.mod
new:
	$(MAKE) clean
	$(MAKE)
