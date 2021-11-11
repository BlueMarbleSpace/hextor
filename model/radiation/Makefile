#
# Makefile for testing radiation module
#

### Edit FC and WDIR to match your system ###
FC 	= /usr/bin/h5pfc
WDIR	= /models/hextor/model/radiation
#############################################

# Flags for gfortran
FLAGS    = -fdefault-real-8 -Wall -frecursive -O3 -g
FLAGSFFH = -cpp -Wall -frecursive -O3 -g

#Flags for ifort
#FLAGS	= -r8 -align all -i-static -warn nogeneral

TARGET	= radtest
OBJECTS	= m_ffhash.o radiation.o main.o

$(TARGET) : $(OBJECTS)
	$(FC) $(FLAGS) -o $(TARGET) $(OBJECTS)


m_ffhash.o : $(WDIR)/m_ffhash.f90
	$(FC) $(FLAGSFFH) -c $(WDIR)/m_ffhash.f90
radiation.o : $(WDIR)/radiation.f90
	$(FC) $(FLAGS) -c $(WDIR)/radiation.f90
main.o : $(WDIR)/main.f90
	$(FC) $(FLAGS) -c $(WDIR)/main.f90

clean:
	rm -rf $(OBJECTS) $(TARGET)
	rm radiation_mod.mod
new:
	$(MAKE) clean
	$(MAKE)
