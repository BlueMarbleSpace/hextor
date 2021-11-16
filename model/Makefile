#
# Makefile for EBM
#

### Edit FC and WDIR to match your system ###
FC 	= /usr/bin/h5pfc
#FC 	= /opt/intel/bin/ifort
WDIR	= /models/hextor/model
#############################################

# Flags for gfortran
FLAGS    = -fdefault-real-8 -Wall -frecursive -O3 -g
FLAGSFFH = -cpp -Wall -frecursive -O3 -g

#Flags for ifort
#FLAGS	= -r8 -align all -i-static -warn nogeneral

#FLAGS	= -parallel -O2 -xHost -ipo -r8 -align all -i-static -warn nogeneral 
#FLAGS	= -O3 -xHost -ipo -r8 -align all -i-static -warn nogeneral
TARGET	= driver
OBJECTS = m_ffhash.o radiation.o driver.o

$(TARGET) : $(OBJECTS)
	$(FC) $(FLAGS) -o $(TARGET) $(OBJECTS)

m_ffhash.o : $(WDIR)/radiation/m_ffhash.f90
	$(FC) $(FLAGSFFH) -c $(WDIR)/radiation/m_ffhash.f90
radiation.o : $(WDIR)/radiation/radiation.f90
	$(FC) $(FLAGS) -c $(WDIR)/radiation/radiation.f90
driver.o : $(WDIR)/driver.f
	$(FC) $(FLAGS) -c $(WDIR)/driver.f

clean:
	rm -rf $(OBJECTS) $(TARGET)
new:
	$(MAKE) clean
	$(MAKE)
