#
# Makefile for EBM
#

### Edit FC and WDIR to match your system ###
FC 	= /usr/bin/h5pfc
#FC 	= /opt/intel/bin/ifort
WDIR	= /models/hextor/model
#############################################

FLAGS	= -parallel -O2 -xHost -ipo -r8 -align all -i-static -warn nogeneral 
#FLAGS	= -O3 -xHost -ipo -r8 -align all -i-static -warn nogeneral
TARGET	= driver
OBJECTS = radiation.o driver.o

$(TARGET) : $(OBJECTS)
	$(FC) $(FLAGS) -o $(TARGET) $(OBJECTS)

radiation.o : $(WDIR)/radiation/radiation.f90
	$(FC) $(FLAGS) -c $(WDIR)/radiation/radiation.f90

driver.o : $(WDIR)/driver.f
	$(FC) $(FLAGS) -c $(WDIR)/driver.f

clean:
	rm -rf $(OBJECTS) $(TARGET)

new:
	$(MAKE) clean
	$(MAKE)
