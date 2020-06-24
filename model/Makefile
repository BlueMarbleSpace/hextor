#
# Makefile for EBM
#

### Edit FC and WDIR to match your system ###
FC 	= /opt/intel/bin/ifort
WDIR	= /models/hextor/model
#############################################

FLAGS	= -r8 -align all -i-static -warn nogeneral
TARGET	= driver
OBJECTS = driver.o

$(TARGET) : $(OBJECTS)
	$(FC) $(FLAGS) -o $(TARGET) $(OBJECTS)

driver.o : $(WDIR)/driver.f
	$(FC) $(FLAGS) -c $(WDIR)/driver.f

clean:
	rm -rf $(OBJECTS) $(TARGET)

new:
	$(MAKE) clean
	$(MAKE)
