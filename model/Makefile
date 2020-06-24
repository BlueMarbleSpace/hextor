#
# Makefile for Darren Williams' EBM
#
FC 	= /opt/intel/bin/ifort
FLAGS	= -r8 -align all -i-static -warn nogeneral
TARGET	= driver
OBJECTS = driver.o
WDIR	= /models/williams_ebm/model

$(TARGET) : $(OBJECTS)
	$(FC) $(FLAGS) -o $(TARGET) $(OBJECTS)

driver.o : $(WDIR)/driver.f
	$(FC) $(FLAGS) -c $(WDIR)/driver.f

clean:
	rm -rf $(OBJECTS) $(TARGET)

new:
	$(MAKE) clean
	$(MAKE)
