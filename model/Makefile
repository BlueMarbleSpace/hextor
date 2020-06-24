#
# Makefile for Stochastic EBM
#
FC 	= /opt/intel/bin/ifort
FLAGS	= -r8 -align all -i-static -warn nogeneral
TARGET	= driver
OBJECTS = driver.o
WDIR	= /models/hextor/model

$(TARGET) : $(OBJECTS)
	$(FC) $(FLAGS) -o $(TARGET) $(OBJECTS)

driver.o : $(WDIR)/driver.f
	$(FC) $(FLAGS) -c $(WDIR)/driver.f

clean:
	rm -rf $(OBJECTS) $(TARGET)

new:
	$(MAKE) clean
	$(MAKE)
