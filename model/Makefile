#
# Makefile for Darren Williams' EBM
#
FC 	= ifort
FLAGS	= -r8 -align all -i-static
TARGET	= driver
OBJECTS = driver.o
WDIR	= /neptune/s1/jhaqqmis/williams_ebm/model

$(TARGET) : $(OBJECTS)
	$(FC) $(FLAGS) -o $(TARGET) $(OBJECTS)

driver.o : $(WDIR)/driver.f
	$(FC) $(FLAGS) -c $(WDIR)/driver.f

clean:
	rm -rf $(OBJECTS) $(TARGET)

new:
	$(MAKE) clean
	$(MAKE)
