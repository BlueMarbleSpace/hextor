#
# Makefile for testing radiation module
#

### Edit FC and WDIR to match your system ###
FC 	= /usr/bin/h5pfc
WDIR	= /models/hextor/model/radiation
#############################################

FLAGS	= -r8 -align all -i-static -warn nogeneral
TARGET	= radtest
OBJECTS = radiation.o main.o

$(TARGET) : $(OBJECTS)
	$(FC) $(FLAGS) -o $(TARGET) $(OBJECTS)


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
