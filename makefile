FC:=gfortran
objects = globals.o initialise.o main.o

# this command links all the files
prog: $(objects)
	$(FC) -o prog $(objects)

#defines the objects
globals.mod: globals.o globals.f90
	$(FC) -c globals.f90
globals.o: globals.f90
	$(FC) -c globals.f90
initialise.mod: globals.mod initialise.o initialise.f90
	$(FC) -c initialise.f90
initialise.o: globals.mod initialise.f90 
	$(FC) -c initialise.f90
main.o: globals.mod initialise.mod main.f90
	$(FC) -c main.f90
clean: 
	rm $(objects) prog *.mod

