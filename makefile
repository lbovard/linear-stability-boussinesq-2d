FC=gfortran
objects = globals.o initialise.o fft_routines.o main.o
FFTW= -I/usr/local/include -L/usr/local/lib -lfftw3 -lm

# this command links all the files
prog: $(objects)
	$(FC) $(FFTW) -o prog $(objects)

#defines the objects
globals.mod: globals.o globals.f90
	$(FC) -c globals.f90
globals.o: globals.f90
	$(FC) -c globals.f90

initialise.mod: globals.mod initialise.o initialise.f90
	$(FC) -c initialise.f90
initialise.o: globals.mod initialise.f90 
	$(FC) -c initialise.f90

fft_routines.mod: globals.mod fft_routines.f90
	$(FC) -c  $(FFTW) fft_routines.f90
fft_routines.o: globals.mod fft_routines.f90 
	$(FC) -c  $(FFTW) fft_routines.f90

main.o: globals.mod initialise.mod fft_routines.mod main.f90
	$(FC) $(FFTW) -c main.f90
clean: 
	rm -f $(objects) prog *.mod

