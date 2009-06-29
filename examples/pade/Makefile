
F90 = f95 -O3

.SUFFIXES : .o .f .f90 

.f90.o:
	$(F90) $(F90FLAGS) -c $*.f90 -o $*.o

pade.x : testpade.o pade.o
	$(F90) -o pade.x testpade.o pade.o

clean: 
	rm *.o *.mod *.x *.F90
