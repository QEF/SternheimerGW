
# x61

#F90 = ifort -- sigsev with the intel compiler 10.1.015
F90 = f95 -O3
LDFLAGS = -lfftw -llapack
CPP = cpp 
CPPFLAGS = -C -ansi

# hbar

#F90 = mpif90 
#LDFLAGS = -llapack -lblas -L/home/giustino/lib -lfftw 
#CPP = cpp
#CPPFLAGS = -P -traditional -D__PARA -D__MPI -D__LAM -D__LINUX64

# slater
#
# if the -C -ansi options are not used then the concatenation
# operator // is taken as a comment and what follows is discarded...
# 
# slater has ifort 9.0

#F90 = mpif90 
#LDFLAGS = -llapack  -lblas -L/home/giustino/code/fftw-2.1.5/lib -lfftw
#CPP = cpp -C -ansi
#CPPFLAGS = -D__PARA -D__MPI -D__LINUX64

# pople (ifort 11)
#
# -i-static -openmp are needed when using MKL libs for ifort>9
# (from Espresso)
#
# at this time the code blows up - need to fix something in the libraries

#F90 = mpif90
#LDFLAGS = -i-static -openmp -lmkl_em64t -L/home/giustino/code/fftw-2.1.5/lib -lfftw
#CPP = cpp 
#CPPFLAGS = -D__PARA -D__MPI -D__LINUX64


#########################################################################


MODULES = \
parameters.o \
modules.o 

PARA = \
parallel_include.o \
mp_global.o \
mp.o \
io_global.o \
para.o \
startup.o \
init_pool.o 

OBJS = \
allowed.o \
date_and_tim.o \
cg_psi.o \
ch_psi_all.o \
ch_psi_all_eta.o \
coulomb.o \
coulomb_q0G0.o \
cgsolve_all.o \
bcgsolve_all.o \
bcgsolve_all_fixed.o \
eigenstates.o \
eigenstates2.o \
eigenstates_all.o \
error.o \
refold.o \
rs.o \
gwhs.o \
ggen.o \
ggens.o \
h_psi.o \
h_psi_c.o \
fftw.o \
green_coeff.o \
green_fraction.o \
green_lib.o \
hpsort_eps.o \
cgdiag.o \
ktokpmq.o \
pade.o \
rw_haydock.o \
solve_linter_dyn.o \
solve_linter_dyn_nonSCF.o \
set_ndnmbr.o \
mix_pot.o \
cryst_to_car.o \
ccgdiagg.o 

.SUFFIXES : .o .f .f90 

.f90.o:
	$(CPP) $(CPPFLAGS) $< -o $*.F90 
	$(F90) $(F90FLAGS) -c $*.F90 -o $*.o

.f.o:
	$(F90) -c $<

#
# serial
#

#gwhs.x : $(MODULES) $(OBJS)
#	$(F90) -o gwhs.x $(MODULES) $(OBJS) $(LDFLAGS) 
#	- ( cd ../bin ; ln -fs ../src/gwhs.x . )

#
# parallel
#

gwhs.x : $(MODULES) $(PARA) $(OBJS)
	$(F90) -o gwhs.x $(MODULES) $(PARA) $(OBJS) $(LDFLAGS) 
	- ( cd ../bin ; ln -fs ../src/gwhs.x . )

clean: 
	rm *.o *.mod *.x *.F90
