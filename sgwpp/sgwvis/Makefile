# Makefile for dielectric postprocessing (PP)
include ../../../make.sys
# location of needed modules
MODFLAGS= $(MOD_FLAG)../../../iotk/src $(MOD_FLAG)../../../Modules \
          $(MOD_FLAG)../../../PW/src $(MOD_FLAG).

SIGVISOBJS = \
sigmapp.o\
sigmadens.o

PWOBJS = ../../../PW/src/libpw.a
QEOBJS = ../../../Modules/libqemod.a
PPOBJS = ../../../PP/src/libpp.a

MODULES = $(PWOBJS) $(QEOBJS) $(PPOBJS)

LIBOBJS = ../../../flib/ptools.a ../../../flib/flib.a ../../../clib/clib.a ../../../iotk/src/libiotk.a 

all: sigvis.x

libsigvis.a : $(SIGVISOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

sigvis.x : sigmapp.o libsigvis.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -g -traceback -o sigma_vis.x sigmapp.o libsigvis.a \
	$(MODULES) $(LIBOBJS) $(LIBS)

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.L
