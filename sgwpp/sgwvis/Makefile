# Makefile for dielectric postprocessing (PP)
include ../../../make.sys

MODFLAGS= $(MOD_FLAG)../../../iotk/src $(MOD_FLAG)../../../Modules \
          $(MOD_FLAG)../../../PW/src $(MOD_FLAG).

IFLAGS=



SIGVISOBJS = \
sigmapp.o\
sigmadens.o

PWOBJS = ../../../PW/src/libpw.a
QEOBJS = ../../../Modules/libqemod.a
PPOBJS = ../../../PP/src/libpp.a
SGWOBJS = ../../src/libgw.a

MODULES = $(PWOBJS) $(QEOBJS) $(PPOBJS) $(SGWOBJS)

LIBOBJS = ../../../flib/ptools.a ../../../flib/flib.a ../../../clib/clib.a ../../../iotk/src/libiotk.a 

all: libsigvis.a sigvis.x

libsigvis.a : $(SIGVISOBJS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

sigvis.x : sigmapp.o libsigvis.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o sigma_vis.x sigmapp.o libsigvis.a \
	$(MODULES) $(LIBOBJS) $(LIBS)

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.L
