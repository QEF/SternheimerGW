# Makefile for postprocessing (PP)
include ../../../make.sys

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)../../../iotk/src $(MOD_FLAG)../../../Modules \
          $(MOD_FLAG)../../../PW/src $(MOD_FLAG).
IFLAGS=

SGWPPOBJS = \
coulmat.o \
ktokpmq.o \
openfil_pp.o\
dosg.o\
stop_pp.o

PWOBJS = ../../../PW/src/libpw.a
QEOBJS = ../../../Modules/libqemod.a

MODULES = $(PWOBJS) $(QEOBJS)

LIBOBJS = ../../../flib/ptools.a ../../../flib/flib.a ../../../clib/clib.a ../../../iotk/src/libiotk.a

all: coulmat.x

libsgwppa.a : $(SGWPPOBJS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

coulmat.x : coulmat.o libsgwppa.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o coulmat.x coulmat.o libsgwppa.a $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../src/coulmat.x . )

clean :
	- /bin/rm -f *.x *.o *~ *.d *.mod *.i *.L libsgwppa.a

#include make.depend
