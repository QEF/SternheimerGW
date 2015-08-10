# Makefile for postprocessing (PP)
include ../../../make.sys

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)../../../iotk/src $(MOD_FLAG)../../../Modules \
          $(MOD_FLAG)../../../PW/src $(MOD_FLAG).
IFLAGS=

SGWPPOBJS = \
pw_restart_coul.o\
mustar.o \
ktokpmq.o \
openfil_pp.o\
coulmatsym.o\
dosband.o\
coulmatstar.o\
find_qg_ibz.o\
rotate.o\
read_file_coul.o\
parallelize.o\
star_q.o\
plotmuk.o\
dosg.o\
stop_pp.o

PWOBJS = ../../../PW/src/libpw.a
QEOBJS = ../../../Modules/libqemod.a

MODULES = $(PWOBJS) $(QEOBJS)

LIBOBJS = ../../../flib/ptools.a ../../../flib/flib.a ../../../clib/clib.a ../../../iotk/src/libiotk.a ../../src/libgw.a

all: coulmat.x

libsgwppa.a : $(SGWPPOBJS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

coulmat.x : mustar.o libsgwppa.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o coulmat.x mustar.o libsgwppa.a $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../src/coulmat.x . )

clean :
	- /bin/rm -f *.x *.o *~ *.d *.mod *.i *.L libsgwppa.a

#coulmatsym.o\
#include make.depend
