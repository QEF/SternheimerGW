# Makefile for dielectric postprocessing (PP)
#
ESPRESSO=../../..
include $(ESPRESSO)/make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)$(ESPRESSO)/iotk/src $(MOD_FLAG)$(ESPRESSO)/Modules \
          $(MOD_FLAG)$(ESPRESSO)/PW/src $(MOD_FLAG).

IFLAGS=

SIGVISOBJS = \
sigmapp.o\
sigmadens.o

PWOBJS = $(ESPRESSO)/PW/src/libpw.a
QEMODS = $(ESPRESSO)/Modules/libqemod.a
LIBOBJS	= $(ESPRESSO)/flib/ptools.a $(ESPRESSO)/flib/flib.a $(ESPRESSO)/clib/clib.a $(ESPRESSO)/iotk/src/libiotk.a
TLDEPS= bindir mods libs pw-lib

all: libsigvis.a sigvis.x

sigvis.x : sigmapp.o libsigvis.a $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -g -traceback -o sigma_vis.x sigmapp.o libsigvis.a \
		$(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)

libsigvis.a : $(SIGVISOBJS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@


clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.L
