# Makefile for postprocessing (PP)

ESPRESSO=../../..
SGW=../..
include $(ESPRESSO)/make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)$(ESPRESSO)/iotk/src $(MOD_FLAG)$(ESPRESSO)/Modules \
          $(MOD_FLAG)$(ESPRESSO)/LAXlib  $(MOD_FLAG)$(ESPRESSO)/FFTXlib \
          $(MOD_FLAG)$(ESPRESSO)/PW/src $(MOD_FLAG)$(ESPRESSO)/LR_Modules \
          $(MOD_FLAG).

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


LRMODS = $(ESPRESSO)/LR_Modules/liblrmod.a
PWOBJS = $(ESPRESSO)/PW/src/libpw.a
QEMODS = $(ESPRESSO)/Modules/libqemod.a $(ESPRESSO)/FFTXlib/libqefft.a \
         $(ESPRESSO)/LAXlib/libqela.a

MODULES = $(PWOBJS) $(QEMODS) $(LRMODS)

LIBOBJS = $(ESPRESSO)/clib/clib.a $(ESPRESSO)/iotk/src/libiotk.a $(SGW)/src/libgw.a

all: coulmat.x

libsgwppa.a : $(SGWPPOBJS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

coulmat.x : mustar.o libsgwppa.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o coulmat.x mustar.o libsgwppa.a $(MODULES) $(LIBOBJS) $(LIBS)
	- ( ln -fs $$(pwd)/coulmat.x $(SGW)/bin/coulmat.x )

clean :
	- /bin/rm -f *.x *.o *~ *.d *.mod *.i *.L libsgwppa.a

#coulmatsym.o\
#include make.depend
