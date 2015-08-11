# Makefile for dielectric postprocessing (PP)
include ./make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)/iotk/src $(MOD_FLAG)/Modules \
          $(MOD_FLAG)/EE $(MOD_FLAG)/PW $(MOD_FLAG).

PPOBJS = \
sigma_vis.o\
extract.o \
cube.o \
openfil_pp.o

HOME =/home/lambert/espresso-4.2
PWOBJS =$(HOME)/PW/libpw.a
EEOBJS=$(HOME)/EE/libee.a
QEMODS =$(HOME)/Modules/libqemod.a


MODULES = $(PWOBJS) $(EEOBJS) $(QEOBJS) $(QEMODS)

TLDEPS= bindir mods libs pw

diel_vis.x : sigma_vis.o $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -g -traceback -o sigma_vis.x \
	$(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	-cp ./sigma_vis.x ./SiFullSig

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.L
