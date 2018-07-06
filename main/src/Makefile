#------------------------------------------------------------------------------
#
# This file is part of the SternheimerGW code.
# 
# Copyright (C) 2010 - 2018
# Henry Lambert, Martin Schlipf, and Feliciano Giustino
#
# SternheimerGW is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SternheimerGW is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SternheimerGW. If not, see
# http://www.gnu.org/licenses/gpl.html .
#
#------------------------------------------------------------------------------ 
# Makefile for GW. 

include ../../module
include $(ESPRESSO)/make.inc

MODFLAGS = $(BASEMOD_FLAGS) $(MOD_FLAG)$(ESPRESSO)/PW/src $(MOD_FLAG)$(ESPRESSO)/LR_Modules \
  $(MOD_FLAG)$(VENDOR_MOD) $(MOD_FLAG)$(UTIL_MOD) $(MOD_FLAG)$(DATA_MOD) $(MOD_FLAG)$(ALGO_MOD) $(MOD_FLAG)$(PHYS_MOD)

EXE = gw.x

GWOBJS = \
gw.o\
gw_input.o\
gw_opening.o\
gwq_readin.o\
version.o

GWLIBS = ../../phys/libphys.a ../../algo/libalgo.a ../../data/libdata.a ../../util/src/libutil.a ../../vendor/libvendor.a
LRMODS = $(ESPRESSO)/LR_Modules/liblrmod.a
PWOBJS = $(ESPRESSO)/PW/src/libpw.a
QEMODS = $(ESPRESSO)/Modules/libqemod.a             $(ESPRESSO)/FFTXlib/libqefft.a \
         $(ESPRESSO)/KS_Solvers/Davidson/libdavid.a $(ESPRESSO)/KS_Solvers/CG/libcg.a \
         $(ESPRESSO)/LAXlib/libqela.a               $(ESPRESSO)/UtilXlib/libutil.a \
         $(ESPRESSO)/dft-d3/libdftd3qe.a
LIBOBJS	= $(ESPRESSO)/clib/clib.a $(ESPRESSO)/iotk/src/libiotk.a

all : version libgw.a $(EXE)

depend: tmp.depend
	awk '!/@/' $< > make.depend

tmp.depend: $(SRC)
	$(ESPRESSO)/install/moduledep.sh > $@
	$(ESPRESSO)/install/includedep.sh >> $@

version : 
	bash update_version

version.f90 : version

gw_input.f90 : gw_input.yml input_generator.awk
	awk -f input_generator.awk $< > $@

tiddler: gw_input.yml input_generator.awk check-wiki
	awk -f input_generator.awk -v tiddler=1 -v path=$(WIKI)/tiddlers $<

check-wiki:
ifndef WIKI
	$(error WIKI must be set to the path to the wiki.)
endif
ifeq ("$(wildcard $(WIKI)/tiddlers)", "")
	$(error directory $(WIKI)/tiddlers does not exist)
endif

$(EXE) : gw.o libgw.a $(GWLIBS) $(PWOBJS) $(LRMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ $< libgw.a $(GWLIBS) \
		 $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( mkdir -p ../bin ; ln -fs ../src/$(EXE) ../bin/$(EXE) )

libgw.a : $(GWOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- rm -f *.x *.o *.a *.mod

.PHONY: all version clean
.INTERMEDIATE: tmp.depend

-include make.depend
