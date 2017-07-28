#------------------------------------------------------------------------------
#
# This file is part of the SternheimerGW code.
# 
# Copyright (C) 2010 - 2017
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

ESPRESSO=../..
include $(ESPRESSO)/make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)$(ESPRESSO)/iotk/src   $(MOD_FLAG)$(ESPRESSO)/UtilXlib \
          $(MOD_FLAG)$(ESPRESSO)/Modules    $(MOD_FLAG)$(ESPRESSO)/LAXlib \
          $(MOD_FLAG)$(ESPRESSO)/FFTXlib    $(MOD_FLAG)$(ESPRESSO)/PW/src \
          $(MOD_FLAG)$(ESPRESSO)/LR_Modules $(MOD_FLAG).
EXE = gw.x

IFLAGS=

GWOBJS = \
allocate_gwq.o\
bicgstab.o\
check_initial_status.o\
clean_pw_gw.o\
close_gwq.o\
construct_w.o\
coulomb.o\
coulomb_q0G0.o\
coulomb_vcut.o\
coulpade.o\
deallocate_gwq.o\
debug.o\
do_stern.o\
dvqpsi_us.o\
exchange.o\
fft6.o\
freqbins.o\
gauleg_grid.o\
gmap_sym.o\
godby_needs.o\
gram_schmidt.o\
green.o\
gw.o\
gw_input.o\
gw_opening.o\
gw_type.o\
gwcom.o\
gwq_init.o\
gwq_readin.o\
gwq_setup.o\
initialize_gw.o\
invert_epsilon.o\
kind_gw.o\
linear_op.o\
linear_solver.o\
mix_pot_c.o\
norm.o\
openfilq.o\
opengwfil.o\
pade.o\
pade_robust.o\
parallel.o\
plot_coulomb.o\
pp_output.o\
prepare_q.o\
prepare_q0.o\
print_matel.o\
print_matel_im.o\
q_points.o\
reorder.o\
run_nscf.o\
save_gw_input.o\
select_solver.o\
setup_nscf.o\
setup_nscf_green.o\
sigma.o\
sigma_expect.o\
sigma_grid.o\
sigma_io.o\
sigma_matel.o\
sigma_pade.o\
sleep.o\
solve_linter.o\
star_q.o\
stern_symm.o\
stop_gw.o\
timing.o\
truncation.o\
unfold_w.o\
version.o

LRMODS = $(ESPRESSO)/LR_Modules/liblrmod.a
PWOBJS = $(ESPRESSO)/PW/src/libpw.a
QEMODS = $(ESPRESSO)/Modules/libqemod.a $(ESPRESSO)/FFTXlib/libqefft.a \
         $(ESPRESSO)/LAXlib/libqela.a $(ESPRESSO)/UtilXlib/libutil.a
LIBOBJS	= $(ESPRESSO)/clib/clib.a $(ESPRESSO)/iotk/src/libiotk.a

TLDEPS= bindir mods libs pw-lib lr-lib

no_debug : .no_debug all
.no_debug:
	-if [ -f .debug ]; then rm debug.o .debug; fi
	-touch $@

debug : DFLAGS += -D__DEBUG
debug : .debug all
.debug:
	-if [ -f .no_debug ]; then rm debug.o .no_debug; fi
	-touch $@

all : version libgw.a $(EXE)

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

$(EXE) : gw.o libgw.a $(PWOBJS) $(LRMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ $< libgw.a \
		 $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( mkdir -p ../bin ; ln -fs ../src/$(EXE) ../bin/$(EXE) )

libgw.a : $(GWOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- rm *.x *.o *.a *.mod

.PHONY: all version clean debug no_debug

# DO NOT DELETE
include make.depend
