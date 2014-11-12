# Makefile for GW. 
#
include ../../make.sys

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../PW/src $(MOD_FLAG).

IFLAGS=

GWOBJS = \
gwcom.o\
gw.o\
save_gw_input.o\
bcast_gw_input.o\
gwq_readin.o\
save_gw_input.o\
q_points.o\
check_initial_status.o\
opengwfil.o\
freqbins.o\
sigma_grids.o\
fft6.o\
initialize_gw.o\
allocate_gwq.o\
openfilq.o\
gwq_setup.o\
prepare_q.o\
do_stern.o\
stern_symm.o\
gmap_sym.o\
gwq_init.o\
set_defaults_pw.o\
run_nscf.o\
coulomb.o\
solve_direct.o\
solve_linter.o\
cft_wave.o\
incdrhoscf_w.o\
dgradcorr.o\
truncate_2D.o\
coul_multishift.o\
cBiCG_solve_coul.o\
apply_dpot.o\
cBiCG_solve.o\
incdrhoscf.o\
orthogonalize.o\
cgsolve_all.o\
mix_pot_c.o\
mix_pot.o\
cg_psi.o\
ch_psi_all.o\
h_psiq.o\
dv_of_drho.o\
dvqpsi_us.o\
unfold_w.o\
invert_epsilon.o\
green_linsys_shift.o\
cBiCG_solve_green.o\
h_psi_all.o\
green_multishift.o\
ch_psi_all_green.o\
sigma_c.o\
find_q_ibz.o\
godby_needs.o\
mod_diel.o\
clean_pw_gw.o\
deallocate_gwq.o\
pade.o\
sigma_exch.o\
sigma_matel.o\
print_matel.o\
print_matel_im.o\
sigma_pade.o\
print_clock_gw.o\
close_gwq.o\
stop_gw.o

PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a
LIBOBJS	= ../../flib/ptools.a ../../flib/flib.a ../../clib/clib.a ../../iotk/src/libiotk.a
TLDEPS= bindir mods libs pw-lib

all : libgw.a gw.x

gw.x : gw.o libgw.a $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o gw.x gw.o libgw.a \
		 $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../src/gw.x . )

libgw.a : $(GWOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- rm *.x *.o *.a *.mod

#include make.depend
# DO NOT DELETE
