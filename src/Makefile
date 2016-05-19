# Makefile for GW. 

ESPRESSO=../..
include $(ESPRESSO)/make.sys

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)$(ESPRESSO)/iotk/src $(MOD_FLAG)$(ESPRESSO)/Modules $(MOD_FLAG)$(ESPRESSO)/LAXlib \
          $(MOD_FLAG)$(ESPRESSO)/FFTXlib $(MOD_FLAG)$(ESPRESSO)/PW/src $(MOD_FLAG) $(ESPRESSO)/LR_Modules $(MOD_FLAG).

IFLAGS=

GWOBJS = \
kind_gw.o\
gw_type.o\
gwcom.o\
version.o\
pp_output.o\
gw.o\
testnan.o\
save_gw_input.o\
bcast_gw_input.o\
gwq_readin.o\
q_points.o\
check_initial_status.o\
opengwfil.o\
freqbins.o\
gauleg_grid.o\
sigma_grids.o\
initialize_gw.o\
allocate_gwq.o\
openfilq.o\
gwq_setup.o\
prepare_q.o\
prepare_q0.o\
do_stern.o\
stern_symm.o\
gmap_sym.o\
gwq_init.o\
setup_nscf_green.o\
run_nscf.o\
coulomb.o\
coulomb_q0G0.o\
solve_direct.o\
solve_linter.o\
incdrhoscf_w.o\
dgradcorr.o\
coul_multishift.o\
cBiCG_solve_green.o\
cBiCG_solve_coul.o\
cBiCG_solve.o\
construct_w.o\
rotate_w.o\
gk_sort_safe.o\
mix_pot_c.o\
mix_pot.o\
truncate_2d.o\
dv_of_drho.o\
dvqpsi_us.o\
unfold_w.o\
invert_epsilon.o\
green_linsys_shift.o\
green_linsys_shift_im.o\
green_linsys_shift_re.o\
h_psi_all.o\
green_multishift.o\
green_multishift_im.o\
ch_psi_all.o\
ch_psi_all_green.o\
find_qg_ibz.o\
godby_needs.o\
mod_diel.o\
clean_pw_gw.o\
deallocate_gwq.o\
pade.o\
sigma_exch.o\
sigma_exchG.o\
sigma_matel.o\
print_matel.o\
para_img.o\
para_pool.o\
print_matel_im.o\
sigma_pade.o\
print_clock_gw.o\
close_gwq.o\
stop_gw.o\
fft6.o\
fft6_g.o\
fft6_c.o\
sigcim.o\
sigcre.o\
rotate.o\
coulpade.o\
star_q.o\
wavecut.o\
find_xk_uniq.o\
coul_multi.o\
sigprod.o\
sgw_opening.o\
sigcim.o\
find_trev.o

LRMODS = $(ESPRESSO)/LR_Modules/liblrmod.a
PWOBJS = $(ESPRESSO)/PW/src/libpw.a
QEMODS = $(ESPRESSO)/Modules/libqemod.a $(ESPRESSO)/FFTXlib/libqefft.a $(ESPRESSO)/LAXlib/libqela.a
LIBOBJS	= $(ESPRESSO)/clib/clib.a $(ESPRESSO)/iotk/src/libiotk.a

TLDEPS= bindir mods libs pw-lib lr-lib

all : version libgw.a gw.x

version : 
	sh update_version

version.f90 : version

gw.x : gw.o libgw.a $(PWOBJS) $(LRMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ gw.o libgw.a \
		 $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../src/gw.x . )

libgw.a : $(GWOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- rm *.x *.o *.a *.mod

.PHONY: all version clean

#include make.depend
# DO NOT DELETE
