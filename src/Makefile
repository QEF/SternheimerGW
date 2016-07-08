# Makefile for GW. 

ESPRESSO=../..
include $(ESPRESSO)/make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)$(ESPRESSO)/iotk/src $(MOD_FLAG)$(ESPRESSO)/Modules $(MOD_FLAG)$(ESPRESSO)/LAXlib \
          $(MOD_FLAG)$(ESPRESSO)/FFTXlib $(MOD_FLAG)$(ESPRESSO)/PW/src $(MOD_FLAG) $(ESPRESSO)/LR_Modules $(MOD_FLAG).

IFLAGS=

GWOBJS = \
allocate_gwq.o\
bcast_gw_input.o\
check_initial_status.o\
cBiCG_solve.o\
cBiCG_solve_coul.o\
cBiCG_solve_green.o\
ch_psi_all.o\
ch_psi_all_green.o\
clean_pw_gw.o\
close_gwq.o\
construct_w.o\
coul_multi.o\
coul_multishift.o\
coulomb.o\
coulomb_q0G0.o\
coulpade.o\
deallocate_gwq.o\
dgradcorr.o\
do_stern.o\
dv_of_drho.o\
dvqpsi_us.o\
fft6.o\
fft6_c.o\
fft6_g.o\
find_qg_ibz.o\
find_trev.o\
find_xk_uniq.o\
freqbins.o\
gauleg_grid.o\
gk_sort_safe.o\
gmap_sym.o\
godby_needs.o\
green_linsys_shift_im.o\
green_linsys_shift_re.o\
green_multishift.o\
green_multishift_im.o\
gw.o\
gw_type.o\
gwcom.o\
gwq_init.o\
gwq_readin.o\
gwq_setup.o\
h_psi_all.o\
incdrhoscf_w.o\
initialize_gw.o\
invert_epsilon.o\
kind_gw.o\
mix_pot.o\
mix_pot_c.o\
mod_diel.o\
openfilq.o\
opengwfil.o\
pade.o\
para_img.o\
para_pool.o\
pp_output.o\
prepare_q.o\
prepare_q0.o\
print_clock_gw.o\
print_matel.o\
print_matel_im.o\
q_points.o\
reorder.o\
rotate.o\
rotate_w.o\
run_nscf.o\
save_gw_input.o\
setup_nscf_green.o\
sgw_opening.o\
sigcim.o\
sigcre.o\
sigprod.o\
sigma_exch.o\
sigma_exchG.o\
sigma_expect.o\
sigma_grids.o\
sigma_matel.o\
sigma_pade.o\
solve_direct.o\
solve_linter.o\
star_q.o\
stern_symm.o\
stop_gw.o\
testnan.o\
truncate_2d.o\
truncation.o\
unfold_w.o\
version.o\
wavecut.o

LRMODS = $(ESPRESSO)/LR_Modules/liblrmod.a
PWOBJS = $(ESPRESSO)/PW/src/libpw.a
QEMODS = $(ESPRESSO)/Modules/libqemod.a $(ESPRESSO)/FFTXlib/libqefft.a $(ESPRESSO)/LAXlib/libqela.a
LIBOBJS	= $(ESPRESSO)/clib/clib.a $(ESPRESSO)/iotk/src/libiotk.a

TLDEPS= bindir mods libs pw-lib lr-lib

all : version libgw.a gw.x

version : 
	bash update_version

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

# DO NOT DELETE
include make.depend
