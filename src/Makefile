# Makefile for GWTC from Makefile for PH. Need to strip out everything related to dynamical matrices and the like my ions don't move.

include make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)/iotk/src $(MOD_FLAG)/Modules \
          $(MOD_FLAG)/EE $(MOD_FLAG)/PW $(MOD_FLAG).

#F90FLAGS = #-g #-traceback -fp-stack-check -check all -fpe0

GWOBJS = \
gw.o \
gwcom.o \
allocate_gwq.o \
gwq_recover.o \
gw_restart.o \
gwq_summary.o \
run_pwscf.o \
set_defaults_pw.o \
solve_linter.o \
coulomb.o \
save_gw_input.o \
write_rec.o \
prepare_q.o \
q_points.o \
print_clock_gw.o \
gwq_readin.o \
initialize_gw.o \
check_initial_status.o \
bcast_gw_input.o \
stop_gw.o \
close_gwq.o \
gw_restart.o \
check_restart_recover.o \
clean_pw_gw.o \
gwq_setup.o \
openfilq.o \
gwq_init.o \
deallocate_gwq.o \
setlocq.o \
compute_vsgga.o \
smallgq.o \
sgam_gw.o \
davcio_drho.o \
apply_dpot.o \
cft_wave.o \
orthogonalize.o \
ch_psi_all.o \
cgsolve_all.o \
incdrhoscf.o \
incdrhoscf_w.o \
incdrhoscf_nc.o \
dv_of_drho.o\
h_psiq.o \
dgradcorr.o \
setqmod.o \
dvqpsi_us.o\
addusddens.o\
addusdbec.o\
set_drhoc.o\
cBiCG_solve.o\
cBiCG_solve_fix.o\
cBiCG_solve_green.o\
cBiCG_solve_coul.o\
cch_psi_all_fix.o\
cch_psi_all_green.o\
cg_psi.o\
green_linsys.o\
green_linsys_shift.o\
green_linsys_serial.o\
green_multishift.o\
pade.o\
freqbins.o\
refold.o\
ktokpmq.o\
prepare_kmq.o\
sigma_r2g_sco.o\
sigma_r2g_ex.o\
sig_fft_g.o\
sigma_exch.o\
sigma_exchG.o\
sigma_c.o\
sigma_c_serial.o\
sigma_matel.o\
setup_nscf_green.o\
run_pwscf_green.o\
find_q_ibz.o\
gmap_sym.o\
stern_symm.o\
unfold_w.o\
mix_pot_c.o\
mix_pot_real.o\
godby_needs.o\
sigma_extra.o\
bicgstabl.o\
solve_direct.o\
invert_epsilon.o\
mod_diel.o\
coul_multishift.o\
truncate_2D.o\
fft6.o\
write_sigma.o


HOME =/home/lambert/espresso-4.2.1
PWOBJS =$(HOME)/PW/libpw.a
EEOBJS=$(HOME)/EE/libee.a
QEMODS =$(HOME)/Modules/libqemod.a

gw.x : gw.o libgw.a $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o gw.x gw.o libgw.a \
		 $(PWOBJS) $(EEOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../src/gw.x . )

q2r.x : q2r.o rigid.o $(PWOBJS)  $(MODULES) $(EEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ q2r.o rigid.o \
		 $(PWOBJS) $(EEOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../src/$@ . )

libgw.a : $(GWOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- rm *.x *.o *.a *.mod

include make.depend

