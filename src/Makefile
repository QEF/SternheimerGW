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
bicgstab.o\
check_initial_status.o\
clean_pw_gw.o\
close_gwq.o\
construct_w.o\
coulomb.o\
coulomb_q0G0.o\
coulpade.o\
deallocate_gwq.o\
do_stern.o\
dvqpsi_us.o\
exchange.o\
fft6.o\
freqbins.o\
gauleg_grid.o\
gmap_sym.o\
godby_needs.o\
green.o\
gw.o\
gw_type.o\
gwcom.o\
gwq_init.o\
gwq_readin.o\
gwq_setup.o\
initialize_gw.o\
invert_epsilon.o\
kind_gw.o\
linear_op.o\
mix_pot_c.o\
openfilq.o\
opengwfil.o\
pade.o\
parallel.o\
pp_output.o\
prepare_q.o\
prepare_q0.o\
print_matel.o\
print_matel_im.o\
q_points.o\
reorder.o\
run_nscf.o\
save_gw_input.o\
setup_nscf.o\
setup_nscf_green.o\
sgw_opening.o\
sigma.o\
sigma_expect.o\
sigma_grid.o\
sigma_io.o\
sigma_matel.o\
sigma_pade.o\
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
QEMODS = $(ESPRESSO)/Modules/libqemod.a $(ESPRESSO)/FFTXlib/libqefft.a $(ESPRESSO)/LAXlib/libqela.a
LIBOBJS	= $(ESPRESSO)/clib/clib.a $(ESPRESSO)/iotk/src/libiotk.a

TLDEPS= bindir mods libs pw-lib lr-lib

all : version libgw.a gw.x

debug : DFLAGS += -D__DEBUG
debug : all

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

.PHONY: all version clean debug

# DO NOT DELETE
include make.depend
