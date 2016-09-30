.PHONY: clean distclean tests all
.DEFAULT_GOAL = tests

TOP_DIR := $(shell pwd)
SRC_DIR = $(TOP_DIR)/../src
UTIL_DIR = $(TOP_DIR)/lib
LIB_NAME= libtest.a
TEST_DIR= $(TOP_DIR)/tests
ESPRESSO=../..

VPATH = . $(SRC_DIR) $(UTIL_DIR) $(TEST_DIR)

include $(ESPRESSO)/make.inc
include $(PFUNIT)/include/base.mk

ifeq ($(USEMPI),YES)
   MPIF90=mpif90
endif

SGWOBJS = $(UTIL_DIR)/$(LIB_NAME) $(SRC_DIR)/libgw.a
LRMODS  = $(ESPRESSO)/LR_Modules/liblrmod.a
PWOBJS  = $(ESPRESSO)/PW/src/libpw.a
QEMODS  = $(ESPRESSO)/Modules/libqemod.a $(ESPRESSO)/FFTXlib/libqefft.a $(ESPRESSO)/LAXlib/libqela.a
LIBOBJS = $(ESPRESSO)/clib/clib.a $(ESPRESSO)/iotk/src/libiotk.a
LIBS += $(SGWOBJS) $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS)
LIBS += $(SCALAPACK_LIBS) $(LAPACK_LIBS) $(FFT_LIBS) $(BLAS_LIBS) $(MPI_LIBS) $(MASS_LIBS) ${HDF5_LIB} $(LD_LIBS)


EXE = tests$(EXE_EXT)

all: $(EXE)
	mpirun -np 4 ./$(EXE)

SUT:
	make -C $(SRC_DIR)
	make -C $(UTIL_DIR)
	make -C $(TEST_DIR) tests

tests: all

$(EXE): testSuites.inc SUT
	$(MPIF90) -o $@ -I$(PFUNIT)/mod -I$(PFUNIT)/include -Itests $(PFUNIT)/include/driver.F90 $(TEST_DIR)/*$(OBJ_EXT) -L$(PFUNIT)/lib -lpfunit $(LIBS) -DUSE_MPI $(FFLAGS) $(FPPFLAGS)

echo: local-E0-echo

local-E0-echo:
	@echo PFUNIT: $(PFUNIT)
	@echo BUILDROBUST: \"$(BUILDROBUST)\"
	@echo COMPILER: $(COMPILER)
	@echo FFLAGS: $(FFLAGS)
	@echo FPPFLAGS: $(FPPFLAGS)
	@echo FC: $(FC)
	@echo MPI: $(MPI)
	@echo CPPFLAGS: $(CPPFLAGS)
	@echo LIBS: $(LIBS)

clean: local-E0-clean

local-E0-clean:
	make -C $(SRC_DIR) clean
	make -C $(UTIL_DIR) clean
	make -C $(TEST_DIR) clean
	rm -f $(EXE) *$(OBJ_EXT)

export FC
export FPPFLAGS
export FFLAGS
export SRC_DIR
export TEST_DIR
export OBJ_EXT
export LIB_EXT
export EXE_EXT
export MPIF90
