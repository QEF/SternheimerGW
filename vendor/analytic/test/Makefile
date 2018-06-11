# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

.PHONY: all lib
.DEFAULT_GOAL = all

include $(PFUNIT)/include/base.mk

FFLAGS += -I$(PFUNIT)/mod -I../src

ifeq ($(USEMPI),YES)
   MPIF90 = mpif90
   FPPFLAGS += -DUSE_MPI
endif

EXE = tests$(EXE_EXT)
SRC = $(wildcard *.pf)
OBJ = $(SRC:.pf=$(OBJ_EXT))

LIB_LAPACK = -llapack -lblas

LIB = ../src/libanalytic.a

all: lib $(EXE)
	mpirun -np 4 ./$(EXE)

$(EXE): testSuites.inc $(OBJ) $(LIB)
	$(MPIF90) -o $@ -I. $(FFLAGS) $(FPPFLAGS) $(PFUNIT)/include/driver.F90 $(OBJ) -L$(PFUNIT)/lib -lpfunit $(LIBS) $(LIB) $(LIB_LAPACK)

%$(OBJ_EXT): %.F90
	$(MPIF90) -c $(FFLAGS) $<

%.F90: %.pf
	$(PFUNIT)/bin/pFUnitParser.py $<  $@

assert.x: testAssert.f90 $(LIB)
	$(MPIF90) -o $@ $< $(FFLAGS) $(FPPFLAGS) $(LIBS) $(LIB)

lib:
	make -C ../src lib
