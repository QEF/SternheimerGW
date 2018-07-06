# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

.PHONY: lib

include ../make.inc

SRC = $(wildcard *.f90)
OBJ = $(SRC:.f90=.o)
LIB = libanalytic.a

lib: $(LIB)

$(LIB) : $(OBJ)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean:
	-rm *.o *.mod $(LIB)

-include make.depend
