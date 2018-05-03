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

include ../../../module
include $(ESPRESSO)/make.inc

SRC = $(wildcard *.f90)
OBJ = $(SRC:.f90=.o)
LIB = libcoul.a

MODFLAGS = $(BASEMOD_FLAGS) $(MOD_FLAG)$(UTIL_MOD) $(MOD_FLAG)$(ALGO_MOD)

all: lib mod

lib: $(LIB)

mod: $(OBJ)
	for file in $$(ls *.mod); do \
	  ln -f -s $(CURDIR)/$$file -t ../../module/; \
	done

$(LIB) : $(OBJ)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean:
	-rm *.o *.mod $(LIB)

.PHONY: all lib mod

-include make.depend
