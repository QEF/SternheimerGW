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

include ../module
include $(ESPRESSO)/make.inc

DEPEND = algebra fft parallel timing
DEPEND_LIB = $(foreach file, $(DEPEND), $(file)/src/lib$(file).a)
LIB = libdata.a

all: depend $(LIB)

$(LIB): $(DEPEND_LIB)
	rm -f $@
	$(AR) -rcT $@ $(DEPEND_LIB)

depend:
	for file in $(DEPEND); do \
	  make -C $$file; \
	done

test:
	for file in $(DEPEND); do \
	  make -C $$file test; \
	done

clean:
	for file in $(DEPEND); do \
	  make -C $$file clean; \
	done

.PHONY: all depend test clean
