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

LIB = libvendor.a

all: .analytic module $(LIB)

.analytic:
	ln -f -s -t analytic $(ESPRESSO)/make.inc
	make -C analytic/src
	touch $@

module:
	mkdir -p $@
	for file in $$(ls */src/*mod); do \
	  export libname=$$(dirname $$(dirname $$file)); \
	  mod=$$(basename $$file); \
	  ln -s $(CURDIR)/$$file; \
	  mv $$mod $@/$$libname"_"$$mod; \
	done

$(LIB):
	rm -f $@
	$(AR) -rcT $@ */src/lib*.a

clean:
	make -C analytic/src clean ; rm -f .analytic
	rm -f module/*.mod
	rm -f $(LIB)
	

.PHONY: all clean module
