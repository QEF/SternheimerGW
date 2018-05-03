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
all: module 
	make -C util
	make -C data
	make -C algo
	make -C phys
	make -C user 

test:
	make -C util test
	make -C data test
	make -C algo test
	make -C phys test
	make -C user test 

clean:
	make -C util clean
	make -C data clean
	make -C algo clean
	make -C phys clean
	make -C user clean

module: Makefile
	echo "ESPRESSO=$(CURDIR)/.." > $@
	echo "UTIL_MOD=$(CURDIR)/util/src" >> $@
	echo "DATA_MOD=$(CURDIR)/data/module" >> $@
	echo "ALGO_MOD=$(CURDIR)/algo/module" >> $@
	echo "PHYS_MOD=$(CURDIR)/phys/module" >> $@

.PHONY: all test clean
