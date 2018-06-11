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
all: install
	make -C util
	make -C data
	make -C algo
	make -C phys
	make -C main

install:
	make -C install
	make -C vendor

depend: install
	make -C util depend
	make -C data depend
	make -C algo depend
	make -C phys depend
	make -C main depend

test: install
	make -C util test
	make -C data test
	make -C algo test
	make -C phys test
	make -C main test

clean: install
	make -C install clean
	make -C vendor clean
	make -C util clean
	make -C data clean
	make -C algo clean
	make -C phys clean
	make -C main clean

.PHONY: all install test clean
