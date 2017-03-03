#------------------------------------------------------------------------------
#
# This file is part of the Sternheimer-GW code.
# 
# Copyright (C) 2010 - 2017
# Henry Lambert, Martin Schlipf, and Feliciano Giustino
#
# Sternheimer-GW is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Sternheimer-GW is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Sternheimer-GW. If not, see
# http://www.gnu.org/licenses/gpl.html .
#
#------------------------------------------------------------------------------ 
# Makefile for utility functions for unit testing

ESPRESSO=../../..
include $(ESPRESSO)/make.inc
LIB_NAME= libtest.a

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)$(ESPRESSO)/iotk/src $(MOD_FLAG)$(ESPRESSO)/Modules $(MOD_FLAG)$(ESPRESSO)/LAXlib \
          $(MOD_FLAG)$(ESPRESSO)/FFTXlib $(MOD_FLAG)$(ESPRESSO)/PW/src $(MOD_FLAG)$(ESPRESSO)/LR_Modules \
          $(MOD_FLAG)$(ESPRESSO)/SGW/src $(MOD_FLAG).

IFLAGS =

SRCS = $(wildcard *.f90)
OBJS = $(SRCS:.f90=.o)

all: $(OBJS) $(LIB_NAME)

$(LIB_NAME): $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean: local-E1-clean
 
local-E1-clean:
	$(RM) *$(OBJ_EXT) *.mod

.PHONY: all clean
