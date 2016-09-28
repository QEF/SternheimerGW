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
