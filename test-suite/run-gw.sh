#!/bin/bash
#------------------------------------------------------------------------------
#
# This file is part of the SternheimerGW code.
# 
# Copyright (C) 2010 - 2017 
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

bash ../ENVIRONMENT

if [[ $QE_USE_MPI == 1 ]]; then
  export PARA_PREFIX="mpirun -np ${TESTCODE_NPROCS}"
  export PARA_SUFFIX="-npool ${TESTCODE_NPROCS}"
else
  unset PARA_PREFIX
  unset PARA_SUFFIX
fi

# determine whether we are in a PW or SternheimerGW step
scf=$(echo $1 | awk '/scf.in/{ print 1 }')
gw=$(echo $1 | awk '/gw.in/{ print 1 }')

if (( scf == 1 ))
then
  echo "Running PW..."
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x ${PARA_SUFFIX} < $1 > $2 2> $3
elif (( gw == 1 ))
then
  echo "Running SternheimerGW..."
  ${PARA_PREFIX} ${SternheimerGW_ROOT}/bin/gw.x ${PARA_SUFFIX} < $1 > $2 2> $3
else
  echo "Unknown input file" > /dev/stderr
fi  
