#!/bin/bash
#------------------------------------------------------------------------------
#
# This file is part of the SternheimerGW code.
# Parts of this file are taken from the Quantum ESPRESSO software
# P. Giannozzi, et al, J. Phys.: Condens. Matter, 21, 395502 (2009)
#
# Copyright (C) 2010 - 2018 Quantum ESPRESSO group,
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

source ./ENVIRONMENT

if test "`which curl`" = "" ; then
   if test "`which wget`" = "" ; then
      echo "### wget or curl not found: will not be able to download missing PP ###"
   else
      DOWNLOADER="wget -O"
      # echo "wget found"
   fi
else
   DOWNLOADER="curl -o"
   # echo "curl found"
fi

inputs=`find $1* -type f -name "*.in" -not -name "test.*" -not -name "benchmark.*"`
pp_files=`for x in ${inputs}; do grep UPF ${x} | awk '{print $3}'; done`

for pp_file in ${pp_files} ; do
if ! test -f ${ESPRESSO_PSEUDO}/${pp_file} ; then 
	#echo -n "Downloading ${pp_file} to ${ESPRESSO_PSEUDO} ... "
	${DOWNLOADER} ${ESPRESSO_PSEUDO}/${pp_file} ${NETWORK_PSEUDO}/${pp_file} 2> /dev/null
	if test $? != 0 ; then 
		echo "FAILED, do it manually -- Testing aborted!"
		exit -1
	#else
		#echo "done."
	fi
#else
	#echo "No need to download ${pp_file}."
fi
done
