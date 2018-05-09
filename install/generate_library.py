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
import os
import module 

os.chdir('..')

def goto_directory_of_library_to_get_string(struct):
  string = ''
  for key in struct:
    os.chdir(key)
    string += generate_library_string(struct[key], key)
    os.chdir('..')
  return string

def generate_library_string(struct, key):
  string = key.upper() + '_LIB = ' + os.getcwd() 
  if module.depend in struct:
    string += '/src/lib' + key.lower() + '.a\n' 
  else:
    string += '/lib' + key.lower() + '.a\n'
    string += goto_directory_of_library_to_get_string(struct)
  return string 

string  = 'BASE_LIB = $(ESPRESSO)/Modules/libqemod.a             $(ESPRESSO)/FFTXlib/libqefft.a \\\n'
string += '           $(ESPRESSO)/KS_Solvers/Davidson/libdavid.a $(ESPRESSO)/KS_Solvers/CG/libcg.a \\\n'
string += '           $(ESPRESSO)/LAXlib/libqela.a               $(ESPRESSO)/UtilXlib/libutil.a \\\n'
string += '           $(ESPRESSO)/dft-d3/libdftd3qe.a $(ESPRESSO)/clib/clib.a $(ESPRESSO)/iotk/src/libiotk.a\n'
string += 'PW_LIB = $(ESPRESSO)/PW/src/libpw.a\n'
string += 'LR_LIB = $(ESPRESSO)/LR_Modules/liblrmod.a\n'
string += goto_directory_of_library_to_get_string(module.structure)
library_file = open('library', 'w')
library_file.write(string)
library_file.close()
