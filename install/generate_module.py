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

def goto_directory_of_module_to_get_string(struct):
  string = ''
  for key in struct:
    os.chdir(key)
    string += generate_module_string(struct[key], key)
    os.chdir('..')
  return string

def generate_module_string(struct, key):
  string = key.upper() + '_MOD=' + os.getcwd() 
  if module.depend in struct:
    string += '/src\n' 
  else:
    string += '/module\n'
    string += goto_directory_of_module_to_get_string(struct)
  return string 

string = 'ESPRESSO=' + os.getcwd() + '/..\n'
string += 'VENDOR_MOD=' + os.getcwd() + '/vendor/module\n'
string += goto_directory_of_module_to_get_string(module.structure)
module_file = open('module', 'w')
module_file.write(string)
module_file.close()
