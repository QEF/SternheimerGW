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

main_makefile = open('main_makefile', 'r').read()
root_makefile = open('root_makefile', 'r').read()
src_makefile = open('src_makefile', 'r').read()

os.chdir('..')

def goto_directory_to_create_makefile(struct, layer):
  for key in struct:
    os.chdir(key)
    create_makefile_or_recurse(struct[key], key, layer)
    os.chdir('..')

def create_makefile_or_recurse(struct, key, layer):
    if module.depend in struct:
      create_makefile(struct[module.depend], key, layer) 
    else:
      create_makefile_root(struct, key)
      goto_directory_to_create_makefile(struct, layer + 1)

def create_makefile(depend_array, key, layer):
  create_makefile_main()
  if layer > 1: create_makefile_src(depend_array, key)

def create_makefile_main():
  makefile = open('Makefile~', 'w') 
  makefile.write(main_makefile)
  makefile.close()

def create_makefile_src(depend_array, key):
  makefile = open('src/Makefile~', 'w')
  src_string = generate_src_string(depend_array, key)
  makefile.write(src_string)
  makefile.close()

def generate_src_string(depend_array, key):
  result = src_makefile
  result = result.replace('@LIB@', key)
  result = result.replace('@MOD@', generate_mod_string(depend_array))
  return result

def generate_mod_string(depend_array):
  result = ''
  for dep in depend_array:
    result += ' ' + mod_string(dep)
  return result

def mod_string(dep):
  if dep == 'base':
    return '$(BASEMOD_FLAGS)'
  elif dep == 'pw':
    return '$(MOD_FLAG)$(ESPRESSO)/PW/src'
  elif dep == 'lrmods':
    return '$(MOD_FLAG)$(ESPRESSO)/LR_Modules'
  else:
    return '$(MOD_FLAG)$(' + dep.upper() + '_MOD)'

def create_makefile_root(struct, key):
  makefile = open('Makefile~', 'w')
  root_string = generate_root_string(struct, key)
  makefile.write(root_string)
  makefile.close()

def generate_root_string(struct, key):
  result = root_makefile
  result = result.replace('@DEP@', generate_depend_string(struct))
  result = result.replace('@LIB@', key)
  return result

def generate_depend_string(struct):
  result = ''
  for key in struct:
    result += ' ' + key
  return result

goto_directory_to_create_makefile(module.structure, 1)
