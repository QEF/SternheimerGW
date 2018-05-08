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
test_makefile = open('test_makefile', 'r').read()

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
  if layer > 1:
    create_makefile_src(depend_array, key)
    create_makefile_test_if_possible(depend_array, key)

def create_makefile_main():
  makefile = open('Makefile', 'w') 
  makefile.write(main_makefile)
  makefile.close()

def create_makefile_src(depend_array, key):
  makefile = open('src/Makefile', 'w')
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

def create_makefile_test_if_possible(depend_array, key):
  if os.path.isdir('test'):
    create_makefile_test(depend_array, key)

def create_makefile_test(depend_array, key):
  makefile = open('test/Makefile', 'w')
  test_string = generate_test_string(depend_array, key)
  makefile.write(test_string)
  makefile.close()

def generate_test_string(depend_array, key):
  result = test_makefile
  result = result.replace('@LIB@', key)
  result = result.replace('@MOD@', generate_mod_string(depend_array))
  result = result.replace('@LINK@', generate_link_string(depend_array))
  return result

def generate_link_string(depend_array):
  result = ''
  for dep in depend_array:
    result += ' ' + link_string(dep)
  return result

def link_string(dep):
  if dep == 'base':
    return '$(BASELIB)'
  elif dep == 'pw':
    return '$(ESPRESSO)/PW/src/libpw.a'
  elif dep == 'lrmods':
    return '$(ESPRESSO)/LR_Modules/liblrmod.a'
  else:
    return ''

def create_makefile_root(struct, key):
  makefile = open('Makefile', 'w')
  root_string = generate_root_string(struct, key)
  makefile.write(root_string)
  makefile.close()

def generate_root_string(struct, key):
  result = root_makefile
  result = result.replace('@MOD@', generate_depend_string(struct))
  result = result.replace('@LIB@', key)
  return result

def generate_depend_string(struct):
  result = ''
  for key in module.sort_dict_by_order(struct):
    result += ' ' + key
  return result

goto_directory_to_create_makefile(module.structure, 1)
