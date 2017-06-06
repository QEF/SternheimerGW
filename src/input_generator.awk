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
# initialize
BEGIN {
  num_namelist = 0
  last_whitespace = 0
  var_active = 0
  spec_active = 0
  # constants that indicate which part we work on
  k_type = 1
  k_default = 2
  k_description = 3
}
# check for comment or empty lines
/^[ \t]*#/ || NF == 0 { skip = 1 }
!/^[ \t]*#/ && NF > 0 { skip = 0 }

# lines introducing new element contain ':'
skip == 0 && /:/ {

  # if there is a :, this is not a continuation line
  cont_line = 0

  # set whitespace in variable RLENGTH
  match($0, /^ */)
  whitespace = RLENGTH

  # no whitespace - new namelist
  if (whitespace == 0) {
    num_namelist++
    split($0, array, ":")
    namelist[num_namelist] = array[1]
    # initialize variable count
    num_variable[num_namelist] = 0 
  }

  # if the whitespace increased, we go to a subelement
  # namelist -> variable -> specification
  if (whitespace > last_whitespace) {
     if (var_active == 0 && spec_active == 0) {
       # namelist -> variable
       var_active = 1
     } else if (var_active == 1 && spec_active == 0) {
       # variable -> specification
       var_active = 0
       spec_active = 1
     } else {
       print "Error: increasing whitespace in line", NR > "/dev/stderr"
       exit 1
     }
  }

  # if whitespace decreased, we go back to parent element
  # specification -> variable -> namelist
  if (whitespace < last_whitespace) {
    if (var_active == 0 && spec_active == 1) {
       # specification -> variable
       var_active = 1
       spec_active = 0
    } else if (var_active == 1 && spec_active == 0) {
       # variable -> namelist
       var_active = 0
    } else {
       print "Error: decreasing whitespace in line", NR > "/dev/stderr"
       exit 1
    }
  }

  last_whitespace = whitespace

}

skip == 0 && !/:/ { 
  # if there is a :, this is not a continuation line
  cont_line = 1
  # continuation lines are only allowed in specification part
  if (spec_active == 0) {
    print "Error: line", NR, "appears to continue the previous line" > "/dev/stderr"
    print "       This is only allowed for the specification part!"  > "/dev/stderr"
    exit 1
  }
}

# store the name of the variable
skip == 0 && var_active {
  num_variable[num_namelist]++
  # introduce abbreviation
  num_var = num_variable[num_namelist]
  split($1, array, ":")
  variable[num_namelist,num_var] = array[1]
}

# store the elements of the variable
skip == 0 && spec_active && cont_line == 0 {

  # determine which element we deal with
  split($0, array, ":")
  gsub(/^[ \t]/, "", array[1])
  element = tolower(array[1])

  if (element == "type") {
    active = k_type 
    type[num_namelist,num_var] = array[2]
  } else if (element == "default") {
    active = k_default
    default[num_namelist,num_var] = array[2]
  } else if (element == "description") {
    active = k_description
    description[num_namelist,num_var] = array[2]
  } else {
    print "Error: element in line", NR, "does not match any know element" > "/dev/stderr"
    exit 1
  }

}

# continuation line
cont_line == 1 {

  # trim whitespace
  gsub(/^[ \t]/, "", $0)

  # add content of line to active element
  if (active == k_type) {
    type[num_namelist,num_var] = type[num_namelist,num_var]" "$0
  } else if (active == k_default) {
    default[num_namelist,num_var] = default[num_namelist,num_var]" "$0
  } else if (active == k_description) {
    description[num_namelist,num_var] = description[num_namelist,num_var]" "$0
  }

  cont_line = 0
}

# generate the reading routine
END {

  print "!------------------------------------------------------------------------------"
  print "!"
  print "! This file is part of the SternheimerGW code."
  print "!"
  print "! Copyright (C) 2010 - 2017"
  print "! Henry Lambert, Martin Schlipf, and Feliciano Giustino"
  print "!"
  print "! SternheimerGW is free software: you can redistribute it and/or modify"
  print "! it under the terms of the GNU General Public License as published by"
  print "! the Free Software Foundation, either version 3 of the License, or"
  print "! (at your option) any later version."
  print "!"
  print "! SternheimerGW is distributed in the hope that it will be useful,"
  print "! but WITHOUT ANY WARRANTY; without even the implied warranty of"
  print "! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the"
  print "! GNU General Public License for more details."
  print "!"
  print "! You should have received a copy of the GNU General Public License"
  print "! along with SternheimerGW. If not, see"
  print "! http://www.gnu.org/licenses/gpl.html ."
  print "!"
  print "!------------------------------------------------------------------------------"
  print "! NOTE: This file is automatically generated - do NOT modify."
  print "! Modify the definition file input_gw.yml instead."
  print "!------------------------------------------------------------------------------"
  print "!> Automatically generated input reader for SternheimerGW"
  print "MODULE gw_input_module"
  print ""
  print "  IMPLICIT NONE"
  print ""
  print "  PRIVATE"
  print "  PUBLIC gw_input_type, gw_input"
  print ""
  print "  !> contains all data read from input"
  print "  TYPE gw_input_type"
  print "  END TYPE gw_input_type"
  print ""
  print "CONTAINS"
  print ""
  print "  !> read the user input from the input file"
  print "  SUBROUTINE gw_input(user_input)"
  print ""
  print "    !> store the user input in this type"
  print "    TYPE(gw_input_type), INTENT(OUT) :: user_input"
  print ""
  print "  END SUBROUTINE gw_input"
  print ""
  print "END MODULE gw_input_module"
#  for (nml = 1; nml <= num_namelist; nml++) {
#    print namelist[nml]
#    for (var = 1; var <= num_variable[nml]; var++) {
#       print " ", variable[nml,var]
#       print "   ", type[nml,var]
#       print "   ", default[nml,var]
#       print "   ", description[nml,var]
#    }
#  }
}
