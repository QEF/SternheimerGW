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
  # check if input should be generated
  if (tiddler == 1) {
    gen_input = 0
  } else {
    gen_input = 1
  }
  # character identifying a list
  list_char = "-"
}
function line_cont(first_field) {
  if (gen_input == 1) {
    return "\n    !!"
  }
  if (tiddler == 1) {
    if (first_field == list_char) {
      return "\n"
    } else {
      return ""
    }
  }
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
    var_active = 0
    spec_active = 0
  }
  # if the whitespace increased, we go to a subelement
  # namelist -> variable -> specification
  else if (whitespace > last_whitespace) {
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
  else if (whitespace < last_whitespace) {
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
  sub(/^[ \t]*/, "", array[1])
  sub(/^[ \t]*/, "", array[2])
  element = tolower(array[1])

  if (element == "type") {
    active = k_type 
    type[num_namelist,num_var] = array[2]
  } else if (element == "default") {
    active = k_default
    default_[num_namelist,num_var] = array[2]
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

  # list start with a -
  list = ($1 == list_char)

  # at the start of a list print a new line
  new_list = (list && !prev_list)

  # trim whitespace
  sub(/^[ \t]*/, "", $0)

  # add content of line to active element
  if (active == k_type) {
    type[num_namelist,num_var] = type[num_namelist,num_var]" "$0
  } else if (active == k_default) {
    default_[num_namelist,num_var] = default_[num_namelist,num_var]" "$0
  } else if (active == k_description) {
    if (new_list) description[num_namelist,num_var] = description[num_namelist,num_var]line_cont(list_char)
    # tiddlywiki uses different list character
    current_line = $0
    if (list && tiddler == 1) sub(list_char, "*", current_line)
    description[num_namelist,num_var] = description[num_namelist,num_var]line_cont($1)" "current_line
  }

  cont_line = 0
  prev_list = list

}

# print input parameters for subroutine
function print_interface(name) {
  printf "  SUBROUTINE "name"("
  for (nml = 1; nml <= num_namelist; nml++) {
    if (nml > 1) printf ", "
    printf namelist[nml]"_t"
  }
  print ")"
}
 
# print the namelist type
function print_type(intent) {
  for (nml = 1; nml <= num_namelist; nml++) {
    print ""
    print "    !> contains the user input in namelist", namelist[nml]
    print "    TYPE("namelist[nml]"_type), INTENT("intent") :: "namelist[nml]"_t"
  }
  print ""
}

# generate the reading routine
END {

  if (gen_input == 1) {
    print "!------------------------------------------------------------------------------"
    print "!"
    print "! This file is part of the SternheimerGW code."
    print "!"
    print "! Copyright (C) 2010 - 2018"
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
    print "  USE kinds, ONLY: DP"
    # for any private type, load the private bcast routine
    for (nml = 1; nml <= num_namelist; nml++) {
      for (var = 1; var <= num_variable[nml]; var++) {
        if (index(toupper(type[nml, var]), "TYPE") != 0) {
          print "  USE", variable[nml, var]"_module, ONLY: "variable[nml, var]"_type"
        }
      }
    }
    print ""
    print "  IMPLICIT NONE"
    # generate one type per input namelist
    for (nml = 1; nml <= num_namelist; nml++) {
      print ""
      print "  !> contains all data read from", namelist[nml]
      print "  TYPE", namelist[nml]"_type"
      for (var = 1; var <= num_variable[nml]; var++) {
        print ""
        print "    !>", description[nml, var]
        print "   ", toupper(type[nml, var]), "::", variable[nml, var]
      }
      print ""
      print "  END TYPE", namelist[nml]"_type"
    }
    print ""
    print "CONTAINS"
    print ""
    print "  !> read the user input from the input file"
    # return one element per input namelist
    print_interface("gw_input_read")
    print ""
    print "    USE io_global, ONLY: stdin"
    # definition of the return types
    print_type("OUT")
    # private types that are used to read the input
    for (nml = 1; nml <= num_namelist; nml++) {
      if (nml > 1) print ""
      print "    !"
      print "    !  variables used for namelist", namelist[nml]
      for (var = 1; var <= num_variable[nml]; var++) {
        print "    !"
        print "    !>", description[nml, var]
        print "   ", toupper(type[nml, var]), "::", variable[nml, var]
      }
    }
    print ""
    print "    !> error flag returned by reading the namelist"
    print "    INTEGER :: ierr"
    for (nml = 1; nml <= num_namelist; nml++) {
      # skip empty namelist
      if (num_variable[nml] == 0) continue
      print ""
      print "    NAMELIST /"toupper(namelist[nml])"/ &"
      old_string = ""
      for (var = 1; var < num_variable[nml]; var++) {
        string = old_string variable[nml, var]", "
        if (length(string) > 80) {
          print "     ", old_string"&"
          old_string = variable[nml, var]", "
        } else {
          old_string = string
        }
      }
      print "     ", old_string variable[nml, var]
    }
    print ""
    print "    CALL input_from_file()"
    print ""
    print "    !"
    print "    ! set the default values"
    print "    !"
    for (nml = 1; nml <= num_namelist; nml++) {
      print "    ! namelist", namelist[nml]
      for (var = 1; var <= num_variable[nml]; var++) {
        if (default_[nml, var] != "") {
          print "   ", variable[nml, var], "=", default_[nml, var]
        }
      }
    }
    print ""
    print "    !"
    print "    ! read the namelist from stdin"
    print "    !"
    for (nml = 1; nml <= num_namelist; nml++) {
      print "    READ(stdin, NML="namelist[nml]", IOSTAT=ierr)"
      print "    CALL errore(__FILE__, \"error reading namelist", namelist[nml],"\", ierr)"
      print ""
    }
    print "    !"
    print "    ! copy the namelist variables to the output type"
    print "    !"
    for (nml = 1; nml <= num_namelist; nml++) {
      print "    ! namelist", namelist[nml]
      for (var = 1; var <= num_variable[nml]; var++) {
        print "   ", namelist[nml]"_t%"variable[nml, var], "=", variable[nml, var]
      }
      print ""
    }
    print "  END SUBROUTINE gw_input_read"
    print ""
    print "  !> broadcast input to all CPU"
    print_interface("gw_input_bcast")
    print ""
    print "    USE io_global, ONLY : meta_ionode_id"
    print "    USE mp,        ONLY : mp_bcast"
    print "    USE mp_world,  ONLY : world_comm"
    # for any private type, load the private bcast routine
    for (nml = 1; nml <= num_namelist; nml++) {
      for (var = 1; var <= num_variable[nml]; var++) {
        if (index(toupper(type[nml, var]), "TYPE") != 0) {
          print "    USE", variable[nml, var]"_module, ONLY : mp_bcast_"variable[nml, var]
        }
      }
    }
    print_type("INOUT")
    print "    !"
    print "    ! broadcast all variables"
    print "    !"
    for (nml = 1; nml <= num_namelist; nml++) {
      print "    ! namelist", namelist[nml]
      for (var = 1; var <= num_variable[nml]; var++) {
        # for default types use mp_bcast
        if (index(toupper(type[nml, var]), "TYPE") == 0) {
          print "    CALL mp_bcast("namelist[nml]"_t%"variable[nml, var]", meta_ionode_id, world_comm)"
        # for user defined type use seperate bcast routine
        } else {
          print "    CALL mp_bcast_"variable[nml, var]"("namelist[nml]"_t%"variable[nml, var]", meta_ionode_id, world_comm)"
        }
      }
      print ""
    }
    print "  END SUBROUTINE gw_input_bcast"
    print ""
    print "END MODULE gw_input_module"

  } else if (tiddler == 1) {

    if (path == "") {
      print "error: path not specified" > "/dev/stderr"
      exit 1
    }

    # main user guide with embedded tiddlers
    main = path"/Input variables.tid"
    print "title: Input variables"                           >  main
    print "type: text/vnd.tiddlywiki"                        >> main
    print ""                                                 >> main
    for (nml = 1; nml <= num_namelist; nml++) {
      print "!", namelist[nml]                               >> main
      print ""                                               >> main
      for (var = 1; var <= num_variable[nml]; var++) {
        print "{{"variable[nml, var]"||$:/template/input}}"  >> main
        print ""                                             >> main
      }
    }
    # extend the input description by parts not automatically documented
    print "{{ext_input}}"                                    >> main

    # individual tiddlers for all variables
    for (nml = 1; nml <= num_namelist; nml++) {
      for (var = 1; var <= num_variable[nml]; var++) {
         individual = path"/"variable[nml, var]".tid"
         print "title:", variable[nml, var]                  >  individual
         print "tags: Input", namelist[nml]                  >> individual
         print "type: text/vnd.tiddlywiki"                   >> individual
         print ""                                            >> individual
         print "|input_class|k"                              >> individual
         print "| ''type'' |", type[nml, var], "|"           >> individual
         if (default_[nml, var] == "''") {
           print "| ''default'' | ' ' |"                     >> individual
         } else if (default_[nml, var] != "") {
           print "| ''default'' |", default_[nml, var], "|"  >> individual
         }
         print description[nml, var]                         >> individual
      }
    }

  }
}
