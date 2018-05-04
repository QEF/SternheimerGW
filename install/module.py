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
depend = "depend"
structure = {
"util": {
  depend: ["base", "lrmods"],
},
"data": {
  "algebra": {
    depend: ["base"],
  },
  "fft": {
    depend: ["base", "timing"],
  },
  "parallel": {
    depend: ["base"],
  },
  "timing": {
    depend: ["base"],
  },
},
"algo": {
  "analytic": {
    depend: ["base", "util", "data", "grid"],
  },
  "grid": {
    depend: ["base", "util", "data"],
  },
  "io": {
    depend: ["base", "pw", "util", "grid"],
  },
  "linear_solver": {
    depend: ["base", "pw", "lrmods", "util", "data"]
  },
  "mesh": {
    depend: ["base", "pw", "util"],
  },
  "nscf": {
    depend: ["base", "pw", "lrmods", "util"],
  },
  "reorder": {
    depend: ["base"],
  },
  "setup": {
    depend: ["base", "pw", "lrmods", "util"],
  },
  "symmetry": {
    depend: ["base", "pw", "lrmods", "util"],
  },
  "teardown": {
    depend: ["base", "pw", "lrmods", "util"],
  },
  "truncation": {
    depend: ["base", "util"],
  },
},
"phys": {
  "corr": {
    depend: ["base", "pw", "util", "data", "algo", "coul", "green"],
  },
  "coul": {
    depend: ["base", "pw", "lrmods", "util", "data", "algo", "postproc"],
  },
  "exch": {
    depend: ["base", "pw", "lrmods", "util", "data", "algo"],
  },
  "green": {
    depend: ["base", "pw", "util", "data", "algo"],
  },
  "matrix_el": {
    depend: ["base", "pw", "lrmods", "util", "data", "algo"],
  },  
  "postproc": {
    depend: ["base", "algo"],
  },
},
}
