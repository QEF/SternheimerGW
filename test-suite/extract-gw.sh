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

fname=$1

# determine whether we are in a PW or SternheimerGW step
scf=$(echo $fname | awk '/scf.in/{ print 1 }')
gw=$(echo $fname | awk '/gw.in/{ print 1 }')

if (( scf == 1 ))
then
  #
  # extract same quantities for SCF as QE test suite
  #
  e1=`grep ! $fname | tail -1 | awk '{printf "%12.6f\n", $5}'`
  n1=`grep 'convergence has' $fname | tail -1 | awk '{print $6}'`
  f1=`grep "Total force" $fname | head -1 | awk '{printf "%8.4f\n", $4}'`
  p1=`grep "P= " $fname | tail -1 | awk '{print $6}'`

  #
  # print the elements if they appear
  #
  if test "$e1" != ""; then
    echo e1
    echo "$e1"
  fi

  if test "$n1" != ""; then
    echo n1
    echo "$n1"
  fi

  if test "$f1" != ""; then
    echo f1
    echo "$f1"
  fi

  if test "$p1" != ""; then
    echo p1
    echo "$p1"
  fi

elif (( gw == 1 ))
then
  #
  # extract quantities from SternheimerGW
  #
  # real and imaginary part of the dielectric function
  re_eps=$(grep 'eps_{GG}(q,w)' $fname | awk '{ print $3 }')
  im_eps=$(grep 'eps_{GG}(q,w)' $fname | awk '{ print $4 }')
  if [[ $re_eps != "" ]]
  then
    echo "re_eps"
    echo "$re_eps"
  fi
  if [[ $im_eps != "" ]]
  then
    echo "im_eps"
    echo "$im_eps"
  fi

  # find the line number with the Re and Im part of Sigma and the spectral function
  line_re_sig=$(grep -n 'REsigma' $fname | awk -F":" '{ print $1 }')
  line_im_sig=$(grep -n 'IMsigma' $fname | awk -F":" '{ print $1 }')
  line_spec=$(grep -n 'ASpec' $fname | awk -F":" '{ print $1 }')

  # loop over elements
  i=0
  for line in $line_re_sig
  do
    # add new line character after first element
    if (( i > 0 ))
    then
      re_sig="$re_sig\n"
      im_sig="$im_sig\n"
      spec="$spec\n"
    fi

    (( i++ ))

    # determine number of frequencies
    num_line=$(echo $line_im_sig | awk "{ split(\$0, arr); print arr[$i] - $line - 2 }")

    # determine end point for data
    last_re_sig=$(echo "$line + $num_line" | bc -l)
    last_im_sig=$(echo $line_im_sig | awk "{ split(\$0, arr); print arr[$i] + $num_line }")
    last_spec=$(echo $line_spec | awk "{ split(\$0, arr); print arr[$i] + $num_line }")

    # extract real part of sigma
    data=$(head -$last_re_sig $fname | tail -$num_line)
    freq="$freq"$(echo "$data" | awk '{ print $1 }')
    re_sig="$re_sig"$(echo "$data" | awk '{ for ( i = 2; i <= NF; i++ ) { print $i } }')

    # extract imag part of sigma
    data=$(head -$last_im_sig $fname | tail -$num_line)
    freq="$freq"$(echo "$data" | awk '{ print $1 }')
    im_sig="$im_sig"$(echo "$data" | awk '{ for ( i = 2; i <= NF; i++ ) { print $i } }')

    # extract spectral function
    data=$(head -$last_spec $fname | tail -$num_line)
    freq="$freq"$(echo "$data" | awk '{ print $1 }')
    spec="$spec"$(echo "$data" | awk '{ for ( i = 2; i <= NF; i++ ) { print $i } }')
  done

  # print the elements if present
  if [[ $re_sig != "" ]]
  then
    echo "re_sig"
    echo -e "$re_sig"
  fi
  if [[ $im_sig != "" ]]
  then
    echo "im_sig"
    echo -e "$im_sig"
  fi
  if [[ $spec != "" ]]
  then
    echo "spec"
    echo -e "$spec"
  fi

  # maximum number of lines to read
  num_line=100

  # read DFT eigenvalues
  data=$(grep -A$num_line 'LDA eigenval (eV)' $fname | grep -B$num_line 'GWKpoint cart' | head -n -1)
  dft_eval=$(echo "$data" | awk '{ for ( i = 4; i <= NF; i++ ) { print $i } }')

  # read the QP eigenvalues
  data=$(grep -A$num_line 'GW qp energy (eV)' $fname | grep -B$num_line 'Vxc expt val (eV)' | head -n -1)
  gw_eval=$(echo "$data" | awk '{ for ( i = 5; i <= NF; i++ ) { print $i } }')

  # read the Vxc expectation value
  data=$(grep -A$num_line 'Vxc expt val (eV)' $fname | grep -B$num_line 'Sigma_ex val (eV)' | head -n -1)
  vxc_exp=$(echo "$data" | awk '{ for ( i = 5; i <= NF; i++ ) { print $i } }')

  # read the HF self energy
  data=$(grep -A$num_line 'Sigma_ex val (eV)' $fname | grep -B$num_line 'QP renorm' | head -n -1)
  sigma_x=$(echo "$data" | awk '{ for ( i = 4; i <= NF; i++ ) { print $i } }')

  # read the QP eigenvalues
  data=$(grep -A$num_line 'QP renorm' $fname | grep -B$num_line 'REsigma' | head -n -1)
  z_factor=$(echo "$data" | awk '{ for ( i = 3; i <= NF; i++ ) { print $i } }')

  # print the elements if present
  if [[ $dft_eval != "" ]]
  then
    echo "dft_eval"
    echo "$dft_eval"
  fi
  if [[ $gw_eval != "" ]]
  then
    echo "gw_eval"
    echo "$gw_eval"
  fi
  if [[ $vxc_exp != "" ]]
  then
    echo "vxc_exp"
    echo "$vxc_exp"
  fi
  if [[ $sigma_x != "" ]]
  then
    echo "sigma_x"
    echo "$sigma_x"
  fi
  if [[ $z_factor != "" ]]
  then
    echo "z_factor"
    echo "$z_factor"
  fi

  # print the screened coulomb of the plotting function

  line_begin_plot=$(grep -n 'Plotting.*approximation' $fname | awk -F":" '{ print $1 }')
  line_end_plot=$(grep -n 'End of.*approximation' $fname | awk -F":" '{ print $1 }')

  # loop over elements
  i=0
  for line in $line_begin_plot
  do
    # add new line character after first element
    if (( i > 0 ))
    then
      re_coul="$re_coul\n"
      im_coul="$im_coul\n"
    fi

    (( i++ ))

    # determine number of frequencies
    num_line=$(echo $line_end_plot | awk "{ split(\$0, arr); print arr[$i] - $line - 2 }")

    # determine end point of data
    last_line=$(echo "$line + $num_line" | bc -l)

    # extract plotted Coulomb potential
    data=$(head -$last_line $fname | tail -$num_line)
    re_coul="$re_coul"$(echo "$data" | awk 'NF == 4 { print $3 }')
    im_coul="$im_coul"$(echo "$data" | awk 'NF == 4 { print $4 }')

  done

  # print the elements if present
  if [[ $re_coul != "" ]]
  then
    echo "re_coul"
    echo -e "$re_coul"
  fi
  if [[ $re_coul != "" ]]
  then
    echo "im_coul"
    echo -e "$im_coul"
  fi

else
  echo "Unknown input file" > /dev/stderr

fi
