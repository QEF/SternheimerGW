<!--
 
  This file is part of the SternheimerGW code.
  
  Copyright (C) 2010 - 2017 
  Henry Lambert, Martin Schlipf, and Feliciano Giustino
 
  SternheimerGW is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
 
  SternheimerGW is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with SternheimerGW. If not, see
  http://www.gnu.org/licenses/gpl.html .
 
-->

SternheimerGW examples
======================

In this directory, we show a few exemplary calculations that demonstrate
features of the SternheimerGW software. In general, these examples are tuned
to compromise between runtime and accuracy. If you want to perform similar
calculations for a publication, you would need to increase the numerical 
cutoff parameters.

All these examples have a common structure, they contain the input files
scf.in, which is required to run the DFT calculation, and gw.in, which is
used to run the SternheimerGW calculation. In addition, they contain the
required pseudopotentials and reference output files (scf.ref and gw.ref).
When you are in the source directory of a particular example, type
~~~~
../../../bin/pw.x < scf.in > scf.out
../../bin/gw.x < gw.in > gw.out
~~~~
to run the example. Note that this assumes the default installation, where
SternheimerGW is a subdirectory of the main Quantum ESPRESSO directory.

Example 1: Silicon
------------------

The aim of this calculation is to calculate the electronic band structure of
silicon along the Gamma-X direction. In this example, we use full-frequency
integration along the imaginary axis and a 4 x 4 x 4 grid for both k point
and q points. The energy cutoffs for exchange and correlation are 15 and 6 Ry,
respectively.

Example 2: Li
-------------

The aim of this calculation is to calculate the spectral function of lithium
at the Gamma point. In this example, we use full-frequency integration along
the imaginary axis and a 6 x 6 x 6 grid for both k points and q points. The
energy cutoffs for exchange and correlation are 20 and 8 Ry, respectively.

Example 3: BN film
-----------------

The aim of this calculation is to calculate the quasi-particle corrections for
a BN film at the Gamma point. In this example, we use a Godby-Needs plasmon-pole
model along the imaginary axis and a 8 x 8 x 1 grid for k and q points. The
energy cutoffs for exchange and correlation are 20 and 5 Ry, respectively.

Example 4: LiCl
---------------

The aim of this calculation is to obtain the dielectric function on the real
frequency axis with a 6 x 6 x 6 k point mesh at the Gamma q point.
