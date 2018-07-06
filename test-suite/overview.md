<!--
 
  This file is part of the SternheimerGW code.
  
  Copyright (C) 2010 - 2018 
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

Test suite for the SternheimerGW code
===========================

Features of the SternheimerGW code
------------------------

The following SternheimerGW features are tested in this suite

* system with a gap
* metallic system
* evaluate the dielectric function on the imaginary axis
* evaluate the dielectric function along the real axis
* select a specific k-point for the dielectric constant
* manual selection of solver for W
* use the direct solver to determine the dielectric constant
* use the iterative solver to determine the dielectric constant
* change the k-point grid compared to the scf calculation
* adjust the value used for the projection onto the valence band
* Godby-Needs plasmon pole model to obtain a denser grid on the imaginary axis
* AAA approximation to obtain a denser grid on the imaginary axis
* pole filter applied to AAA approximation
* full frequency integration using Pade approximation to interpolate
* symmetrize the frequency mesh to increase accuracy of Pade approximation
* obtain the Green's function on the imaginary axis
* use the multishift solver to obtain the Green's function
* convolute G and W on the imaginary axis to obtain the self energy
* convolute G and W on the real axis to obtain the self energy
* obtain the exchange self energy
* overcome the divergence of the Coulomb potential with spherical truncation
* overcome the divergence of the Coulomb potential with 2d truncation
* overcome the divergence of the Coulomb potential with Wigner-Seitz truncation
* evaluate the quasi particle eigenvalues using the Z factor
* obtain the frequency dependent spectral function for all bands
* restart from previously calculated W
* plotting the result of the Godby-Needs plasmon pole model
* plotting the result of the AAA approximation

Test case: gw\_bn
------------------

Evaluate the GW correction for a film of BN using the direct solver, plasmon
pole model, and imaginary frequency integration. We test both the 2d and the
Wigner-Seitz truncation.

The following SternheimerGW features are tested by this case

* system with a gap
* evaluate the dielectric function on the imaginary axis
* use the direct solver to determine the dielectric constant
* Godby-Needs plasmon pole model to obtain a denser grid on the imaginary axis
* obtain the Green's function on the imaginary axis
* use the multishift solver to obtain the Green's function
* obtain the exchange self energy
* overcome the divergence of the Coulomb potential with 2d truncation
* overcome the divergence of the Coulomb potential with Wigner-Seitz truncation
* evaluate the quasi particle eigenvalues using the Z factor
* obtain the frequency dependent spectral function for all bands
* restart from previously calculated W

Test case: gw\_c
-----------------

Evaluate the GW correction for C using the direct solver, imaginary frequency
integration, and full frequency integration for W.

The following SternheimerGW features are tested by this case

* system with a gap
* evaluate the dielectric function on the imaginary axis
* use the direct solver to determine the dielectric constant
* full frequency integration using Pade approximation to interpolate
* AAA approximation to obtain a denser grid on the imaginary axis
* pole filter applied to AAA approximation
* symmetrize the frequency mesh to increase accuracy of Pade approximation
* obtain the Green's function on the imaginary axis
* use the multishift solver to obtain the Green's function
* convolute G and W on the imaginary axis to obtain self energy
* convolute G and W on the real axis to obtain the self energy
* overcome the divergence of the Coulomb potential with spherical truncation
* evaluate the quasi particle eigenvalues using the Z factor
* obtain the frequency dependent spectral function for all bands
* restart from previously calculated W
* plotting the result of the AAA approximation

Test case: gw\_li
------------------

Evaluate the GW correction for Li using the direct solver, imaginary frequency
integration, and a plasmon-pole model for W.

The following SternheimerGW features are tested by this case

* metallic system
* evaluate the dielectric function on the imaginary axis
* use the direct solver to determine the dielectric constant
* change the k-point grid compared to the scf calculation
* adjust the value used for the projection onto the valence band
* Godby-Needs plasmon pole model to obtain a denser grid on the imaginary axis
* obtain the Green's function on the imaginary axis
* use the multishift solver to obtain the Green's function
* convolute G and W on the imaginary axis to obtain the self energy
* obtain the exchange self energy
* overcome the divergence of the Coulomb potential with spherical truncation
* evaluate the quasi particle eigenvalues using the Z factor
* obtain the frequency dependent spectral function for all bands

Test case: gw\_licl
--------------------

Evaluate the dielectric constant for LiCl using the direct solver. Manually
select the specialized solver for SternheimerGW.

The following SternheimerGW features are tested by this case

* system with a gap
* evaluate the dielectric function along the real axis
* use the direct solver to determine the dielectric constant
* select a specific k-point for the dielectric constant
* manual selection of solver for W

Test case: gw\_si
------------------

Evaluate the GW correction for Si using the direct solver, imaginary frequency
integration, and a plasmon-pole model for W. Then restart the calculation to
evaluate the self energy for the X point.
As a separate calculation, evaluate the GW correction with the iterative solver.

The following SternheimerGW features are tested by this case

* system with a gap
* evaluate the dielectric function on the imaginary axis
* use the direct solver to determine the dielectric constant
* use the iterative solver to determine the dielectric constant
* Godby-Needs plasmon pole model to obtain a denser grid on the imaginary axis
* obtain the Green's function on the imaginary axis
* use the multishift solver to obtain the Green's function
* convolute G and W on the imaginary axis to obtain the self energy
* obtain the exchange self energy
* overcome the divergence of the Coulomb potential with spherical truncation
* evaluate the quasi particle eigenvalues using the Z factor
* obtain the frequency dependent spectral function for all bands
* restart from previously calculated W
* plotting the result of the Godby-Needs plasmon pole model

