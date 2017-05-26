project: SternheimerGW
project_dir: ./src
project_website: http://www.sternheimergw.org/
output_dir: ./Doc/ford
summary: SternheimerGW is a code for calculating the dielectric screening and quasiparticle properties in a range of materials.
author: Martin Schlipf, Henry Lambert, and Feliciano Giustino
author_description: The Materials Modelling and Design Group.
website: http://giustino.materials.ox.ac.uk/index.php
email: martin.schlipf@materials.ox.ac.uk
predocmark: >
media_dir: ./media
docmark_alt: #
predocmark_alt: <
display: public
         private
source: false
graph: true
macro: TEST
       LOGIC=.true.
license: by-nc-sa
extra_filetypes: sh #

The SternheimerGW code evaluates the quasi-particle correction to the DFT
results in the GW approximation. The Greens function G and the screened
Coulomb W interaction are evaluated solving Sternheimer equations (linear
response problem) instead of summing over the virtual, unoccupied subspace.

The main target of SternheimerGW are very accurate calculations for crystals
and surfaces. We can evaluate the full band structure along any user selected
k-point path. In addition to the quasi-particle corrections to the
eigenvalues, we provide the full frequency dependent spectral function,
which can be compared to angular-resolved photo-emission spectroscopy (ARPES)
experiments.
