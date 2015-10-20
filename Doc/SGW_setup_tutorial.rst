=====
Sternheimer-GW
=====

The Sternheimer-GW code can be used to obtain the dielectric and quasiparticle
properties of materials in a linear response formalism. The methodology
has been described here [Giustino10]_ and here [Lambert13]_.

Getting Started
======

Checkout the latest version of espresso-5.x.x/
  svn checkout http://qeforge.qe-forge.org/svn/q-e/trunk/espresso ./espresso-5.1/
  username: anonymous
  password: 

Make all:
  cd ./espresso-5.x.x/
  ./config
  make pw

Checkout the latest version of SGW:
  svn co https://rothery3.materials.ox.ac.uk/svn/sgw/branches/sgw-5.0/ ./SGW
  cd ./SGW/src
  make

