Release notes
=============

Version 0.9
-----------

Create a new module to write/read Sigma in xml format. The advantage of this
module is that it contains some metadata that makes processing the self-energy
much easier. However, the iotk module is significantly slower than the regular
file I/O so that for performance reasons, we currently keep the old reading
routines in SGW where the metadata is not needed.
This version also contains a few changes to accomodate recent changes in QE.
The igk is removed from the wave function module, so that we create our own
igk arrays where needed. In the future, this should be replaced by igk_k.
The API of the FFT routines changed so that we had to adjust a few routines
relying on them.
> compatible QE revision 12824


Version 0.8
-----------

Create a new module to evaluate the expectation value of Sigma.
The idea is to collect parts that are similar for Sigma_x and Sigma_c into
a single file, to avoid the duplication of code.
It also makes it easier to evaluate the expectation value in other parts of
the code should it be necessary.
> compatible QE revision 12630


Version 0.7
-----------

Added new routine to reorder arrays according to an igk array.
Added automatic documentation options (Doxygen/Ford).
Added some checks if the solver converged.
Cleanup some unused files and variables within files.
> compatible QE revision 12581


Version 0.6
-----------

This is the initial version from which we started documenting the changes.
From now on, we plan to document all changes in the form of release notes.
The plan is commit new versions frequently, so that the changes between
versions are small and easy to document.
> compatible QE revision 12534
