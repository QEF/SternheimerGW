Release notes
=============

Version 0.8
-----------

Create a new module to evaluate the expectation value of Sigma.
The idea is to collect parts that are similar for Sigma_x and Sigma_c into
a single file, to avoid the duplication of code.
It also makes it easier to evaluate the expectation value in other parts of
the code should it be necessary.


Version 0.7
-----------

Added new routine to reorder arrays according to an igk array.
Added automatic documentation options (Doxygen/Ford).
Added some checks if the solver converged.
Cleanup some unused files and variables within files.


Version 0.6
-----------

This is the initial version from which we started documenting the changes.
From now on, we plan to document all changes in the form of release notes.
The plan is commit new versions frequently, so that the changes between
versions are small and easy to document.
