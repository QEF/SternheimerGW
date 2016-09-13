Release notes
=============

Version 0.10
------------

Add a new timing module that tracks a bit better the work that is actually done
by the code so that it is more easily visible which parts of the code use how
much time.

Provide a new linear_operator module that acts as a wrapper for the h_psi call.
This module can be used from both the Coulomb and the Green linear solver.

Implement a module that wraps the MPI allgatherv and gatherv routines. This
allows to use those routines instead of the MPI sum routine that is currently
used so that we do not communicate more than necessary.

Implement a BiCGstab(l) multishift solver and use it to evaluate the Greens
function and the screened Coulomb interaction replacing the BiCG solver used
before. It is more stable and converges in less iterations. However it is a
bit slower for the Greens function because the multishift solver is not
optimized. As the preconditioning and restarting from a previous solution is
not implemented at the moment the iterative solver is a bit slower, too.

Rewrite the 6 dimensional Fourier transform in a new module that follows a
bit closer the interface of the conventional 3 dimensional Fourier transforms.

Rewrite the frequency integration of the self-energy. The real frequency
integration should work now, but is a bit unstable, so the imaginary frequency
integration is still the preferred choice.

The Wigner-Seitz truncation from QE is linked. Their implementation is not
very efficient, but if turns out to give good results, we may implement our
own Wigner-Seitz truncation.

Change Sigma_x calculation from real to reciprocal space which is more in line
with the implementation of the hybrid functionals and a bit faster.

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
