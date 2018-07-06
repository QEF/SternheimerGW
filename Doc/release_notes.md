Release notes
=============

Version 0.15
------------

Update FFT type to recent change in Quantum Espresso.
Adjust to renaming of commits in Quantum Espresso.

Further improvement to real frequency integration. Symmetrize enforcing
that +w and -w are the same. Fix bug where frequencies are in the wrong
quadrant.

Employ external analytic continuation library.

> compatible QE commit 6.3 7357cdb

Version 0.14
------------

New analytic continuation scheme: AAA approximation which expands a given
set of points onto a barycentric representation. Add test case for this new
feature.

> compatible QE version 6.2 revision 13949

Version 0.13
------------

New input reader with more consistent variable names and automatic
generation from the documentation of the user variables in a YAML
file. Cleanup of unused variables.

Implement the possiblity to plot W on the real axis using the Pade
or plasmon-pole approximation.

New examples for new input format: silicon, lithium, BN film, and
LiCl (dielectric constant).

Adopt new XML scheme and FFT layout from PW. Update test-suite.

> compatible QE version 6.2-beta revision 13814

Version 0.12
------------

Implement metallic systems. For metallic systems the projector on the valence
band should be converged to 0. In the future a frequency dependent projector
might be implemented.

From now on the code is referred to as SternheimerGW to differentiate it more
clearly from QSGW and other GW implementations. To align to this new name, a
few routines had to be renamed that still used the old handle.

Print timings during the run of the code.
> compatible QE revision 13531 (note: compile with old XML)

Version 0.11
------------

Reinterpret the real frequency integration. We choose the frequencies for the
Green's function above or below the real axis depending on whether we are above
or below the Fermi energy. This leads directly includes the nonanalytic part,
fixes the inconsistencies between real and imaginary frequency integration, and
improves the convergence of the real frequency integration with respect to the
number of necessary frequencies. To improve the stability of the frequency
integration, we use the symmetry of the screened Coulomb interaction to extend
the number of points used for the Pade continuation.

Add new linear solver that construct a Krylov subspace that is used for all
frequencies but might be extended if one requires a tighter convergence setting.
Dynamically change the linear solver if the first choice of the solver does
not converge. Change the linear solver to return error codes and switch to
a different solver if nonzero error code is returned.

Add G parallelization to the second index in the convolution of G and W. This
is more expensive, because the FFT now requires communication, but the only
solution to calculate large memory systems. In the new approach the memory
requirement for the Greens function is distributed across the processes leading
to a reduced overall memory consumption.

Add the debug module that allows to set debugging options. If the code is
compiled with the \_\_DEBUG flag it will examine the requested parts in more
detail. The only current option is to check (H - w) G = -delta. The unit
test can be used to investigate the behavior of the linear problem.

Create the SternheimerGW testsuite that tests all major features implemented
in SternheimerGW. The test cases are not converged, but any change made to the
code should reproduce the results obtained with the test. If an error larger
than the threshold pops up, the commit should be reviewed with care.
> compatible QE version 6.1 revision 13374

Version 0.10.1
--------------

Bugfix to make systems without inversion symmetry work. There where some arrays
that were the complex conjugate of the correct value, which did not matter with
inversion symmetry but resulted in incorrect results for systems without it.
> compatible QE version 6.0 revision 13079


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
> compatible QE version 6.0 revision 13079

Version 0.9
-----------

Create a new module to write/read Sigma in xml format. The advantage of this
module is that it contains some metadata that makes processing the self-energy
much easier. However, the iotk module is significantly slower than the regular
file I/O so that for performance reasons, we currently keep the old reading
routines in SternheimerGW where the metadata is not needed.
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
