!
! Copyright (C) 2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
  !----------------------------------------------------------------------------
SUBROUTINE setup_green (xq)
  !----------------------------------------------------------------------------
  ! ... This routine initializes variables for the non-scf calculations required to generate the 
  ! ... wave functions and solve the green linear system it requires 
  ! ... all the points in the small group of G_{"k-q"}.
  ! ... from the list of symmetry operations of this group we can then generate \sum_{R} W(R^{-1}q).
  ! ... required for the convolution to form sigma_k.
  ! ... Needed on input (read from data file):
  ! ... "nsym" crystal symmetries s, ftau, t_rev, "nrot" lattice symetries "s"
  ! ... "nkstot" k-points in the irreducible BZ wrt lattice symmetry
  ! ... Produced on output:
  ! ... symmetries ordered with the "nsym of k-q" GW symmetries first
  ! ... "nkstot" k- and k+q-points in the IBZ calculated for the GW symmetries.)
  ! ... Misc. data needed for running the non-scf calculation

  USE kinds,              ONLY : DP
  USE constants,          ONLY : eps8
  USE parameters,         ONLY : npk
  USE io_global,          ONLY : stdout
  USE constants,          ONLY : pi, degspin
  USE cell_base,          ONLY : at, bg, alat, tpiba, tpiba2, ibrav, omega
  USE ions_base,          ONLY : nat, tau, ntyp => nsp, ityp, zv
  USE force_mod,          ONLY : force
  USE basis,              ONLY : natomwfc
  USE gvect,              ONLY : nr1, nr2, nr3
  USE klist,              ONLY : xk, wk, nks, nelec, degauss, lgauss, &
                                 nkstot, qnorm
  USE lsda_mod,           ONLY : lsda, nspin, current_spin, isk, &
                                 starting_magnetization
  USE symm_base,          ONLY : s, t_rev, irt, ftau, nrot, nsym, &
                                 time_reversal, sname, d1, d2, d3, &
                                 copy_sym, s_axis_to_cart
  USE wvfct,              ONLY : nbnd, nbndx
  USE control_flags,      ONLY : ethr, isolve, david, &
                                 noinv, modenum, use_para_diag
  USE mp_global,          ONLY : kunit
  USE spin_orb,           ONLY : domag
  USE noncollin_module,   ONLY : noncolin
  USE start_k,            ONLY : nks_start, xk_start, wk_start
  USE paw_variables,      ONLY : okpaw
  USE modes,              ONLY : nsymq, invsymq !, gi, gimq, irgq, irotmq, minus_q
  !
  IMPLICIT NONE
  !
  REAL (DP), INTENT(IN) :: xq(3)
  !
  REAL (DP), ALLOCATABLE :: rtau (:,:,:)
  LOGICAL  :: minus_q, magnetic_sym, sym(48)
  !
  INTEGER, EXTERNAL :: n_atom_wfc
  !
  IF ( .NOT. ALLOCATED( force ) ) ALLOCATE( force( 3, nat ) )
  !
  ! ... threshold for diagonalization ethr - should be good for all cases
  !
  ethr= 1.0D-9 / nelec
  !
  ! ... variables for iterative diagonalization (Davidson is assumed)
  !

  isolve = 0
  david = 4
  nbndx = david*nbnd
  natomwfc = n_atom_wfc( nat, ityp )

  !
#ifdef __PARA
  IF ( use_para_diag )  CALL check_para_diag( nelec )
#else
  use_para_diag = .FALSE.
#endif

  ! ... Symmetry and k-point section
  ! ... time_reversal = use q=>-q symmetry for k-point generation

  magnetic_sym = noncolin .AND. domag 
  time_reversal = .NOT. noinv .AND. .NOT. magnetic_sym

  !
  ! ... smallg_q flags to false the  symmetry operations of the crystal
  ! ... that are not symmetry operations of the small group of q
  !

  sym(1:nsym)=.true.

! xq should = k find all symmetry operations that leave k invariant.
  call smallg_q (xq, modenum, at, bg, nsym, s, ftau, sym, minus_q)

  IF ( .not. time_reversal ) minus_q = .false.

  ! ... for single-mode calculation: find symmetry operations
  ! ... that leave the chosen mode unchanged. Note that array irt 
  ! ... must be available: it is allocated and read from xml file 
  ! Here we re-order all rotations in such a way that true sym.ops.
  ! are the first nsymq; rotations that are not sym.ops. follow
  !

  nsymq = copy_sym ( nsym, sym )

  !
  ! check if inversion (I) is a symmetry. If so, there should be nsymq/2
  ! symmetries without inversion, followed by nsymq/2 with inversion
  ! Since identity is always s(:,:,1), inversion should be s(:,:,1+nsymq/2)
  !

  invsymq = ALL ( s(:,:,nsymq/2+1) == -s(:,:,1) )

  !  Since the order of the s matrices is changed we need to recalculate:

  call s_axis_to_cart () 

  !
  ! ... Input k-points are assumed to be  given in the IBZ of the Bravais
  ! ... lattice, with the full point symmetry of the lattice.
  !

  nkstot = nks_start
  xk(:,1:nkstot) = xk_start(:,1:nkstot)

  !write(6,'(3f11.7)') xk(:,:)

  wk(1:nkstot)   = wk_start(1:nkstot)

  !
  ! ... If some symmetries of the lattice are missing in the crystal,
  ! ... "irreducible_BZ" computes the missing k-points.
  !

  CALL irreducible_BZ (nrot, s, nsymq, minus_q, at, bg, npk, nkstot, xk, wk, &
                       t_rev)
  
  ! NOW xk, wk should have the appropriate kpoints and weights for all the points required in the green's 
  ! function. 

  IF ( lsda ) THEN
     !
     ! ... LSDA case: two different spin polarizations,
     ! ...            each with its own kpoints
     !
     if (nspin /= 2) call errore ('setup','nspin should be 2; check iosys',1)
     !
     CALL set_kup_and_kdw( xk, wk, isk, nkstot, npk )
     !
  ELSE IF ( noncolin ) THEN
     !
     ! ... noncolinear magnetism: potential and charge have dimension 4 (1+3)
     !
     if (nspin /= 4) call errore ('setup','nspin should be 4; check iosys',1)
     current_spin = 1
     !
  ELSE
     !
     ! ... LDA case: the two spin polarizations are identical
     !
     wk(1:nkstot)    = wk(1:nkstot) * degspin
     current_spin = 1
     !
     IF ( nspin /= 1 ) &
        CALL errore( 'setup', 'nspin should be 1; check iosys', 1 )
     !
  END IF
  !
  IF ( nkstot > npk ) CALL errore( 'setup', 'too many k points', nkstot )
  !
  ! ...notice: qnorm is used by allocate_nlpot to determine
  ! the correct size of the interpolation table "qrad"
  !
  qnorm = sqrt(xq(1)**2 + xq(2)**2 + xq(3)**2)
  !
#ifdef __PARA
  !
  ! ... set the granularity for k-point distribution
  !
  IF ( ABS( xq(1) ) < eps8 .AND. ABS( xq(2) ) < eps8 .AND. &
       ABS( xq(3) ) < eps8 ) THEN
     !
     kunit = 1
     !
  ELSE
     !
     kunit = 2
     !
  ENDIF
  !
  ! ... distribute k-points (and their weights and spin indices)
  !
  !HL this should be cause for concern...
  CALL divide_et_impera( xk, wk, isk, lsda, nkstot, nks )
  !
#else
  !
  nks = nkstot
  !
#endif
  !
  RETURN
  !
END SUBROUTINE setup_nscf
