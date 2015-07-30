!
! Copyright (C) 2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------

SUBROUTINE setup_nscf_green(xq)

  !----------------------------------------------------------------------------
  !
  ! ... This routine initializes variables for the non-scf calculations at k 
  ! ... and k+q required by the linear response calculation at finite q.

  ! I think I can just use the monkhorst pack generated q-mesh to do a run.
  ! This approach means doing an nscf step for each of the unique q-points to generate
  ! the full k, k+q grid with all the weights and nice symmetry reduction properties built in.
  ! ... Here we find the symmetry group of the crystal that leaves
  ! ... the GW q-vector (xq) unchanged; determines the k- and k+q points in the irreducible BZ
  ! ... Needed on input (read from data file):
  ! ... "nsym" crystal symmetries s, ftau, t_rev, "nrot" lattice symetries "s"
  ! ... "nkstot" k-points in the irreducible BZ wrt lattice symmetry
  ! ... Produced on output:
  ! ... symmetries ordered with the "nsymq" GW symmetries first
  ! ... "nkstot" k- and k+q-points in the IBZ calculated for the GW symmetries.)
  ! ... Misc. data needed for running the non-scf calculation

  USE kinds,              ONLY : DP
  USE parameters,         ONLY : npk
  USE io_global,          ONLY : stdout
  USE constants,          ONLY : pi, degspin, eps8
  USE cell_base,          ONLY : at, bg
  USE ions_base,          ONLY : nat, tau, ityp, zv
  USE basis,              ONLY : natomwfc
  USE klist,              ONLY : xk, wk, nks, nelec, degauss, lgauss, &
                                 nkstot, qnorm
  USE lsda_mod,           ONLY : lsda, nspin, current_spin, isk, &
                                 starting_magnetization
  USE wvfct,              ONLY : nbnd, nbndx
  USE control_flags,      ONLY : ethr, isolve, david, &
                                 noinv, modenum, use_para_diag, max_cg_iter
  USE mp_pools,           ONLY : kunit
  USE spin_orb,           ONLY : domag
  USE noncollin_module,   ONLY : noncolin
  USE start_k,            ONLY : nks_start, xk_start, wk_start, nk1, nk2, nk3,&
                                 k1, k2, k3
  USE paw_variables,      ONLY : okpaw
  USE modes,              ONLY : nsymq, invsymq
  USE disp,               ONLY : xk_kpoints
  USE uspp_param,         ONLY : n_atom_wfc
  USE symm_base,          ONLY : s, t_rev, irt, ftau, nrot, nsym, &
                                 time_reversal, copy_sym, inverse_s, s_axis_to_cart
  USE control_gw,         ONLY : newgrid, lgamma

  !
  IMPLICIT NONE
  !
  REAL (DP), INTENT(IN) :: xq(3)
  !
  LOGICAL  :: minus_q, magnetic_sym, sym(48)
  !
  INTEGER   :: ik
  !
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
  max_cg_iter=20
  natomwfc = n_atom_wfc( nat, ityp, noncolin )
  !
#ifdef __PARA
  IF ( use_para_diag )  CALL check_para_diag( nbnd )
#else
  use_para_diag = .FALSE.
#endif

  ! ... Symmetry and k-point section
  ! ... time_reversal = use q=>-q symmetry for k-point generation

  magnetic_sym = noncolin .AND. domag 
  time_reversal = .NOT. noinv .AND. .NOT. magnetic_sym
  sym(1:nsym)=.true.
  minus_q=.false.

  call smallg_q (xq, modenum, at, bg, nsym, s, ftau, sym, minus_q)
  IF ( .not. time_reversal ) minus_q = .false.
  ! Here we re-order all rotations in such a way that true sym.ops.
  ! are the first nsymq; rotations that are not sym.ops. follow
   nsymq = copy_sym ( nsym, sym )

   call inverse_s ( )

  ! check if inversion (I) is a symmetry. If so, there should be nsymq/2
  ! symmetries without inversion, followed by nsymq/2 with inversion
  ! Since identity is always s(:,:,1), inversion should be s(:,:,1+nsymq/2)
    invsymq = ALL ( s(:,:,nsymq/2+1) == -s(:,:,1) )
    if (invsymq)      WRITE(6,'(\5x, "qpoint HAS inversion symmetry")')
    if (.not.invsymq) WRITE(6,'(\5x, "qpoint does NOT have inversion symmetry")')

    WRITE(6,'(\5x, "nsym, nsymq, nrot ", i4, i4, i4)') nsym,  nsymq, nrot 
    WRITE(6,'(\5x, "nsym, nsymq, nrot ", i4, i4)') nks_start, nsymq, nkstot
    WRITE(6,*) newgrid 

  ! Since the order of the s matrices is changed we need to recalculate:
    call s_axis_to_cart () 
  ! ... Input k-points are assumed to be  given in the IBZ of the Bravais
  ! ... lattice, with the full point symmetry of the lattice.
  if( nks_start > 0 .AND. .NOT. newgrid ) then
     !
     !  In this case I keep the same points of the Charge density
     !  calculations
     !
     nkstot = nks_start
     xk(:,1:nkstot) = xk_start(:,1:nkstot)
     wk(1:nkstot)   = wk_start(1:nkstot)
  else
     !In this case I generate a new set of k-points
     CALL kpoint_grid ( nsym, time_reversal, .false., s, t_rev, &
                        bg, nk1*nk2*nk3, k1,k2,k3, nk1,nk2,nk3, nkstot, xk, wk)
     !CALL kpoint_grid ( nrot, time_reversal, .false., s, t_rev, &
     !                   bg, nk1*nk2*nk3, k1,k2,k3, nk1,nk2,nk3, nkstot, xk, wk)
  endif

  ! ... If some symmetries of the lattice no longer apply for this kpoint
  ! ... "irreducible_BZ" generates the missing k-points with the reduced number of
  ! ... symmetry operations. 

  CALL irreducible_BZ (nsym, s, nsymq, minus_q, magnetic_sym, &
                       at, bg, npk, nkstot, xk, wk, t_rev)

  !wk(contains weights 
  CALL set_kplusq( xk, wk, xq, nkstot, npk )
  !
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
  IF (lgamma) THEN
  !
       kunit = 1
  !
  ELSE
  !
     kunit = 2
  !
  ENDIF
  !

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
END SUBROUTINE setup_nscf_green

!
!-----------------------------------------------------------------------
subroutine smallg_q (xq, modenum, at, bg, nrot, s, ftau, sym, minus_q)
  !-----------------------------------------------------------------------
  !
  ! This routine selects, among the symmetry matrices of the point group
  ! of a crystal, the symmetry operations which leave q unchanged.
  ! Furthermore it checks if one of the above matrices send q --> -q+G.
  ! In this case minus_q is set true.
  !
  !  input-output variables
  !
  USE kinds, ONLY : DP
  implicit none

  real(DP), parameter :: accep = 1.e-5_dp

  real(DP), intent(in) :: bg (3, 3), at (3, 3), xq (3)
  ! input: the reciprocal lattice vectors
  ! input: the direct lattice vectors
  ! input: the q point of the crystal

  integer, intent(in) :: s (3, 3, 48), nrot, ftau (3, 48), modenum
  ! input: the symmetry matrices
  ! input: number of symmetry operations
  ! input: fft grid dimension (units for ftau)
  ! input: fractionary translation of each symmetr
  ! input: main switch of the program, used for
  !        q<>0 to restrict the small group of q
  !        to operation such that Sq=q (exactly,
  !        without G vectors) when iswitch = -3.
  logical, intent(inout) :: sym (48), minus_q
  ! input-output: .true. if symm. op. S q = q + G
  ! output: .true. if there is an op. sym.: S q = - q + G
  !
  !  local variables
  !

  real(DP) :: aq (3), raq (3), zero (3)
  ! q vector in crystal basis
  ! the rotated of the q vector
  ! the zero vector

  integer :: irot, ipol, jpol
  ! counter on symmetry op.
  ! counter on polarizations
  ! counter on polarizations

  logical :: eqvect
  ! logical function, check if two vectors are equa
  !
  ! return immediately (with minus_q=.true.) if xq=(0,0,0)
  !
  minus_q = .true.
  if ( (xq (1) == 0.d0) .and. (xq (2) == 0.d0) .and. (xq (3) == 0.d0) ) &
       return
  !
  !   Set to zero some variables
  !
  minus_q = .false.
  zero(:) = 0.d0
  !
  !   Transform xq to the crystal basis
  !
  aq = xq
  call cryst_to_cart (1, aq, at, - 1)
  !
  !   Test all symmetries to see if this operation send Sq in q+G or in -q+G
  !
  do irot = 1, nrot
     if (.not.sym (irot) ) goto 100
     raq(:) = 0.d0
     do ipol = 1, 3
        do jpol = 1, 3
           raq(ipol) = raq(ipol) + DBLE( s(ipol,jpol,irot) ) * aq( jpol)
        enddo
     enddo
     sym (irot) = eqvect (raq, aq, zero, accep)
     !
     !  if "iswitch.le.-3" (modenum.ne.0) S must be such that Sq=q exactly !
     !
     if (modenum.ne.0 .and. sym(irot) ) then
        do ipol = 1, 3
           sym(irot) = sym(irot) .and. (abs(raq(ipol)-aq(ipol)) < 1.0d-5)
        enddo
     endif
!     if (.not.minus_q) then
     if (sym(irot).and..not.minus_q) then
        raq = - raq
        minus_q = eqvect (raq, aq, zero, accep)
     endif
100  continue
  enddo
  !
  !  if "iswitch.le.-3" (modenum.ne.0) time reversal symmetry is not included !
  !
  if (modenum.ne.0) minus_q = .false.
  !
  return
end subroutine smallg_q

