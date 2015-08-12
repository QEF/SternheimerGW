  !-----------------------------------------------------------------------
  ! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !-----------------------------------------------------------------------
!
!----------------------------------------------
SUBROUTINE grid_shuffle()
!----------------------------------------------
!
!  In PH they have a very simple set of meshes:
!     nksq = nks / 2
!     ALLOCATE(ikks(nksq), ikqs(nksq))
!     DO ik=1,nksq
!        ikks(ik) = 2 * ik - 1
!        ikqs(ik) = 2 * ik
!     ENDDO
!And they just generate the eigenfunctions on those grids with an nscf step for each new qpoint.
!I think there will come a time when we need to do the full grid shuffle. One NSCF step to generate all the k's.
!This should be sufficient to do a full GW dispersion calculation. Since all we need are all the ks, k+qs
!and k-qs. The current refold routine RELIES on a uniform and complete q-mesh. The trick is now to reformulate
!the folding so that it applies when we become more selective about which k and q points we want to use and we
!start exploiting symmetry reductions. 

!generate uniform {k} and {k \pm q} grids.
!The {k} grid is taken to coincide with the {q} grid generated
!in gwhs.f90, the {k+q} grid is obtained by folding
!FG:In this way we calculate the occupied states only once in gwhs.f90
 
  USE kinds,           ONLY : DP
  USE units_gw,        ONLY : iuwfc, lrwfc
  USE qpoint,          ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE disp,            ONLY : nqs, g0vec, eval_occ, gmap, x_q
  USE klist,           ONLY : xk, wk, nkstot, nks
  USE constants,       ONLY : eps8
  USE wvfct,           ONLY : et, nbnd
  USE gvect,           ONLY : ngm

  IMPLICIT NONE

  INTEGER :: iq, count, i, j, k, ik, ipol, ikk, ikq, ig0, igmg0, nw
  INTEGER :: ig, iw, igp, ierr, ibnd, ios, recl, unf_recl, ikf, fold(nkstot), icounter
  REAL    :: g0(3)
  INTEGER, PARAMETER :: ng0vec = 27
  INTEGER :: shift(nkstot)
  !complex(kind=DP) :: aux (ngm, nbnd_occ) , evq (ngm, nbnd_occ)
  COMPLEX(DP) :: aux (ngm, nbnd) , evq (ngm, nbnd)

! SGW integer, parameter :: nksq = nq, nks = 2 * nksq
! In SGW the xk grid is defined to be the q grid here
! However we already have the k grid from the nscf step
! do ik = 1, nksq
!   xk (:, ik) = xq (:, ik)
!   include spin degeneracy
!   wk = two / float ( nksq )
! enddo
! find folding index and G-vector map
! "nkstot" k- and k+q-points in the IBZ calculated for the phonon sym.
!SGW Pilot  integer, parameter :: nksq = nq, nks = 2 * nksq
! HL gmap appears to be fine here...
!  do ig=1,ngm
!     do ig0=1, ng0vec
!         write(6,*), gmap (ig,ig0)
!     enddo
!  enddo
!  do ik = 1, nksq
!  do ik = 1, nqs

   do ik = 1, nks
     call ktokpmq (xk(:,ik), xq, +1, fold(ik) )
    !
    !  the folding G-vector
    !
     g0 = xk(:, fold(ik)) - ( xk(:,ik) + xq )
    !
     write(6,'("folding index ", 3f7.3)') fold(ik)
     write(6,'("xk original ", 3f7.3)') xk(:, ik)
     write(6,'("xq ", 3f7.3)') xq
     write(6,'("xk folded ", 3f7.3)') xk(:, fold(ik))
     write(6,'("g0 ", 3f7.3)') g0
    !
     shift(ik) = 0
     if(ik.eq.1) then
       write(6,'("g0vec and ig0")')
        do ig0 = 1, ng0vec
        if ( ( abs(g0vec(1,ig0)-g0(1)).lt.eps8 ) .and. &
             ( abs(g0vec(2,ig0)-g0(2)).lt.eps8 ) .and. &
             ( abs(g0vec(3,ig0)-g0(3)).lt.eps8 ) ) then
             shift(ik) = ig0
             write(6,'(3f7.3)') g0vec(:,ig0)
             write(6,'(3f7.3)') g0
         endif
        enddo
     endif
!   if (shift(ik).eq.0) call error ('coulomb','folding vector not found',0)
   enddo
!stop
! double the grid and add k+q in the even positions
! PH.x set up the list like this in initialize_gw:
!     DO ik=1,nksq
!        ikks(ik) = 2 * ik - 1
!        ikqs(ik) = 2 * ik
!     ENDDO
 
   do ik = nksq, 1, -1
     ikk = 2 * ik - 1
     ikq = 2 * ik
     xk(:,ikk) = xk(:,ik)
     xk(:,ikq) = xk(:,ik) + xq
     wk(ikk) = wk(ik)
     wk(ikq) = 0.d0
   enddo

  do ik = 1, nksq
     !
     ikk = 2 * ik - 1
     ikq = 2 * ik
     !
     !  the folded k+q
     !
     ikf = fold(ik)
     write(6,*) fold(ik)
     !
     !  read occupied wavefunctions for the folded k+q point
     !  c_{k+q} (G) = c_{k+q+G0} (G-G0)

     !read ( iuwfc, rec = 2 * ikf - 1, iostat = ios) aux

     call davcio(aux, lrwfc, iuwfc, 2 * ikf - 1, -1)

     !  set the phase factor of evq for the folding
     !  WARNING: here we loose some information since
     !  the sphere G-G0 has half the surface outside the
     !  cutoff. It is important to make sure that the cutoff
     !  used is enough to guarantee that the lost half-surface
     !  does not create problems. In the EPW code I solved
     !  the issue by translating the G-sphere in the fft
     !  mapping. It's painful, so I will do it only in extremis.
     !  FG Aug 2008

     write(6,*) shift(ik)
     do ibnd = 1, nbnd
        do ig = 1, ngm

         ! problem with shift(ik)=0 not finding G-vector which maps k+q back on to k. 
         ! the mesh may not be self-contained. 

          igmg0 = gmap (ig, shift(ik))
          if (igmg0.eq.0) then
             evq (ig,ibnd) = (0.0d0, 0.0d0) 
          else
             write(6,*) igmg0
             evq (ig,ibnd) = aux ( igmg0, ibnd)
          endif
        enddo
     enddo
     !stop
     !  and write it again in the right place
     !  HL write ( iuwfc, rec = ikq, iostat = ios) evq

      call davcio (evq, lrwfc, iuwfc, ikq, 1)

     !SGW:
     !FG: DEBUG: here I checked that by running
     ! call eigenstates ( xk(:,ikq), vr, g2kin, evq, eval_occ(:,ikf) )
     ! the evq (wfs at k+q) are the same as those obained above
     ! (within a phase factor and gauge in degenerate cases -
     ! I actually checked sum_ibnd |evq(:,ibnd)|^2 )
     ! the eigenvalues

     !HL: need to switch around these indices. In sgw eval_occ contains the eigenvalues on the initially
     ! generated q grid. In gw.x et contains the NSCF eigenvalues.
     ! SGW:
     ! et(:,ikk) = eval_occ(:,ik)
     ! et(:,ikq) = eval_occ(:,ikf)
     ! Maybe it would be best to define a new temporary array which gets passed to solve 
     ! linter and contains the shifted eigenvalues. since the grid shuffle needs to take place for each q.  
     ! i.e. call gridshuffle(iq) -> etq(:,:) -> call solve_linter (evc,eqv, etq()).

     eval_occ(:, ikk) = et(:,ik)
     eval_occ(:, ikq) = et(:,ikf)
    
  enddo !on k
END SUBROUTINE grid_shuffle

