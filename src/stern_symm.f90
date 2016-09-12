!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! 
! Copyright (C) 2010 - 2016 
! Henry Lambert, Martin Schlipf, and Feliciano Giustino
!
! Sternheimer-GW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Sternheimer-GW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Sternheimer-GW. If not, see
! http://www.gnu.org/licenses/gpl.html .
!
!------------------------------------------------------------------------------ 
SUBROUTINE stern_symm()
!Finds a list of unique G vectors to run Sternheimer linear system on.
!I still can't find an easy way to prove that from the irreducible set
!of G vectors obtained here we uniquely reconstruct the full set of G vectors
!even though I'm pretty sure that's the case...

USE kinds,         ONLY : DP
USE symm_base,     ONLY : nsym, s, time_reversal, t_rev, ftau, invs
USE gwsigma,       ONLY : sigma_c_st, gcutcorr
USE gwsymm,        ONLY : ngmunique, ig_unique, sym_ig, sym_friend
USE gvect,         ONLY : g, ngm
USE control_gw,    ONLY : loqua
USE cell_base,     ONLY : at, bg
USE qpoint,        ONLY : xq
USE io_global,     ONLY : stdout
USE symm_base,     ONLY : s, t_rev, irt, ftau, nrot, nsym, &
                          time_reversal, copy_sym, inverse_s, s_axis_to_cart

IMPLICIT NONE

INTEGER      :: ig, igp, npe, irr, icounter, ir, irp
INTEGER      :: isym
INTEGER      :: gmapsym(ngm,48)
INTEGER      :: nsymq
COMPLEX(DP)  :: eigv(ngm,48)
LOGICAL      :: unique_g, invsymq
LOGICAL      :: minus_q, magnetic_sym, sym(48)

  ig_unique(:)  = 0
  gmapsym(:,:)  = 0
  sym_ig(:)     = 0
  sym_friend(:) = 0


  ngmunique = 0
  minus_q=.false.
  sym(1:nsym)=.true.
  call smallg_q (xq, 1, at, nsym, s, sym, minus_q)
  IF ( .not. time_reversal ) minus_q = .false.
 ! Here we re-order all rotations in such a way that true sym.ops.
 ! are the first nsymq; rotations that are not sym.ops. follow
  nsymq = copy_sym ( nsym, sym )
  call inverse_s ( )
 !check if inversion (I) is a symmetry. If so, there should be nsymq/2
 !symmetries without inversion, followed by nsymq/2 with inversion
 !Since identity is always s(:,:,1), inversion should be s(:,:,1+nsymq/2)
  invsymq = ALL ( s(:,:,nsymq/2+1) == -s(:,:,1) )
  if (invsymq)      WRITE(stdout,'(/5x, "qpoint HAS inversion symmetry")')
  if (.not.invsymq) WRITE(stdout,'(/5x, "qpoint does NOT have inversion symmetry")')
  WRITE(stdout,'(/5x, "nsym, nsymq, nrot ", i4, i4)') nsym,  nsymq
  CALL s_axis_to_cart () 
  do isym = 1, nsymq
   WRITE(6,'(3i4)') s(:,:,isym)
   WRITE(6,*)
   WRITE(6,'(3i4)') s(:,:,invs(isym))
   WRITE(6,*)
   WRITE(6,*)
  enddo
  CALL gmap_sym(nsym, s, ftau, gmapsym, eigv, invs)
!Find number of unique vectors:
ngmunique = 1
ig_unique(1) = 1
DO ig = 2, gcutcorr
   unique_g = .true.
!Loop over symmetry operations in small group of q.
   DO isym = 1, nsymq
      DO igp = 1, ngmunique
         IF (ig.eq.gmapsym(ig_unique(igp), invs(isym))) then
             unique_g = .false.
           !Rotation R^{-1}(ig) = ig_{un}
             sym_ig(ig) = isym
           !Corresponding unique vector
             sym_friend(ig) = ig_unique(igp)
         ENDIF
      ENDDO
   ENDDO
   IF(unique_g) then 
      !increment number of unique vectors by one and, keep track of its index.
      ngmunique = ngmunique + 1
      ig_unique(ngmunique) = ig
   ENDIF
ENDDO
write(6,'(/5x, "Number of symmops in Small G_q: ", i4)') nsymq
write(6,'(5x,  "ngmpol ", i4, " and ngmunique ", i4)') gcutcorr, ngmunique
END SUBROUTINE stern_symm
