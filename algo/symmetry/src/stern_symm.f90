!------------------------------------------------------------------------------
!
! This file is part of the SternheimerGW code.
! 
! Copyright (C) 2010 - 2018
! Henry Lambert, Martin Schlipf, and Feliciano Giustino
!
! SternheimerGW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! SternheimerGW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with SternheimerGW. If not, see
! http://www.gnu.org/licenses/gpl.html .
!
!------------------------------------------------------------------------------ 
SUBROUTINE stern_symm(num_g_corr)
!Finds a list of unique G vectors to run Sternheimer linear system on.
!I still can't find an easy way to prove that from the irreducible set
!of G vectors obtained here we uniquely reconstruct the full set of G vectors
!even though I'm pretty sure that's the case...

USE cell_base,     ONLY : at
USE gwsymm,        ONLY : ngmunique, ig_unique, sym_ig, sym_friend
USE io_global,     ONLY : stdout
USE kinds,         ONLY : DP
USE qpoint,        ONLY : xq
USE symm_base,     ONLY : nsym, s, time_reversal, ftau, invs, &
                          copy_sym, inverse_s, s_axis_to_cart

IMPLICIT NONE

  !> the number of G vectors in the correlation grid
  INTEGER, INTENT(IN) :: num_g_corr

INTEGER      :: ig, igp
INTEGER      :: isym
INTEGER      :: gmapsym(num_g_corr,nsym)
INTEGER      :: nsymq
COMPLEX(DP)  :: eigv(num_g_corr,nsym)
LOGICAL      :: unique_g, invsymq
LOGICAL      :: minus_q, sym(48)

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
  WRITE(stdout,'(/5x, "nsym, nsymq ", i4, i4)') nsym,  nsymq
  CALL s_axis_to_cart () 
  do isym = 1, nsymq
   WRITE(6,'(3i4)') s(:,:,isym)
   WRITE(6,*)
   WRITE(6,'(3i4)') s(:,:,invs(isym))
   WRITE(6,*)
   WRITE(6,*)
  enddo
  CALL gmap_sym(num_g_corr, nsym, s, ftau, gmapsym, eigv, invs)
!Find number of unique vectors:
ngmunique = 1
ig_unique(1) = 1
DO ig = 2, num_g_corr
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
write(6,'(5x,  "ngmpol ", i4, " and ngmunique ", i4)') num_g_corr, ngmunique
END SUBROUTINE stern_symm
