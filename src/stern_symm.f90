  !-----------------------------------------------------------------------
  ! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !-----------------------------------------------------------------------
SUBROUTINE stern_symm()
!Finds a list of unique G vectors to run Sternheimer linear system on.
!I still can't find an easy way to prove that from the irreducible set
!of G vectors obtained here we uniquely reconstruct the full set of G vectors
!even though I'm pretty sure that's the case...

USE kinds,         ONLY : DP
USE symm_base,     ONLY : nsym, s, time_reversal, t_rev, ftau, invs
USE gwsigma,       ONLY : sigma_c_st
USE gwsymm,        ONLY : ngmunique, ig_unique, sym_ig, sym_friend
USE gvect,         ONLY : g, ngm
USE modes,         ONLY : nsymq, invsymq 
USE control_gw,    ONLY : loqua

IMPLICIT NONE

INTEGER      :: ig, igp, npe, irr, icounter, ir, irp
INTEGER      :: isym
INTEGER      :: gmapsym(ngm,48)
COMPLEX(DP)  :: eigv(ngm,48)
LOGICAL      :: unique_g


ig_unique(:)  = 0
gmapsym(:,:)  = 0
sym_ig(:)     = 0
sym_friend(:) = 0

CALL gmap_sym(nsym, s, ftau, gmapsym, eigv, invs)

ngmunique = 0

if(loqua) then
  do isym = 1, nsymq
   WRITE(6,'(3i4)') s(:,:,isym)
   WRITE(6,*)
   WRITE(6,'(3i4)') s(:,:,invs(isym))
   WRITE(6,*)
   WRITE(6,*)
  enddo
endif


!Find number of unique vectors:

ngmunique = 1
ig_unique(1) = 1

DO ig = 2, sigma_c_st%ngmt
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
write(6,'(/5x, "Number of symmops in Small G_q: ", i4)'), nsymq
write(6,'(5x,  "ngmpol ", i4, " and ngmunique ", i4)'), sigma_c_st%ngmt, ngmunique
END SUBROUTINE stern_symm
