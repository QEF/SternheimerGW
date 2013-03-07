SUBROUTINE stern_symm()
!Finds a list of unique G vectors to run Sternheimer linear system on.
!I still can't find an easy way to prove that from the irreducible set
!of G vectors obtained here we uniquely reconstruct the full set of G vectors
!even though I'm pretty sure that's the case...

USE kinds,         ONLY : DP
USE symm_base,     ONLY : nsym, s, time_reversal, t_rev, ftau, invs
USE gwsigma,       ONLY : sigma, sigma_g, nrsco, nlsco, fft6_g2r, ecutsco, ngmpol
USE gwsymm,        ONLY : ngmunique, ig_unique, sym_ig, sym_friend
USE gvect,         ONLY : g, ngm, ecutwfc, nl
USE modes,         ONLY : nsymq, invsymq !, gi, gimq, irgq, irotmq, minus_q

IMPLICIT NONE
INTEGER      :: ig, igp, npe, irr, icounter, ir, irp
INTEGER      :: isym
INTEGER      :: gmapsym(ngm,48)
COMPLEX(DP)  :: eigv(ngm,48)
LOGICAL      :: unique_g


ig_unique(:) = 0
gmapsym(:,:) = 0

CALL gmap_sym(nsym, s, ftau, gmapsym, eigv, invs)

ngmunique = 0

!Find number of unique vectors:
write(6,'("Number of symmops in Small G_q: ", i4)'), nsymq
DO ig = 1, ngmpol
   unique_g = .true.
!Loop over symmetry operations in small group of q.
   DO isym = 1, nsymq
      DO igp = 1, ngmunique
         IF (gmapsym(ig,invs(isym)).eq.ig_unique(igp)) then
         !IF (gmapsym(ig,isym).eq.ig_unique(igp)) then
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
write(6,'("ngmpol and ngmunique", i4, i4)'), ngmpol, ngmunique

END SUBROUTINE stern_symm
