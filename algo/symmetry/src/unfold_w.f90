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
SUBROUTINE unfold_w(num_g_corr, scrcoul_in, scrcoul_out)

USE cell_base,        ONLY : at
USE freq_gw,          ONLY : nfs
USE gwsymm,           ONLY : ig_unique, ngmunique, use_symm, sym_ig, sym_friend
USE io_global,        ONLY : stdout
USE kinds,            ONLY : DP
USE qpoint,           ONLY : xq
USE symm_base,        ONLY : nsym, s, time_reversal, ftau, invs, &
                             copy_sym, inverse_s, s_axis_to_cart

IMPLICIT NONE

!> the number of G vectors in the correlation grid
INTEGER,     INTENT(IN)  :: num_g_corr

!> the screened coulomb interaction only symmetrized elements
COMPLEX(DP), INTENT(IN)  :: scrcoul_in(num_g_corr, nfs, ngmunique)

!> the screened coulomb interaction all elements
COMPLEX(DP), INTENT(OUT) :: scrcoul_out(num_g_corr, num_g_corr, nfs)

COMPLEX(DP)  :: scrcoul_tmp(num_g_corr, nfs)
COMPLEX(DP)  :: phase
INTEGER      :: ig, igp
INTEGER      :: iwim
INTEGER      :: done, ngmdone
INTEGER      :: ngmdonelist(num_g_corr)
INTEGER      :: gmapsym(num_g_corr,nsym)
COMPLEX(DP)  :: eigv(num_g_corr,nsym)
REAL(DP)     :: xq_loc(3)
INTEGER      :: nsymq
LOGICAL      :: sym(48), minus_q, invsymq

!Unpacks the symmetry reduced list of G vectors to fill the whole W 
!matrix before writing this to file, alternatively could just 
!write the symmetry reduced matrix to file... but right now this isn't necessary.
!Again a local smallg_q so it doesn't conflict with solver_linter.
  minus_q=.false.
  sym(1:nsym)=.true.
  call smallg_q (xq, 1, at, nsym, s, sym, minus_q)
  IF ( .not. time_reversal ) minus_q = .false.
!Here we re-order all rotations in such a way that true sym.ops.
!are the first nsymq; rotations that are not sym.ops. follow
  nsymq = copy_sym ( nsym, sym )
  call inverse_s ( )
!check if inversion (I) is a symmetry. If so, there should be nsymq/2
!symmetries without inversion, followed by nsymq/2 with inversion
!Since identity is always s(:,:,1), inversion should be s(:,:,1+nsymq/2)
  invsymq = ALL ( s(:,:,nsymq/2+1) == -s(:,:,1) )
  if (invsymq)      WRITE(stdout,'(/5x, "qpoint HAS inversion symmetry")')
  if (.not.invsymq) WRITE(stdout,'(/5x, "qpoint does NOT have inversion symmetry")')
  WRITE(stdout,'(/5x, "nsym, nsymq ", i4, i4)') nsym,  nsymq
  ! Since the order of the s matrices is changed we need to recalculate:
  CALL s_axis_to_cart () 
  gmapsym(:,:) = 0
  CALL gmap_sym(num_g_corr, nsym, s, ftau, gmapsym, eigv, invs)
!
! reorder the input to the output array
!
  FORALL (ig = 1:ngmunique, igp = 1:num_g_corr, iwim=1:nfs)
    scrcoul_out(ig_unique(ig), igp, iwim) = CONJG(scrcoul_in(igp, iwim, ig))
  END FORALL
! trivial case
! no unfolding needs to be done
  IF (.NOT.use_symm .OR. nsymq == 1) THEN
    RETURN
  END IF ! no unfolding
!end trivial case
!
!stack ngmdone list with vectors that aren't unique:
  xq_loc = xq
  CALL cryst_to_cart(1, xq_loc(:), at, -1)
  ngmdonelist(:) = 0
  ngmdone = 0
  do ig = 1, ngmunique
     ngmdone = ngmdone + 1
     ngmdonelist(ngmdone) = ig_unique(ig)
  enddo
    DO ig = 1, num_g_corr
        DO done = 1, ngmdone
           if (ig.eq.ngmdonelist(done)) then
!               write(6,'("Cycling: unique or already unfolded.")') 
               GOTO 128
           endif
        ENDDO 
!still need to unfold this vector so we append it to the done list:
        ngmdone = ngmdone + 1
        ngmdonelist(ngmdone) = ig
!and unfold it with the correct correspondence ig has sym_friend(ig) where R^{-1} ig = ig_unique:
        DO iwim = 1, nfs
            DO igp = 1, num_g_corr
               scrcoul_tmp(igp, iwim) = scrcoul_out(sym_friend(ig), igp, iwim)
            ENDDO
        ENDDO
!the relationship R between ig and sym_friend(ig) is given by sym_ig.
          DO iwim = 1, nfs
             DO igp = 1, num_g_corr
!For symmetry operations with fraction translations we need to include:
!the \tau_{r} part which applies to the original G, G' rotation on R
!e^{-i2\pi(G - G')\cdot\tau_{R}} = eigv(G)*conjg(eigv(G'))
                phase = eigv(sym_friend(ig), sym_ig(ig))*conjg(eigv(igp, sym_ig(ig)))
                scrcoul_out(ig, gmapsym(igp, invs(sym_ig(ig))), iwim) = scrcoul_tmp(igp, iwim)*phase
             ENDDO
         ENDDO
128 CONTINUE
    ENDDO

END SUBROUTINE unfold_w
