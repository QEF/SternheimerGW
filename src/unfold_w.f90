  !-----------------------------------------------------------------------
  ! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !-----------------------------------------------------------------------
SUBROUTINE unfold_w(scrcoul_g_in, iq)
USE kinds,         ONLY : DP
USE symm_base,     ONLY : nsym, s, time_reversal, t_rev, ftau, invs
USE gwsymm,        ONLY : ig_unique, ngmunique, use_symm, sym_ig, sym_friend
USE gvect,         ONLY : g, ngm, nl
USE freq_gw,       ONLY : fpol, fiu, nfs, nfsmax, nwcoul, wcoul
USE control_gw,    ONLY : zue, convt, rec_code, modielec, eta, godbyneeds, padecont
USE qpoint,        ONLY : xq
USE cell_base,     ONLY : at
USE gwsigma,       ONLY : sigma_c_st
USE noncollin_module, ONLY : noncolin, nspin_mag
USE io_global,        ONLY : stdout
USE symm_base,        ONLY : s, t_rev, irt, ftau, nrot, nsym, &
                             time_reversal, copy_sym, inverse_s, s_axis_to_cart
USE cell_base,        ONLY : at, bg

IMPLICIT NONE

COMPLEX(DP)  :: scrcoul_g_in(sigma_c_st%ngmt, sigma_c_st%ngmt, nfs, nspin_mag)
COMPLEX(DP)  :: scrcoul_g_tmp(sigma_c_st%ngmt, nfs)
COMPLEX(DP)  :: phase
INTEGER      :: ig, igp, npe, irr, icounter, ir, irp
INTEGER      :: isym, iwim, iq, iw
INTEGER      :: done, ngmdone, isp
INTEGER      :: ngmdonelist(sigma_c_st%ngmt)
INTEGER      :: gmapsym(ngm,48)
COMPLEX(DP)  :: eigv(ngm,48)
LOGICAL      :: not_unique 
REAL(DP)     :: xq_loc(3)
INTEGER      :: nsymq
LOGICAL      :: sym(48), minus_q, invsymq

!Unpacks the symmetry reduced list of G vectors to fill the whole W 
!matrix before writing this to file, alternatively could just 
!write the symmetry reduced matrix to file... but right now this isn't necessary.
!Again a local smallg_q so it doesn't conflict with solver_linter.
  minus_q=.false.
  sym(1:nsym)=.true.
  call smallg_q (xq, 1, at, bg, nsym, s, ftau, sym, minus_q)
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
  WRITE(stdout,'(/5x, "nsym, nsymq, nrot ", i4, i4)') nsym,  nsymq
  ! Since the order of the s matrices is changed we need to recalculate:
  CALL s_axis_to_cart () 
  gmapsym(:,:) = 0
  CALL gmap_sym(nsym, s, ftau, gmapsym, eigv, invs)
!Cases where no unfolding needs to be done:
  if(.not.use_symm)GOTO 126
  if(nsymq.eq.1)GOTO 126
!end Cases
!stack ngmdone list with vectors that aren't unique:
  xq_loc = xq
  CALL cryst_to_cart(1, xq_loc(:), at, -1)
  ngmdonelist(:) = 0
  ngmdone = 0
  do ig = 1, ngmunique
     ngmdone = ngmdone + 1
     ngmdonelist(ngmdone) = ig_unique(ig)
  enddo
IF(modielec) then
!only diagonal needs unfolding:
      DO ig = 1, ngmunique
         DO done = 1, ngmdone
            if (ig.eq.ngmdonelist(done)) then
!               write(6,'("Cycling: unique or already unfolded.")')
                CYCLE
            endif
         ENDDO
         DO iwim = 1, nfs
            DO isym = 1, nsymq
               scrcoul_g_in(gmapsym(ig_unique(ig),invs(isym)), gmapsym(ig_unique(ig),invs(isym)),iwim,1) = &
                 scrcoul_g_in(ig_unique(ig), ig_unique(ig), iwim, 1)
            ENDDO
         ENDDO
      ENDDO
ELSE
    DO ig = 1, sigma_c_st%ngmt
        DO done = 1, ngmdone
           if (ig.eq.ngmdonelist(done)) then
!               write(6,'("Cycling: unique or already unfolded.")') 
               GOTO 128
           endif
        ENDDO 
!still need to unfold this vector so we append it to the done list:
        ngmdone = ngmdone + 1
        ngmdonelist(ngmdone) = ig
        DO isp = 1, nspin_mag
!and unfold it with the correct correspondence ig has sym_friend(ig) where R^{-1} ig = ig_unique:
        DO iwim = 1, nfs
            DO igp = 1, sigma_c_st%ngmt
               scrcoul_g_tmp(igp,iwim) = scrcoul_g_in(sym_friend(ig), igp, iwim, isp)
            ENDDO
        ENDDO
!the relationship R between ig and sym_friend(ig) is given by sym_ig.
          DO iwim = 1, nfs
             DO igp = 1, sigma_c_st%ngmt
!For symmetry operations with fraction translations we need to include:
!the \tau_{r} part which applies to the original G, G' rotation on R
!e^{-i2\pi(G - G')\cdot\tau_{R}} = eigv(G)*conjg(eigv(G'))
                phase = eigv(sym_friend(ig), sym_ig(ig))*conjg(eigv(igp, sym_ig(ig)))
                scrcoul_g_in(ig, gmapsym(igp, invs(sym_ig(ig))), iwim, isp) = scrcoul_g_tmp(igp, iwim)*phase
             ENDDO
         ENDDO
        ENDDO
128 CONTINUE
    ENDDO
ENDIF

126 CONTINUE
END SUBROUTINE unfold_w
