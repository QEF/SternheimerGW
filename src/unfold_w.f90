SUBROUTINE unfold_w(scrcoul_g_in, iq)
USE kinds,         ONLY : DP
USE symm_base,     ONLY : nsym, s, time_reversal, t_rev, ftau, invs
USE gwsymm,        ONLY : ig_unique, ngmunique, use_symm, sym_ig, sym_friend
USE gvect,         ONLY : g, ngm, nl
USE modes,         ONLY : nsymq, invsymq 
USE freq_gw,       ONLY : fpol, fiu, nfs, nfsmax, nwcoul, wcoul
USE control_gw,    ONLY : zue, convt, rec_code, modielec, eta, godbyneeds, padecont
USE qpoint,        ONLY : xq
USE cell_base,     ONLY : at
USE gwsigma,      ONLY : sigma_c_st

IMPLICIT NONE

COMPLEX(DP)  :: scrcoul_g_in(sigma_c_st%ngmt, sigma_c_st%ngmt, nfs, 1)
COMPLEX(DP)  :: scrcoul_g_tmp(sigma_c_st%ngmt, nfs)
COMPLEX(DP)  :: phase

INTEGER      :: ig, igp, npe, irr, icounter, ir, irp
!counter
INTEGER      :: isym, iwim, iq, iw
INTEGER      :: done, ngmdone
INTEGER      :: ngmdonelist(sigma_c_st%ngmt)
INTEGER      :: gmapsym(ngm,48)
COMPLEX(DP)  :: eigv(ngm,48)
LOGICAL      :: not_unique
REAL(DP)     :: xq_loc(3)

!unpacks the symmetry reduced list of G vectors to fill the whole W 
!matrix before writing this to file, alternatively could just 
!write the symmetry reduced matrix to file... but right now this isn't necessary.

gmapsym(:,:) = 0
CALL gmap_sym(nsym, s, ftau, gmapsym, eigv, invs)
!do isym = 1, nsymq
!   WRITE(6,'(3i4)') s(:,:,isym) 
!   WRITE(6,*)
!   WRITE(6,'(3i4)') s(:,:,invs(isym)) 
!   WRITE(6,*)
!   WRITE(6,*)
!enddo
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
               scrcoul_g_in(gmapsym(ig_unique(ig),invs(isym)), gmapsym(ig_unique(ig),invs(isym)),iwim,1) = scrcoul_g_in(ig_unique(ig), ig_unique(ig), iwim,1)
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
!and unfold it with the correct correspondence ig has sym_friend(ig) where R^{-1} ig = ig_unique:
        DO iwim = 1, nfs
            DO igp = 1, sigma_c_st%ngmt
               scrcoul_g_tmp(igp,iwim) = scrcoul_g_in(sym_friend(ig), igp, iwim, 1)
            ENDDO
        ENDDO
!the relationship R between ig and sym_friend(ig) is given by sym_ig.
        DO iwim = 1, nfs
            DO igp = 1, sigma_c_st%ngmt
            !For symmetry operations with fraction translations we need to include:
            !the \tau_{r} part which applies to the original G, G' rotation on R
            !e^{-i2\pi(G - G')\cdot\tau_{R}} = eigv(G)*conjg(eigv(G'))
             phase = eigv(sym_friend(ig), sym_ig(ig))*conjg(eigv(igp, sym_ig(ig)))
             scrcoul_g_in(ig, gmapsym(igp, invs(sym_ig(ig))), iwim, 1) = scrcoul_g_tmp(igp, iwim)*phase
            ENDDO
        ENDDO
128 CONTINUE
    ENDDO
ENDIF

126 CONTINUE
END SUBROUTINE unfold_w
