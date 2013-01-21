SUBROUTINE unfold_w(scrcoul_g_in, iq)
USE kinds,         ONLY : DP
USE symm_base,     ONLY : nsym, s, time_reversal, t_rev, ftau, invs
USE gwsymm,        ONLY : ig_unique, ngmunique, use_symm, sym_ig, sym_friend
USE gvect,         ONLY : g, ngm, ecutwfc, nl
USE modes,         ONLY : nsymq, invsymq !, gi, gimq, irgq, irotmq, minus_q
USE gwsigma,       ONLY : ngmsco, sigma, sigma_g, nrsco, nlsco, fft6_g2r, ecutsco, ngmpol
USE freq_gw,       ONLY : fpol, fiu, nfs, nfsmax, nwcoul, wcoul
USE control_gw,    ONLY : zue, convt, rec_code, modielec, eta, godbyneeds, padecont
USE qpoint,        ONLY : xq

IMPLICIT NONE

COMPLEX(DP)  :: scrcoul_g_in(ngmpol, ngmpol, nfs, 1)
COMPLEX(DP)  :: scrcoul_g_tmp(ngmpol, nfs)
COMPLEX(DP)  :: phase

INTEGER      :: ig, igp, npe, irr, icounter, ir, irp
!counter
INTEGER      :: isym, iwim, iq
INTEGER      :: done, ngmdone
!INTEGER      :: ngmdonelist(nsymq+1)
INTEGER      :: ngmdonelist(ngmpol)

INTEGER      :: gmapsym(ngm,48)
COMPLEX(DP)  :: eigv(ngm,48)

LOGICAL   :: not_unique

!unpacks the symmetry reduced list of G vectors to fill the whole W 
!matrix before writing this to file, alternatively could just 
!write the symmetry reduced matrix to file... but right now this isn't necessary.

gmapsym(:,:) = 0
CALL gmap_sym(nsym, s, ftau, gmapsym, eigv, invs)

do isym = 1, nsymq
   WRITE(6,'(3i4)') s(:,:,isym) 
   WRITE(6,*)
   WRITE(6,'(3i4)') s(:,:,invs(isym)) 
   WRITE(6,*)
   WRITE(6,*)
enddo

!Cases where no unfolding needs to be done:
if(.not.use_symm)GOTO 126
if(nsymq.eq.1)GOTO 126
!end Cases

ngmdonelist(:) = 0
ngmdone = 0

!stack ngmdone list with vectors that aren't unique:

do ig = 1, ngmunique
   ngmdone = ngmdone + 1
   ngmdonelist(ngmdone) = ig_unique(ig)
enddo

IF(modielec) then
!only diagonal needs unfolding:
      DO ig = 1, ngmunique
         DO iwim = 1, nfs
            DO isym = 1, nsymq
               scrcoul_g_in(gmapsym(ig_unique(ig),invs(isym)), gmapsym(ig_unique(ig),invs(isym)),iwim,1) = scrcoul_g_in(ig_unique(ig), ig_unique(ig), iwim,1)
            ENDDO
         ENDDO
      ENDDO
ELSE
    DO ig = 1, ngmpol
        DO done = 1, ngmdone
           if (ig.eq.ngmdonelist(done)) then
               write(6,'("Cycling: unique or already unfolded.")') 
               GOTO 128
           endif
        ENDDO 
!still need to unfold this vector so we append it to the done list:
            ngmdone = ngmdone + 1
            ngmdonelist(ngmdone) = ig

!and unfold it with the correct correspondence ig has sym_friend(ig) where R^{-1} ig = ig_unique:
        DO iwim = 1, nfs
            DO igp = 1, ngmpol
               scrcoul_g_tmp(igp,iwim) = scrcoul_g_in(sym_friend(ig), igp, iwim, 1)
            ENDDO
        ENDDO

!the relationship R between ig and sym_friend(ig) is given by sym_ig.
        DO iwim = 1, nfs
            DO igp = 1, ngmpol
            !For symmetry operations with fraction translations we need to include:
            !e^{-i2\pi(G - G')\cdot\tau_{r}} = eigv(G)*conjg(eigv(G'))
              phase = eigv(sym_friend(ig), sym_ig(ig))*conjg(eigv(igp, sym_ig(ig)))
              scrcoul_g_in(ig, gmapsym(igp, invs(sym_ig(ig))), iwim, 1) = scrcoul_g_tmp(igp, iwim)*phase
            ENDDO
        ENDDO

128 CONTINUE
    ENDDO
!    write(6,*) ngmdonelist(:)
ENDIF

!Populate tmp array with the unique row of perturbations 
!tmp array required or we fold on top of each other ...
!            DO iwim = 1, nfs
!               DO igp = 1, ngmpol
!                   scrcoul_g_tmp(igp,iwim) = scrcoul_g_in(ig_unique(ig), igp, iwim, 1)
                  !scrcoul_g_tmp(gmapsym(igp,invs(isym)), iwim) = scrcoul_g_in(ig_unique(ig), igp, iwim, 1)
!               ENDDO
!            ENDDO
!Rotate the unique vector to the W_{q} (R^{-1}G, WR^{-1}G').
!I should have a completed shell of unique G vectors so there shouldn't be problems with mapping outside of ngmpol
!however I'm sure a number of pathological cases exist for different bravais lattices, etc...
 !           DO iwim = 1, nfs
 !              DO igp = 1, ngmpol
 !                scrcoul_g_in(gmapsym(ig, invs(sym_ig(ig, ))), gmapsym(igp, invs(isym), iwim, 1) = scrcoul_g_tmp(igp, iwim)
!               ENDDO
!            ENDDO

126 CONTINUE

!Diagonal
!Zero wings of W:
IF(iq.eq.1) then
  Write(6, '("Zeroing Wings of W.")')
  if(godbyneeds) then
     do igp = 2, ngmpol
        scrcoul_g_in(1,igp,1,1) = (10.0d0, 0.0d0)
        scrcoul_g_in(1,igp,2,1) = ( 0.0d0, 0.0d0)
     enddo
     do igp = 2, ngmpol
        scrcoul_g_in(igp,1,1,1) = (10.0d0, 0.0d0)
        scrcoul_g_in(igp,1,2,1) = ( 0.0d0, 0.0d0)
     enddo
  endif
!How to zero for pade continuation?
  if(padecont) then
     do igp = 2, ngmpol
        do iwim = 1, nfs
           scrcoul_g_in(1,igp,iwim,1) = (0.0d0, 0.0d0)
           scrcoul_g_in(igp,1,iwim,1) = (0.0d0, 0.0d0)
        enddo
     enddo
  endif
ENDIF
         do ig = 20, 50
             write(6,'(i4, f14.7)')ig, real(scrcoul_g_in(ig,ig,1,1))
         enddo

         do ig = 1, 14
             write(6,'(14f14.7)')real(scrcoul_g_in(ig,1:14,1,1))
         enddo

         do ig = 1, 14
             write(6,'(14f14.7)')aimag(scrcoul_g_in(ig,1:14,1,1))
         enddo
END SUBROUTINE unfold_w
