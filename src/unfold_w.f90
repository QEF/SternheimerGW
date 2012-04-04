SUBROUTINE unfold_w(scrcoul_g_in)
USE kinds,         ONLY : DP
USE symm_base,     ONLY : nsym, s, time_reversal, t_rev, ftau, invs
USE gwsymm,        ONLY : ig_unique, ngmunique, use_symm
USE gvect,         ONLY : g, ngm, ecutwfc, nl
USE modes,         ONLY : nsymq, invsymq !, gi, gimq, irgq, irotmq, minus_q
USE gwsigma,       ONLY : ngmsco, sigma, sigma_g, nrsco, nlsco, fft6_g2r, ecutsco, ngmsig
USE freq_gw,       ONLY : fpol, fiu, nfs, nfsmax, nwcoul, wcoul
USE control_gw,    ONLY : zue, convt, rec_code, modielec, eta
USE qpoint,        ONLY : xq

IMPLICIT NONE

COMPLEX(DP)  :: scrcoul_g_in(ngmsig, ngmsig, nfs, 1)
INTEGER      :: ig, igp, npe, irr, icounter, ir, irp
!counter
INTEGER      :: isym, iwim
INTEGER      :: gmapsym(ngm,48)
COMPLEX(DP)  :: eigv(ngm,48)

!unpacks the symmetry reduced list of G vectors to fill the whole W 
!matrix before writing this to file, alternatively could just 
!write the symmetry reduced matrix to file much lower input/output...

CALL gmap_sym(nsym, s, ftau, gmapsym, eigv, invs)

if(.not.use_symm)GOTO 126
write(6,*) nsymq
IF(modielec) then
!only diagonal needs unfolding
      DO ig = 1, ngmunique
            DO iwim = 1, nfs
               DO isym = 1, nsymq
               !if gmap sym(ig,invsym) = ig_unique(ig) then that location will just get 
               !over written with the value that's already there... this shouldn't be a 
               !problem. 
                scrcoul_g_in(gmapsym(ig_unique(ig),invs(isym)), gmapsym(ig_unique(ig),invs(isym)),iwim, 1) = scrcoul_g_in(ig_unique(ig), ig_unique(ig), iwim,1) 
               ENDDO
            ENDDO
      ENDDO
ELSE
      DO ig = 1, ngmunique
         DO igp = 1, ngmunique
            DO iwim = 1, nfs
               DO isym = 1, nsymq
                !if gmap sym(ig,invsym) = ig_unique(ig) then that location will just get 
                !over written with the value that's already there... this shouldn't be a 
                !problem. 
                !scrcoul_g_in(gmapsym(ig,invs(isym)), gmapsym(igp,invs(isym)),iwim, 1) = scrcoul_g_in(ig_unique(ig), ig_unique(igp), iwim,1) 
                  scrcoul_g_in(gmapsym(ig_unique(ig),invs(isym)), gmapsym(ig_unique(igp), invs(isym)), iwim, 1) = scrcoul_g_in(ig_unique(ig), ig_unique(igp), iwim, 1)
               ENDDO
            ENDDO
         ENDDO
       ENDDO
ENDIF
126 CONTINUE
         write(6,'(3f14.7)')xq(:)
         do ig = 1, 50
             write(6,'(i4, f14.7)')ig, real( scrcoul_g_in(ig,ig,1,1))
         enddo
END SUBROUTINE unfold_w
