SUBROUTINE write_sigma(sigma, iw0)
  USE io_global,     ONLY : stdout, ionode_id, ionode
  USE kinds,         ONLY : DP
  USE klist,         ONLY : wk, xk
  USE wvfct,         ONLY : nbnd, npw, npwx, igk, g2kin, et
  USE eqv,           ONLY : evq, eprec
  USE freq_gw,       ONLY : fpol, fiu, nfs, nfsmax,&
                            nwcoul, nwgreen, nwalloc, nwsigma, wtmp, wcoul,&
                            wgreen, wsigma, wsigmamin, wsigmamax,&
                            deltaw, wcoulmax, ind_w0mw, ind_w0pw
  USE units_gw,      ONLY : iuncoul, iungreen, iunsigma, lrsigma, lrcoul, lrgrn, iuwfc, lrwfc
  USE qpoint,        ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE gwsigma,       ONLY : ngmsco, nrsco, nlsco, fft6_g2r, ecutsco, ngmsig,&
                            nr1sco, nr2sco, nr3sco, ngmgrn, ngmpol
  USE gvect,         ONLY : g, ngm, ecutwfc, nl
  USE cell_base,     ONLY : tpiba2, tpiba, omega, alat, at
  USE symm_base,     ONLY : nsym, s, time_reversal, t_rev, ftau, invs, nrot
  USE mp_global,     ONLY : inter_pool_comm, intra_pool_comm, mp_global_end, mpime, npool, &
                            nproc_pool, me_pool, my_pool_id, nproc
  USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
  USE fft_scalar,    ONLY : cfft3d
  USE modes,         ONLY : nsymq, invsymq, gi, gimq, irgq, irotmq, minus_q
  USE wavefunctions_module, ONLY : evc
  USE control_flags,        ONLY : noinv

IMPLICIT NONE

COMPLEX (DP) :: sigma_g (ngmsco, ngmsco)
COMPLEX(DP)  :: sigma   (nrsco, nrsco) 
INTEGER      :: ig, igp, irr, icounter, ir, irp, iw0

   sigma_g = (0.0d0,0.0d0)
   call fft6(sigma_g(1,1), sigma(1,1), -1)

   do igp = ngmsco + 1, nrsco
      do ig = ngmsco + 1, nrsco
              sigma (ig, igp) = (0.0d0, 0.0d0)
      enddo 
   enddo

   do igp = 1, ngmsco
     do ig = 1, ngmsco
        sigma_g(ig,igp)  = sigma(ig,igp)
     enddo
   enddo
   CALL davcio (sigma_g, lrsigma, iunsigma, iw0, 1)
END SUBROUTINE write_sigma
