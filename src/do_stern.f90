SUBROUTINE do_stern()
  USE io_global,  ONLY : stdout
  USE kinds,      ONLY : DP
  USE disp,       ONLY : nqs, num_k_pts, xk_kpoints, w_of_q_start
  USE gwsigma,    ONLY : sigma_c_st
  USE gwsymm,     ONLY : ngmunique, ig_unique, use_symm, sym_friend, sym_ig
  USE control_gw, ONLY : done_bands, reduce_io, recover, tmp_dir_gw,&
                          ext_restart, bands_computed, bands_computed, nbnd_occ, lgamma,&
                          do_q0_only, solve_direct, tinvert, lrpa, do_epsil
  USE freq_gw,    ONLY : nfs


IMPLICIT NONE

  INTEGER :: iq, ik, ig, igstart, igstop
  COMPLEX(DP), ALLOCATABLE :: scrcoul_g(:,:,:,:)
  LOGICAL :: do_band, do_iq, setup_pw, exst, do_matel

  ALLOCATE ( scrcoul_g( sigma_c_st%ngmt, sigma_c_st%ngmt, nfs, 1))
  ALLOCATE ( ig_unique( sigma_c_st%ngmt) )
  ALLOCATE ( sym_ig(sigma_c_st%ngmt))
  ALLOCATE ( sym_friend(sigma_c_st%ngmt))

  do_iq=.TRUE.
  setup_pw = .TRUE.
  do_band  = .TRUE.
  do_matel = .TRUE.

  DO iq = w_of_q_start, nqs
    scrcoul_g(:,:,:,:) = (0.0d0, 0.0d0)
    CALL prepare_q(do_band, do_iq, setup_pw, iq)
    do_matel = .FALSE.
    CALL run_nscf(do_band, do_matel, iq)
    CALL initialize_gw()

    IF(use_symm) THEN
      WRITE(6,'("")')
      WRITE(6,'("SYMMETRIZING COULOMB Perturbations")')
      WRITE(6,'("")')
      CALL stern_symm()
    ELSE
      ngmunique = sigma_c_st%ngmt
      DO ig = 1, sigma_c_st%ngmt
         ig_unique(ig) = ig
      ENDDO
    ENDIF
!
!Need to distribute G vectors according to a sensible algorithm based
!on images which maintains, pool, and plane wave parallelism!
!      CALL distribute_pert()
!
       igstart = 1
       igstop = ngmunique
       print*, "iq, igstart, igstop"
       print*, iq, igstart, igstop
       CALL coulomb(iq, igstart, igstop, scrcoul_g)
!       CALL mp_sum ( scrcoul_g, inter_pool_comm )
100  CALL clean_pw_gw(iq)
    if(do_q0_only) GOTO 124
    if(do_epsil.and.(iq.eq.nqs)) GOTO 126
  ENDDO
124 CONTINUE
126 CONTINUE
WRITE(stdout, '("Finished Calculating Screened Coulomb")')
   DEALLOCATE( scrcoul_g )
   DEALLOCATE( ig_unique )
   DEALLOCATE( sym_ig )
   DEALLOCATE( sym_friend )
!Write W_{q}(G,G';iw) to file:
!       IF (ionode) THEN
!           CALL unfold_w(scrcoul_g,iq)
!If solve direct now need to invert epsilon:
!           IF(solve_direct.and.tinvert) CALL invert_epsilon(scrcoul_g, iq)
!           CALL davcio(scrcoul_g, lrcoul, iuncoul, iq, +1, ios)
!       ENDIF
END SUBROUTINE do_stern
