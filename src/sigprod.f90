subroutine sigprod(isymop, dz, scrcoul_g, greenf_g, sigma_g, gmapsym)
  use kinds,          only : DP
  use cell_base,      only : tpiba2, tpiba, omega, alat, at
  use gwsigma,        only : sigma_c_st, gcutcorr
  use fft_interfaces, only : invfft, fwfft
  use fft_custom,     only : fft_cus, set_custom_grid, ggent, gvec_init
  use gvect,          only : nl, ngm, g, nlm, gstart, gl, igtongl
  use symm_base,      only : nsym, s, time_reversal, t_rev, ftau, invs, nrot
  use mp_world,       only : nproc, mpime

  implicit none

  complex(DP)              :: dz
  complex(DP)              :: sigma_r   (sigma_c_st%dfftt%nnr, gcutcorr)
  complex(DP)              :: greenf_r  (sigma_c_st%dfftt%nnr, gcutcorr)
  complex(DP)              :: scrcoul_r (sigma_c_st%dfftt%nnr, gcutcorr)
  complex(DP)              :: greenfr   (sigma_c_st%dfftt%nnr)
  complex(DP)              :: scrcoulr  (sigma_c_st%dfftt%nnr)
  complex(DP)              :: aux       (sigma_c_st%dfftt%nnr)
  complex(DP)              :: greenf_g  (gcutcorr, gcutcorr)
  complex(DP)              :: scrcoul_g (gcutcorr, gcutcorr)
  complex(DP)              :: sigma_g   (gcutcorr, gcutcorr)
  integer                  :: gmapsym(ngm, nrot), isymop, isym 
  integer                  :: ig, igp, ir
  integer                  :: ig1, ig1p, igw
  real(DP), parameter      :: eps=1.0d-5

! G_{1}
  do igp = 1, gcutcorr
     scrcoulr(:) = (0.0d0,0.0d0)
     greenfr(:)  = (0.0d0,0.d0)
     !scrcoulr (sigma_c_st%nlt(gmapsym(1:gcutcorr,invs(isymop)))) = conjg(scrcoul_g(1:gcutcorr,igp))
     scrcoulr (sigma_c_st%nlt(gmapsym(1:gcutcorr,isymop))) = conjg(scrcoul_g(1:gcutcorr,igp))
     !scrcoulr (sigma_c_st%nlt(1:gcutcorr)) = conjg(scrcoul_g(1:gcutcorr,igp))
     greenfr  (sigma_c_st%nlt(1:gcutcorr)) = conjg(greenf_g(1:gcutcorr,igp))
     call invfft('Custom', scrcoulr, sigma_c_st%dfftt)
     call invfft('Custom', greenfr, sigma_c_st%dfftt)
     scrcoul_r(1:sigma_c_st%dfftt%nnr, igp)  = conjg(scrcoulr(1:sigma_c_st%dfftt%nnr))
     greenf_r(1:sigma_c_st%dfftt%nnr, igp)   = conjg(greenfr(1:sigma_c_st%dfftt%nnr))
  enddo
  do ir = 1, sigma_c_st%dfftt%nnr
     scrcoulr(:) = (0.0d0,0.0d0)
     greenfr(:)  = (0.0d0,0.d0)
     aux(:)      = (0.0d0,0.d0)
     !scrcoulr(sigma_c_st%nlt(gmapsym(1:gcutcorr,invs(isymop)))) = scrcoul_r(ir, 1:gcutcorr)
     scrcoulr(sigma_c_st%nlt(gmapsym(1:gcutcorr,isymop))) = scrcoul_r(ir, 1:gcutcorr)
     !scrcoulr(sigma_c_st%nlt(1:gcutcorr)) = scrcoul_r(ir, 1:gcutcorr)
     greenfr (sigma_c_st%nlt(1:gcutcorr)) = greenf_r (ir, 1:gcutcorr)
     call invfft('Custom', scrcoulr, sigma_c_st%dfftt)
     call invfft('Custom', greenfr, sigma_c_st%dfftt)
     aux         = dz*greenfr(:)*scrcoulr(:)/omega
     call fwfft('Custom', aux, sigma_c_st%dfftt)
     greenfr(:)  = (0.0d0,0.d0)
     greenfr(1:gcutcorr)  = aux(sigma_c_st%nlt(1:gcutcorr))
     sigma_r(ir, 1:gcutcorr) = greenfr(1:gcutcorr)
  enddo
  do igp = 1, gcutcorr
     aux(:)  = (0.0d0,0.d0)
     aux(:) = conjg(sigma_r(:, igp))
     call fwfft('Custom', aux, sigma_c_st%dfftt)
     sigma_r(1:gcutcorr, igp) = conjg(aux(sigma_c_st%nlt(1:gcutcorr)))
  enddo
  sigma_g(1:gcutcorr,1:gcutcorr) = sigma_g(1:gcutcorr,1:gcutcorr) + sigma_r(1:gcutcorr,1:gcutcorr)
end subroutine sigprod
