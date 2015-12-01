subroutine sigprod(isymop, dz, scrcoul_g, greenf_g, sigma_g, gmapsym)
  USE kinds,          ONLY : DP
  USE cell_base,      ONLY : tpiba2, tpiba, omega, alat, at
  USE gwsigma,        ONLY : sigma_c_st, gcutcorr
  USE fft_interfaces, ONLY : invfft, fwfft
  USE fft_custom,     ONLY : fft_cus, set_custom_grid, ggent, gvec_init
  USE gvect,          ONLY : ngm
  USE symm_base,      ONLY : nsym, s, time_reversal, t_rev, ftau, invs, nrot

  IMPLICIT NONE

  COMPLEX(DP)              :: greenfr   (sigma_c_st%nnr), scrcoulr(sigma_c_st%nnr)
  COMPLEX(DP)              :: greenf_g  (sigma_c_st%ngmt, sigma_c_st%ngmt, 2*nwcoul)
  COMPLEX(DP)              :: scrcoul_g (sigma_c_st%ngmt, sigma_c_st%ngmt, nfs)
  COMPLEX(DP)              :: sigma_g   (sigma_c_st%ngmt, sigma_c_st%ngmt)
  COMPLEX(DP)              :: dz
  INTEGER                  :: gmapsym  (ngm, nrot) 
  INTEGER                  :: ig, igp
!G_{1}
    do ig = 1, gcutcorr
       greenfr  (sigma_c_st%nlt(1:gcutcorr)) = (0.0d0, 0.0d0)
       scrcoulr (sigma_c_st%nlt(1:gcutcorr)) = (0.0d0, 0.0d0)
       greenfr  (sigma_c_st%nlt(1:gcutcorr)) =   greenf_g(ig,1:gcutcorr) 
       scrcoulr (sigma_c_st%nlt(gmapsym(1:gcutcorr,isymop))) =  scrcoul_g(ig,1:gcutcorr) 
       call invfft('Custom', greenfr, sigma_c_st%dfftt)
       call invfft('Custom', scrcoulr, sigma_c_st%dfftt)
       aux(:) = dz*greenfr(:)*scrcoulr(:)
       call fwfft('Custom', aux, sigma_c_st%dfftt)
       sigma_g(ig, 1:gcutcorr) = sigma_g(ig, 1:gcutcorr) +  aux(sigma_c_st%nlt(1:gcutcorr))
    enddo
!G_{2}
    do igp = 1, gcutcorr
       greenfr  (sigma_c_st%nlt(1:gcutcorr)) = (0.0d0, 0.0d0)
       scrcoulr (sigma_c_st%nlt(1:gcutcorr)) = (0.0d0, 0.0d0)
       greenfr  (sigma_c_st%nlt(1:gcutcorr)) = greenf_g(1:gcutcorr, igp) 
       scrcoulr (sigma_c_st%nlt(gmapsym(1:gcutcorr,isymop))) =  scrcoul_g(ig,1:gcutcorr) 
       call invfft('Custom', greenfr, sigma_c_st%dfftt)
       call invfft('Custom', scrcoulr, sigma_c_st%dfftt)
       aux(:) = dz*greenfr(:)*scrcoulr(:)
       call fwfft('Custom', aux, sigma_c_st%dfftt)
       sigma_g(1:gcutcorr, igp) = sigma_g(1:gcutcorr, igp) +  aux(sigma_c_st%nlt(1:gcutcorr))
    enddo
end subroutine sigprod
