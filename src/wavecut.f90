subroutine wavecut(npwq, counter, igk, igkq_tmp, igkq_ig)

  USE gwsigma,       ONLY : sigma_x_st, nbnd_sig
  USE wvfct,         ONLY : nbnd, npw, npwx, g2kin, et, ecutwfc

IMPLICIT NONE

  INTEGER,  INTENT(INOUT)   :: counter
  INTEGER,  INTENT(INOUT)  :: igkq_ig(npwx) 
  INTEGER,  INTENT(INOUT)  :: igkq_tmp(npwx) 
  INTEGER,  INTENT(IN)     :: igk(npwx)
  INTEGER   :: ig, npwq  

  counter  = 0
  igkq_tmp(:) = 0
  igkq_ig(:)  = 0
  do ig = 1, npwq
     if((igk(ig).le.sigma_x_st%ngmt).and.((igk(ig)).gt.0)) then
               counter = counter + 1
               igkq_tmp (counter) = igk(ig)
               igkq_ig  (counter) = ig
     endif
  enddo

end subroutine wavecut
