SUBROUTINE invert_epsilon(scrcoul_g_in, iq)
USE kinds,         ONLY : DP
USE gwsigma,       ONLY : ngmsco, sigma, sigma_g, nrsco, nlsco, fft6_g2r, ecutsco, ngmpol
USE freq_gw,       ONLY : fpol, fiu, nfs, nfsmax, nwcoul, wcoul

IMPLICIT NONE    

COMPLEX(DP)       :: scrcoul_g_in(ngmpol, ngmpol, nfs, 1)
INTEGER           :: ig, igp, npe, irr, icounter, ir, irp
INTEGER           :: isym, iwim, iq, iw
INTEGER           :: iwork(ngmpol), info
complex(kind(DP)) :: work(ngmpol)
!NOTE: the polar and w operators are HERMITEAN when calculated for imaginary frequencies
!Need block inversion routine if iq is gamma.
!The factorisation seems to go screwy here if we enforce the hermiticity
!going to try it with a generalized inversion routine in case the matrix is no longer hermitian.
!might have to enforce this...

do iw = 1, nfs
  call ZHETRF ('U', ngmpol, scrcoul_g_in(:,:,iw,1), ngmpol, iwork, work, ngmpol, info)
  call errore ('invert epsilon', 'factorization', info)

  call ZHETRI ('U', ngmpol, scrcoul_g_in(:,:,iw,1), ngmpol, iwork, work, info)
  call errore ('invert epsilon', 'inversion', info)

  write(6,*)
  do ig = 1, 14
        write(6,'(14f14.7)')real(scrcoul_g_in(ig,1:14,1,1))
  enddo
enddo

END SUBROUTINE invert_epsilon
