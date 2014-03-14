SUBROUTINE invert_epsilon(scrcoul_g_in, iq)
USE kinds,         ONLY : DP
USE gwsigma,       ONLY : ngmsco, sigma, sigma_g, nrsco, nlsco, fft6_g2r, ecutsco, ngmpol
USE freq_gw,       ONLY : fpol, fiu, nfs, nfsmax, nwcoul, wcoul

IMPLICIT NONE    

COMPLEX(DP)       :: scrcoul_g_in(ngmpol, ngmpol, nfs, 1)
COMPLEX(DP)       :: work(ngmpol)
INTEGER           :: ig, igp, npe, irr, icounter, ir, irp
INTEGER           :: isym, iwim, iq, iw
INTEGER           :: iwork(ngmpol), info

!Need block inversion routine if iq is gamma.
do iw = 1, nfs
 ! call ZHETRF ('U', ngmpol, scrcoul_g_in(1:ngmpol,1:ngmpol,iw,1), ngmpol, iwork, work, ngmpol, info)
   call ZGETRF (ngmpol, ngmpol, scrcoul_g_in(1:ngmpol,1:ngmpol,iw,1), ngmpol, iwork, info)
   call errore ('invert epsilon', 'factorization', info)
 ! call ZHETRI ('U', ngmpol, scrcoul_g_in(1:ngmpol,1:ngmpol,iw,1), ngmpol, iwork, work, info)
   call ZGETRI (ngmpol, scrcoul_g_in(1:ngmpol,1:ngmpol,iw,1), ngmpol, iwork, work, ngmpol, info)
   call errore ('invert epsilon', 'inversion', info)
enddo

write(6,*)
write(6,'("Done epsilon inversion.")') 


!at Gamma wings of W are 0.
if(iq.eq.1) then
do iw = 1, nfs
    do ig = 2, ngmpol
       scrcoul_g_in(ig,1,iw,1) = dcmplx(0.0d0,0.0d0)
    enddo
    do igp = 2, ngmpol
       scrcoul_g_in(1,igp,iw,1) = dcmplx(0.0d0,0.0d0)
    enddo
enddo
endif


!do iw=1,nfs
!    write(6,'(15f12.7)') real(scrcoul_g_in(1:15,1:15,iw,1))
!    print*,""
!enddo

!We store epsilon-1 to disk:
do iw = 1, nfs
   do ig = 1, ngmpol
      scrcoul_g_in(ig,ig,iw,1) = scrcoul_g_in(ig,ig,iw,1) - dcmplx(1.0d0,0.0d0)
   enddo
enddo

END SUBROUTINE invert_epsilon
