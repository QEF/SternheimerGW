SUBROUTINE pade_eps()
  USE io_global,     ONLY : stdout, ionode_id, ionode
  USE kinds,         ONLY : DP
  USE units_gw,      ONLY : iuncoul, iungreen, iunsigma, lrsigma, lrcoul, lrgrn, iuwfc, lrwfc
  USE freq_gw,       ONLY : fpol, fiu, nfs, nfsmax
  USE constants,     ONLY : e2, fpi, RYTOEV, tpi, pi
  USE gwsigma,       ONLY : ngmsco, sigma, sigma_g, nrsco, nlsco, fft6_g2r, ecutsco, ngmsig,&
                            nr1sco, nr2sco, nr3sco, ngmgrn, ngmpol
  USE qpoint,        ONLY : xq
  USE gvect,         ONLY : ngm, nrxx, g, nr1, nr2, nr3, nrx1, nrx2, nrx3, nl
  USE disp,            ONLY : x_q, done_iq, rep_iq, done_rep_iq, comp_iq,&
                              xk_kpoints

IMPLICIT NONE

  COMPLEX(DP), ALLOCATABLE :: scrcoul_g (:,:,:)
  COMPLEX(DP), ALLOCATABLE :: scrcoul_pade_g (:,:)
  COMPLEX(DP), ALLOCATABLE :: z(:), u(:), a(:)

!For Coulomb grid
  REAL(DP)                  ::   wcoulmax, wcoulmin, deltaws
  INTEGER                   ::   nwcoul
  REAL(DP), ALLOCATABLE     ::   wcoul(:), w_ryd(:)

  REAL(DP)                 :: eta
  INTEGER :: ig, igp, iqrec, iw, iwim, counter

  ALLOCATE  (z(nfs), a(nfs), u(nfs))
  ALLOCATE ( scrcoul_g      (ngmpol, ngmpol, nfs)     )
  ALLOCATE ( scrcoul_pade_g (ngmpol, ngmpol)          )

  scrcoul_g(:,:,:)   = (0.0d0, 0.0d0)
  scrcoul_pade_g(:,:) = (0.0d0, 0.0d0)

!Choose q-point you wish to analytically continue
!epsilon_{\q}(\G,\G'):

   counter = 0

   eta = 0.007

   wcoulmin   =   0.00
   wcoulmax   = 100.00
   deltaws    =   0.2
   nwcoul  = 1 + ceiling((wcoulmax - wcoulmin)/deltaws)

   allocate (wcoul(nwcoul))
   allocate (w_ryd(nwcoul))

   do iw = 1, nwcoul
      wcoul(iw) = wcoulmin + (wcoulmax-wcoulmin)/float(nwcoul-1)*float(iw-1)
   enddo
   w_ryd = wcoul/RYTOEV
   do iqrec = 1,6,3

      scrcoul_g(:,:,:)   = (0.0d0, 0.0d0)
      scrcoul_pade_g(:,:) = (0.0d0, 0.0d0)
      CALL davcio(scrcoul_g, lrcoul, iuncoul, iqrec, -1)
    do ig    = 1,8,4
!     do igp   = 1,8,4
      igp   = ig
     write(6,'("All good so far.")')
!    write(6, '("#iq ", 3f12.7, "ig ", 3f12.7, "igp ", 3f12.7)') x_q(1:3,iqrec), g(1:3,ig), g(1:3,igp)
!    write(200+counter, '("#iq ", 3f12.7, "ig ", 3f12.7, "igp ", 3f12.7)') x_q(1:3,iqrec), g(1:3,ig), g(1:3,igp)
        do iw = 1, nfs
           z(iw) = fiu(iw)
           u(iw) = scrcoul_g(ig, igp, iw)
           if (ig.eq.igp) then
              write(200+counter,'(4f15.8)') z(iw)*RYTOEV, 1.0d0 + scrcoul_g(ig,igp,iw)
           else
              write(200+counter,'(4f15.8)') z(iw)*RYTOEV, scrcoul_g(ig,igp,iw)
           endif
        enddo
        write(6,'("All good so far.")')
        call pade_coeff ( nfs, z, u, a)
!    write(100+counter,'("#iq", 3f12.7, "ig", 3f12.7, "igp", 3f12.7)') x_q(1:3,iqrec), g(1:3,ig), g(1:3,igp)
        do iw = 1, nwcoul
           call pade_eval (nfs, z, a, dcmplx(w_ryd(iw), eta), scrcoul_pade_g (ig,igp))
           if (ig.eq.igp) then
              write(100+counter,'(3f15.8)') wcoul(iw), 1.00+scrcoul_pade_g(ig,igp)
           else
              write(100+counter,'(3f15.8)') wcoul(iw), scrcoul_pade_g(ig,igp)
           endif
        enddo
        counter=counter+1
!!     end do!igp
    end do!ig
   end do!iqrec
end subroutine pade_eps
