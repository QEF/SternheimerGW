SUBROUTINE pade_eps()
  USE io_global,     ONLY : stdout, ionode_id, ionode
  USE kinds,         ONLY : DP
  USE units_gw,      ONLY : iuncoul, iungreen, iunsigma, lrsigma, lrcoul, lrgrn, iuwfc, lrwfc
  USE freq_gw,       ONLY : fpol, fiu, nfs, nfsmax
  USE constants,     ONLY : e2, fpi, RYTOEV, tpi, pi
  USE gwsigma,       ONLY : ngmsco, sigma, sigma_g, nrsco, nlsco, fft6_g2r, ecutsco, ngmsig,&
                            nr1sco, nr2sco, nr3sco, ngmgrn, ngmpol

IMPLICIT NONE

  COMPLEX(DP), ALLOCATABLE :: scrcoul_g (:,:,:)
  COMPLEX(DP), ALLOCATABLE :: scrcoul_pade_g (:,:)
  COMPLEX(DP), ALLOCATABLE :: z(:), u(:), a(:)

!For Coulomb grid
  REAL(DP)                  ::   wcoulmax, wcoulmin, deltaws
  INTEGER                   ::   nwcoul
  REAL(DP), ALLOCATABLE     ::   wcoul(:), w_ryd(:)

  REAL(DP)                 :: eta
  INTEGER :: ig, igp, iqrec, iw,iwim

  ALLOCATE  (z(nfs), a(nfs), u(nfs))
  ALLOCATE ( scrcoul_g      (ngmpol, ngmpol, nfs)     )
  ALLOCATE ( scrcoul_pade_g (ngmpol, ngmpol)          )

  scrcoul_g(:,:,:)   = (0.0d0, 0.0d0)
  scrcoul_pade_g(:,:) = (0.0d0, 0.0d0)

!Choose q-point you wish to analytically continue
!epsilon_{\q}(\G,\G'):

   iqrec = 1
   ig    = 10
   igp   = 10

   eta = 0.01

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

   CALL davcio(scrcoul_g, lrcoul, iuncoul, iqrec, -1)

   do iw = 1, nfs
      z(iw) = fiu(iw)
      u(iw) = scrcoul_g (ig, igp, iw)
   enddo

   call pade_coeff ( nfs, z, u, a)

   do iw = 1, nfs 
      scrcoul_g (ig, igp, iw) = a(iw)
   enddo

   do iwim = 1, nfs
      z(iwim) = fiu(iwim)
      a(iwim) = scrcoul_g (ig,igp,iwim)
   enddo

   do iw = 1, nwcoul
      call pade_eval (nfs, z, a, dcmplx(w_ryd(iw), eta), scrcoul_pade_g (ig,igp))
      if (ig.eq.igp) then
         write(100,*) wcoul(iw), eta, scrcoul_pade_g(ig,igp)
         write(stdout,'(3f15.8)') wcoul(iw), 1.00+scrcoul_pade_g(ig,igp)
      else
         write(100,*) wcoul(iw), eta, scrcoul_pade_g(ig,igp)
         write(stdout,'(3f15.8)') wcoul(iw), scrcoul_pade_g(ig,igp)
      endif
   enddo

!Could put in a plotting function:
!i.e. time dependence of the plasmon modes!
end subroutine pade_eps
