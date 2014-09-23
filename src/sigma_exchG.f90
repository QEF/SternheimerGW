SUBROUTINE sigma_exchg(ik0)
  USE kinds,         ONLY : DP
  USE gvect,         ONLY : ngm, nrxx, g, nr1, nr2, nr3, nrx1, nrx2, nrx3, nl
  USE gsmooth,       ONLY : nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nls, ngms
  USE ions_base,     ONLY : nat
  USE lsda_mod,      ONLY : nspin
  USE constants,     ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE disp,          ONLY : nqs, nq1, nq2, nq3, wq, x_q, xk_kpoints, num_k_pts
  USE control_gw,    ONLY : lgamma, eta, nbnd_occ
  USE klist,         ONLY : wk, xk, nkstot, nks
  USE io_files,      ONLY : prefix, iunigk
  USE wvfct,         ONLY : nbnd, npw, npwx, igk, g2kin, et
  USE cell_base,     ONLY : omega, tpiba2, at, bg, tpiba, alat
  USE eqv,           ONLY : evq, eprec
  USE units_gw,      ONLY : iuncoul, iungreen, iunsigma, lrsigma, lrcoul, lrgrn, iuwfc, lrwfc,&
                            iunsex, lrsex
  USE qpoint,        ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE gwsigma,       ONLY : ngmsig, nr1sex, nr2sex, nr3sex, nrsex, nlsex, ngmsex, fft6_g2r, &
                            sigma_ex, sigma_g_ex, ecutsex, ngmsco, nbnd_sig
  USE fft_scalar,    ONLY : cfft3d
  USE io_global,     ONLY : stdout, ionode_id, ionode

 
IMPLICIT NONE

!ARRAYS to describe exchange operator.
  LOGICAL  :: do_band, do_iq, setup_pw, exst, limit, single_line
  COMPLEX(DP), ALLOCATABLE :: sigma_band_ex(:,:), barcoul(:), miv(:), mvj(:)
  COMPLEX(DP), ALLOCATABLE ::  psi_ij(:,:), psi(:), dpsic(:)
  COMPLEX(DP) :: dipole(nrxxs), matel
  REAL(DP) :: rcut, spal, dvoxel, wgt, sigma_ex_tr
  INTEGER  :: ikmq, ik0, ik, igkdim
  INTEGER  :: ig, igp, npe, irr, icounter, ir, irp
  INTEGER  :: iq, ipol, ibnd, jbnd, vbnd, counter
  INTEGER  :: rec0, ios
  REAL(DP) :: qg2, qg 
  COMPLEX(DP) :: ZDOTC
  COMPLEX(DP) :: czero, exch_element
!q-vector of coulomb potential xq_coul := k_{0} - xk(ik)
  REAL(DP) :: xq_coul(3)
  REAL(DP) :: voxel
!Arrays to handle case where nlsco does not contain all G vectors required for |k+G| < ecut
  INTEGER     :: igkq_ig(npwx) 
  INTEGER     :: igkmat(npwx) 
  INTEGER     :: igkq_tmp(npwx) 
  REAL(DP)    :: sigma_ex_diag(nbnd_sig)
  INTEGER     :: iman, nman, ndeg(nbnd_sig), ideg, ikq

  CALL start_clock('sigma_exch')

if (nksq.gt.1) rewind (unit = iunigk)
!set appropriate weights for points in the brillouin zone. 
!weights of all the k-points are in odd positions in list.  
wq(:) = 0.0d0

ALLOCATE (sigma_band_ex(nbnd_sig, nbnd_sig))
ALLOCATE (psi_ij(npwx, nbnd_sig), miv(nrxxs), mvj(nrxxs))
ALLOCATE (psi(nrxxs))
ALLOCATE (dpsic(nrxxs))
ALLOCATE (barcoul(npwx))

sigma_band_ex(:,:) = (0.0d0, 0.d0)
do iq = 1, nksq
   if(lgamma) then
      wq(iq) = 0.5d0*wk(iq)
   else
      wq(iq) = 0.5d0*wk(2*iq-1) 
   endif
enddo

wgt = 1.d0/omega

if (lgamma) then
   ikq = ik0
else
   ikq = 2*ik0
endif

WRITE(6," ")
WRITE(6,'(4x,"Sigma exchange for k",i3, 3f12.7)') ik0, (xk(ipol, ikq), ipol=1,3)
WRITE(6,'(4x,"Sigma exchange for k",i3, 3f12.7)') ik0, (xk_kpoints(ipol, ik0), ipol=1,3)

czero = (0.0d0, 0.0d0)
sigma_ex(:,:) = (0.0d0, 0.0d0)
limit =.false.

DO iq = 1, nksq
    write(6,'(4x,"q ",i3, 3f12.7)') iq, (xk(ipol, iq), ipol=1,3)
!odd numbers xk
    if (lgamma) then
       ikq = iq
    else
!even k + q 
       ikq = 2*iq
    endif
!\psi_{k+q}(r)\psi_{k+q}(r')
    CALL davcio (evq, lrwfc, iuwfc, ikq, -1)
    if (nksq.gt.1) then
         read (iunigk, err = 100, iostat = ios) npw, igk
100      CALL errore ('green_linsys', 'reading igk', abs (ios) )
    endif
    if (lgamma)  npwq = npw
    if (.not.lgamma.and.nksq.gt.1) then
          read (iunigk, err = 200, iostat = ios) npwq, igkq
200       call errore ('green_linsys', 'reading igkq', abs (ios) )
    endif
    if(iq.eq.1) then
       igkmat(:) = igk(:)
       psi_ij(:,:) = evq(:,:)
       igkdim = npwq
    endif
    rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))
! q = (k0 +q) - k0
! xq_coul(:) = xk(:,ikq) - xk_kpoints(:,ik0)
! -q = k0 - (k0 + q)
   xq_coul(:) = xk_kpoints(:,ik0) - xk(:,ikq)
   do ig = 1, npwq
      qg  = sqrt((g(1,ig)  + xq_coul(1))**2.d0 + (g(2,ig) + xq_coul(2))**2.d0   &
               +  (g(3,ig ) + xq_coul(3))**2.d0)

      qg2 =     (g(1,ig)   + xq_coul(1))**2.d0 + (g(2,ig) + xq_coul(2))**2.d0   &
               + ((g(3,ig)) + xq_coul(3))**2.d0 
      if (qg < eps8) limit=.true.
      if(.not.limit) then
          spal = 1.0d0 - cos (rcut * tpiba * qg)
          barcoul(ig) = e2 * fpi / (tpiba2*qg2) * dcmplx(spal, 0.0d0)
          limit=.false.
      else
          write(6,'("Taking Limit.")')
          barcoul(ig) = (fpi*e2*(rcut**2))/2 
          limit=.false.
      endif
   enddo
!EXX
!EXX debug
   do vbnd = 1, nbnd_occ(iq)
  !do vbnd = 1, 1 
     !FFT[psi_i(ig)]
      psi (:) = (0.d0, 0.d0)
      do ig = 1, npwq
         psi (nls(igkq(ig))) = evq(ig, vbnd)
      enddo
      call cft3s (psi, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)
      do ibnd = 1, nbnd_sig
     !do ibnd = 1, 1
     !FFT[psi_i(ig)]
         dpsic(:) = (0.d0, 0.d0)
         do ig = 1, igkdim
            dpsic (nls(igkmat(ig))) = psi_ij(ig, ibnd)
         enddo
         call cft3s (dpsic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)
         dipole(:) = (0.0d0, 0.0d0)
         do ir = 1, nrxxs
            dipole (ir) = dipole (ir) + CONJG(psi (ir) ) * dpsic (ir)
         enddo
         call cft3s (dipole, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2)
         miv(:) = (0.0d0,0.0d0)
         do ig = 1, npwq
            miv(ig) = dipole(nls(ig))
         enddo
         do jbnd = 1, nbnd_sig
!FFT[psi_j(ig)]
            dpsic(:) = (0.d0, 0.d0)
            do ig = 1, igkdim
               dpsic (nls (igkmat (ig) ) ) = psi_ij (ig, jbnd)
            enddo
            call cft3s (dpsic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)
            dipole(:) = (0.0d0, 0.0d0)
            do ir = 1, nrxxs
               dipole (ir) = dipole (ir) + CONJG(psi (ir) ) * dpsic (ir)
            enddo
            call cft3s (dipole, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -1)
            mvj(:) = (0.0d0, 0.0d0)
            do ig = 1, npwq
               mvj(ig) = conjg(dipole(nl(ig)))
            enddo
            do ig = 1, npwq
                sigma_band_ex(ibnd, jbnd) = sigma_band_ex(ibnd, jbnd) - wgt*wq(iq)*miv(ig)*mvj(ig)*barcoul(ig)
            enddo 
        enddo !jbnd
      enddo !ibnd
   enddo !v\inocc
enddo ! on q


 WRITE(6,*) 
 write(stdout,'(4x,"Sigma_ex (eV)")')
 write(stdout,'(8(1x,f7.3))') real(sigma_band_ex(:,:))*RYTOEV
 write(stdout,*)
 write(stdout,'(8(1x,f7.3))') aimag(sigma_band_ex(:,:))*RYTOEV
 write(stdout,'(4x,"Sigma_ex val (eV)",8(1x,f7.2))') real(sigma_band_ex(1:nbnd_sig,1:nbnd_sig))*RYTOEV


!choose appropriate k point
  if (lgamma) then
     ikq = ik0
  else
     ikq = 2*ik0
  endif

  nman = 1
  ndeg = 1
  if(nbnd_sig.le.8) single_line=.true.
  if(nbnd_sig.gt.8) single_line=.false.

  do ibnd = 2, nbnd_sig
     if ( abs( et (ibnd, ikq) - et (ibnd-1, ikq)  ) .lt. 1.d-4 ) then
        ndeg (nman) = ndeg(nman) + 1
     else
        nman = nman + 1
     endif
  enddo


  ibnd = 0
  jbnd = 0
  do iman = 1, nman
     sigma_ex_tr = 0.0d0
     do ideg = 1, ndeg(iman)
        ibnd = ibnd + 1
        sigma_ex_tr = sigma_ex_tr + real(sigma_band_ex(ibnd,ibnd))
     enddo
     do ideg = 1, ndeg(iman)
        jbnd = jbnd + 1
        sigma_ex_diag(jbnd)= sigma_ex_tr/float(ndeg(iman))
     enddo
  enddo

  if(single_line) then 
     write(stdout,'(4x,"Sigma_ex val (eV)",8(1x,f7.2))') sigma_ex_diag(1:8)*RYTOEV
  else
     write(stdout,'(4x,"Sigma_ex val (eV)",8(1x,f7.2))', advance='no') sigma_ex_diag(1:8)*RYTOEV
  endif

  if(nbnd_sig.gt.8) then
  do ideg = 9, nbnd_sig, 8  
     if(ideg+7.lt.nbnd_sig) write(stdout,9000, advance='no') sigma_ex_diag(ideg:ideg+7)*RYTOEV
     if(ideg+7.ge.nbnd_sig) write(stdout,9000) sigma_ex_diag(ideg:nbnd_sig)*RYTOEV
  enddo
  endif

 CALL stop_clock('sigma_exch')

 9000 format(8(1x,f7.2))
 9005 format(8(1x,f14.7))

END SUBROUTINE sigma_exchg
