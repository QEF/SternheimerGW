SUBROUTINE sigma_extra(ik0)
!Sigma Extra
!now could add extra ananlytic piece: -\inf to wsigmamin-wcoul
!and from wsigmamin+wcoul to +\inf
!G(r,r';w'-w0) =  FFT[\delta(\G,G')/(w-w0-w1)]
!W(r,r';w')
!\Sigma(r, r'; w0) = \Sigma(r,r'; w0) + \int_{\-inf}^{wcoul) G(r,r'; w0 - w')W(r,\r';w')dw'
!                                     + \int_{wcoul}^{\inf) G(r,r'; w0 - w')W(r,\r';w')dw'
!
  USE kinds,                ONLY : DP
  USE klist,         ONLY : wk, xk, nkstot, nks
  USE constants,     ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE disp,          ONLY : nqs, nq1, nq2, nq3, wq, x_q, xk_kpoints, num_k_pts
  USE cell_base,     ONLY : omega, tpiba2
  USE units_gw,      ONLY : iuncoul, iungreen, iunsigma, lrsigma, lrcoul, lrgrn, iuwfc, lrwfc,iunsigext,lrsigext 
  USE gwsigma,       ONLY : ngmsco, nrsco, nlsco, fft6_g2r, ecutsco, ngmsig,&
                            nr1sco, nr2sco, nr3sco, ngmgrn, ngmpol
  USE eqv,           ONLY : evq, eprec
  USE qpoint,        ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE wvfct,         ONLY : nbnd, npw, npwx, igk, g2kin, et
  USE io_files,      ONLY : prefix, iunigk
  USE fft_scalar,    ONLY : cfft3d
  USE io_global,     ONLY : stdout, ionode_id, ionode
  USE freq_gw,       ONLY : fpol, fiu, nfs, nfsmax,&
                            nwcoul, nwgreen, nwalloc, nwsigma, wtmp, wcoul,&
                            wgreen, wsigma, wsigmamin, wsigmamax,&
                            deltaw, wcoulmax, ind_w0mw, ind_w0pw, plasmon
  USE control_gw,    ONLY : lgamma, eta
  USE gvect,         ONLY : ngm, nrxx, g, nr1, nr2, nr3, nrx1, nrx2, nrx3, nl

IMPLICIT NONE

  INTEGER :: rec0, ios
  INTEGER     :: ikq
!Integration Variable:
  COMPLEX(DP) :: cprefac
!V arrays:
  COMPLEX(DP), ALLOCATABLE :: scrcoul(:,:)
  INTEGER, ALLOCATABLE  :: gmapsym(:,:)
  COMPLEX(DP), ALLOCATABLE ::  barcoul(:,:), barcoulr(:,:)
  REAL(DP) :: rcut, spal
  REAL(DP) :: qg2, xq0s(3), qg, xxq(3)
!G arrays:
  COMPLEX(DP), ALLOCATABLE ::  green_id (:,:), green_idr(:,:)
  COMPLEX(DP), ALLOCATABLE :: sigma_g(:,:), sigma(:,:)
!Constants:
  COMPLEX(DP)  :: cone
  COMPLEX(DP)  :: ci, czero
!Counters:
  INTEGER :: ig, igp, npe, irr, icounter, ir, irp
  INTEGER :: iq, ipol, ibnd, jbnd, counter
  INTEGER :: ikmq, ik0, ik
  INTEGER     :: igkq_ig(npwx) 
  INTEGER     :: igkq_tmp(npwx) 
  COMPLEX(DP) :: aux(nrsco)
  COMPLEX(DP) :: aux1(nrsco)
  REAL(DP) :: xq_coul(3)
  REAL(DP) :: dirac, delta, x

    ALLOCATE ( barcoul     (ngmsco, ngmsco) )
    ALLOCATE ( barcoulr    (nrsco, nrsco)   )

    ALLOCATE ( green_id    (ngmgrn, ngmgrn) )
    ALLOCATE ( green_idr   (nrsco, nrsco)   )

    if(.not.allocated(sigma_g)) ALLOCATE(sigma_g (ngmsco, ngmsco))
    ALLOCATE ( sigma (nrsco, nrsco) )

    cone = DCMPLX(1.0d0,0.0d0)
    ci = (0.0d0, 1.d0)
    czero = (0.0d0, 0.0d0)
    sigma(:,:) = (0.0d0, 0.0d0)

!Need a loop to find all plane waves below ecutsco when igkq takes us outside of this sphere.  
    counter  = 0
    igkq_tmp = 0
    igkq_ig  = 0

if (nksq.gt.1) rewind (unit = iunigk)

wq(:) = 0.0d0

do iq = 1, nksq
   if(lgamma) then
      wq(iq) = 0.5d0*wk(iq)
   else
      wq(iq) = 0.5d0*wk(2*iq-1) 
   endif
enddo

if (lgamma) then
   ikq = ik0
else
   ikq = 2*ik0
endif

DO iq = 1, nksq
     if (lgamma) then
        ikq = iq
     else
! even k + q 
        ikq = 2*iq
     endif
     if (nksq.gt.1) then
          read (iunigk, err = 100, iostat = ios) npw, igk
100       CALL errore ('green_linsys', 'reading igk', abs (ios) )
     endif
     if (lgamma)  npwq = npw
     if (.not.lgamma.and.nksq.gt.1) then
           read (iunigk, err = 200, iostat = ios) npwq, igkq
200        call errore ('green_linsys', 'reading igkq', abs (ios) )
     endif
!Should have a data_type which keeps track of these indices...
    counter = 0
    igkq_tmp(:) = 0
    igkq_ig(:) = 0
    do ig = 1, npwq
       if((igkq(ig).le.ngmsco).and.((igkq(ig)).gt.0)) then
           counter = counter + 1
           igkq_tmp (counter) = igkq(ig)
           igkq_ig  (counter) = ig
       endif
    enddo

!Green \delta(G,G')=
    green_id(:,:) = (0.0d0, 0.0d0)
    do ig = 1, counter
           green_id(igkq_tmp(ig),igkq_tmp(ig)) = (1.0d0, 0.0d0)
    enddo
!   WRITE(6, '("FFT_GREEN")')
    green_idr(:,:) = czero
    do ig = 1, ngmgrn
        aux(:) = czero
        do igp = 1, ngmgrn
           aux(nlsco(igp)) = green_id(ig,igp)
        enddo
        call cfft3d (aux, nr1sco, nr2sco, nr3sco, nr1sco, nr2sco, nr3sco, +1)
        do irp = 1, nrsco
           green_idr(ig, irp) = aux(irp) / omega
        enddo
    enddo
!   the conjg/conjg is to calculate sum_G f(G) exp(-iGr)
!   following the convention set in the paper
!  [because the standard transform is sum_G f(G) exp(iGr) ]
    do irp = 1, nrsco
        aux = czero
        do ig = 1, ngmsco
           aux(nlsco(ig)) = conjg( green_idr(ig,irp) )
        enddo
        call cfft3d (aux, nr1sco, nr2sco, nr3sco, nr1sco, nr2sco, nr3sco, +1)
        green_idr(1:nrsco,irp) = conjg ( aux )
    enddo
!Now have G^A(\r,\r')...
!Calculate v(r,r'):
    rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))
    ! -q
    xq0s = (/ -0.01 , 0.00, 0.00 /) ! this should be set from input
    barcoul(:,:) = (0.0d0,0.0d0)
    ! -q = k0 - (k0 + q)
    xq_coul(:) = xk_kpoints(:,ik0) - xk(:,ikq)
    do ig = 1, ngmgrn
        qg = sqrt((g(1,ig)   + xq_coul(1))**2.d0 + (g(2,ig) + xq_coul(2))**2.d0   &
                + (g(3,ig )  + xq_coul(3))**2.d0)
        qg2 =   (g(1,ig)     + xq_coul(1))**2.d0 + (g(2,ig) + xq_coul(2))**2.d0   &
                + ((g(3,ig)) + xq_coul(3))**2.d0 

      ! These if conditions need to go...
        if (qg < eps8) qg =  sqrt((g(1, ig) + xq0s(1))**2.d0 + (g(2, ig) + xq0s(2))**2.d0 & 
                                                             + (g(3, ig) + xq0s(3))**2.d0)

        if (qg2 < eps8) qg2 =  ((g(1, ig) + xq0s(1))**2.d0   + (g(2, ig) + xq0s(2))**2.d0 & 
                                                             + (g(3, ig) + xq0s(3))**2.d0)

        spal = 1.0d0 - cos (rcut * sqrt(tpiba2) * qg)
        barcoul (ig, ig) = e2 * fpi / (tpiba2*qg2) * dcmplx(spal, 0.0d0)
    enddo

!   WRITE(6, '("FFT_COULOMB")')
    barcoulr(:,:) = (0.0d0, 0.0d0)
    do ig = 1, ngmgrn
       aux(:) = czero
       do igp = 1, ngmgrn
          aux(nlsco(igp)) = barcoul(ig,igp)
       enddo
       call cfft3d (aux, nr1sco, nr2sco, nr3sco, nr1sco, nr2sco, nr3sco, +1)
       do irp = 1, nrsco
          barcoulr(ig, irp) = aux(irp) / omega
       enddo
    enddo

! the conjg/conjg is to calculate sum_G f(G) exp(-iGr)
! following the convention set in the paper
! [because the standard transform is sum_G f(G) exp(iGr) ]

    do irp = 1, nrsco
       aux = czero
       do ig = 1, ngmgrn
          aux(nlsco(ig)) = conjg( barcoulr(ig,irp) )
       enddo
       call cfft3d (aux, nr1sco, nr2sco, nr3sco, nr1sco, nr2sco, nr3sco, +1)
       barcoulr(1:nrsco,irp) = conjg ( aux )
    enddo

!Calculate \Sigma_extra(\r,\r') = \frac{i}{2\pi}(\frac{\omega_p}{\omega_c})^2 G(\r,\r')v(\r,\r')
    !cprefac = (DCMPLX(plasmon, 0.0d0)/DCMPLX((wcoulmax), 0.0d0))**2 &
    !         * DCMPLX(0.0d0,1.0d0)/tpi
    cprefac = (DCMPLX(plasmon, 0.0d0)/DCMPLX((wcoulmax), 0.0d0))**2/tpi

   !sigma = sigma + cprefac * green_idr * barcoulr
   !sigma = sigma + cprefac * green_idr
    sigma(:,:) = sigma(:,:) + cprefac*wq(iq)*green_idr(:,:)*barcoulr(:,:)
ENDDO ! on q

    sigma_g = (0.0d0,0.0d0)
    do ir = 1, nrsco
      aux = (0.0d0, 0.0d0)
      do irp = 1, nrsco
         aux(irp) = sigma(ir,irp)
      enddo
      call cfft3d (aux, nr1sco, nr2sco, nr3sco, nr1sco, nr2sco, nr3sco, -1)
      do igp = 1, ngmsco
         sigma (ir, igp) = aux(nlsco(igp))
      enddo
    enddo

    do igp = 1, ngmsco
      aux = czero
      do ir = 1, nrsco
        aux(ir) = conjg ( sigma(ir,igp) )
      enddo
      call cfft3d (aux, nr1sco, nr2sco, nr3sco, nr1sco, nr2sco, nr3sco, -1)
      do ig = 1, ngmsco
         sigma (ig,igp) = conjg ( aux( nlsco( ig )) ) * omega
      enddo
    enddo

    do ig = ngmsco + 1, nrsco
       do igp = ngmsco + 1, nrsco
             sigma (ig, igp) = (0.0d0, 0.0d0)
       enddo
    enddo

!sigma_g = sigma(1:ngmsco,1:ngmsco,:)
    do ig = 1, ngmsco
     do igp = 1, ngmsco
        sigma_g(ig,igp)  = sigma(ig,igp)
       !if(ig.eq.igp) sigma_g(ig,igp)  = ((plasmon)/(wcoulmax))**2
     enddo 
    enddo

    CALL davcio (sigma_g, lrsigext, iunsigext, 1, 1)

    DEALLOCATE ( sigma    )
    DEALLOCATE ( sigma_g  )
    RETURN
END SUBROUTINE sigma_extra

