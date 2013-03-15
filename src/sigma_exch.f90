SUBROUTINE sigma_exch(ik0)
  USE kinds,         ONLY : DP
  USE gvect,         ONLY : ngm, nrxx, g, nr1, nr2, nr3, nrx1, nrx2, nrx3, nl
  USE gsmooth,       ONLY : nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nls, ngms
  USE ions_base,     ONLY : nat
  USE lsda_mod,      ONLY : nspin
  USE constants,     ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE disp,          ONLY : nqs, nq1, nq2, nq3, wq, x_q, xk_kpoints, num_k_pts
  USE control_gw,    ONLY : lgamma, eta
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
  LOGICAL :: do_band, do_iq, setup_pw, exst, limit
  COMPLEX(DP), ALLOCATABLE ::  greenf_na (:,:), greenf_nar(:,:)
  COMPLEX(DP), ALLOCATABLE ::  barcoul(:,:), barcoulr(:,:)
  REAL(DP) :: rcut, spal, dvoxel
  INTEGER :: ikmq, ik0, ik
  INTEGER :: ig, igp, npe, irr, icounter, ir, irp
  INTEGER :: iq, ipol, ibnd, jbnd, counter
  INTEGER :: rec0, ios
  REAL(DP) :: qg2, qg 
  COMPLEX(DP) :: ZDOTC
  COMPLEX(DP) :: czero, exch_element
  COMPLEX(DP) :: aux(nrsex)
  COMPLEX(DP) :: sigma_band_ex(nbnd_sig, nbnd_sig)
 
!q-vector of coulomb potential xq_coul := k_{0} - xk(ik)
  REAL(DP) :: xq_coul(3)
  REAL(DP) :: voxel

!Arrays to handle case where nlsco does not contain all G vectors required for |k+G| < ecut
  INTEGER     :: igkq_ig(npwx) 
  INTEGER     :: igkq_tmp(npwx) 
  INTEGER     :: ikq

! Self-Energy grid:
! iGv
  ALLOCATE ( sigma_ex    (nrsex, nrsex)   ) 

  CALL start_clock('sigma_exch')

if (nksq.gt.1) rewind (unit = iunigk)
!set appropriate weights for points in the brillouin zone. 
!weights of all the k-points are in odd positions in list.  


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

WRITE(6," ")
WRITE(6,'(4x,"Sigma exchange for k",i3, 3f12.7)') ik0, (xk(ipol, ikq), ipol=1,3)
WRITE(6,'(4x,"Sigma exchange for k",i3, 3f12.7)') ik0, (xk_kpoints(ipol, ik0), ipol=1,3)

czero = (0.0d0, 0.0d0)
sigma_ex(:,:) = (0.0d0, 0.0d0)
limit =.false.

DO iq = 1, nksq

WRITE(6,'(4x,"q ",i3, 3f12.7)') iq, (xk(ipol, iq), ipol=1,3)
! odd numbers xk
     if (lgamma) then
        ikq = iq
     else
! even k + q 
        ikq = 2*iq
     endif

!\psi_{k+q}(r)\psi_{k+q}(r')

     CALL davcio (evq, lrwfc, iuwfc, ikq, -1)

     if (nksq.gt.1) then
          read (iunigk, err = 100, iostat = ios) npw, igk
100       CALL errore ('green_linsys', 'reading igk', abs (ios) )
     endif

     if (lgamma)  npwq = npw

     if (.not.lgamma.and.nksq.gt.1) then
           read (iunigk, err = 200, iostat = ios) npwq, igkq
200        call errore ('green_linsys', 'reading igkq', abs (ios) )
     endif

!Need a loop to find all plane waves below ecutsco when igkq takes us outside of this sphere.  
     counter  = 0
     igkq_tmp = 0
     igkq_ig  = 0

!Should have a data_type which keeps track of these indices...
    !do ig = 1, npwx
    !HL be careful about going past end of npwq?
    do ig = 1, npwq
       if((igkq(ig).le.ngmsex).and.((igkq(ig)).gt.0)) then
           counter = counter + 1
           igkq_tmp (counter) = igkq(ig)
           igkq_ig  (counter) = ig
       endif
    enddo

ALLOCATE ( greenf_na   (ngmsex, ngmsex) )
    greenf_na = (0.0d0, 0.0d0)
!   psi_{k+q}(r)psi^{*}_{k+q}(r')
    do ig = 1, counter
      do igp = 1, counter
        do ibnd = 1, nbnd
           greenf_na(igkq_tmp(ig),igkq_tmp(igp)) = greenf_na(igkq_tmp(ig), igkq_tmp(igp)) + &
                                            tpi * (0.0d0, 1.0d0) * (evq(igkq_ig(ig),ibnd))* &
                                            conjg((evq(igkq_ig(igp), ibnd)))

!           greenf_na(igkq_tmp(ig),igkq_tmp(igp)) = greenf_na(igkq_tmp(ig), igkq_tmp(igp)) + &
!                                            tpi * (0.0d0, 1.0d0) * (evq(igkq_ig(ig),ibnd))* &
!                                            conjg((evq(igkq_ig(igp), ibnd)))
        enddo
      enddo
    enddo

!Fourier transform of green's function
!   WRITE(6, '("FFT_GREEN")')
ALLOCATE ( greenf_nar  (nrsex,  nrsex)  )
    greenf_nar(:,:) = czero
    do ig = 1, ngmsex
        aux(:) = czero
        do igp = 1, ngmsex
           aux(nlsex(igp)) = greenf_na(ig,igp)
        enddo
        call cfft3d (aux, nr1sex, nr2sex, nr3sex, nr1sex, nr2sex, nr3sex, +1)
        do irp = 1, nrsex
           greenf_nar(ig, irp) = aux(irp) / omega
        enddo
    enddo
DEALLOCATE(greenf_na)
!   the conjg/conjg is to calculate sum_G f(G) exp(-iGr)
!   following the convention set in the paper
!  [because the standard transform is sum_G f(G) exp(iGr) ]
    do irp = 1, nrsex
        aux = czero
        do ig = 1, ngmsex
           aux(nlsex(ig)) = conjg( greenf_nar(ig,irp) )
        enddo
        call cfft3d (aux, nr1sex, nr2sex, nr3sex, nr1sex, nr2sex, nr3sex, +1)
        greenf_nar(1:nrsex,irp) = conjg ( aux )
    enddo

ALLOCATE ( barcoul     (ngmsex, ngmsex) )

    rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))
!Sax for silicon:
!     rcut = 21.329d0
!HL using  sax cutoff
!    rcut = 0.50d0*minval(sqrt(sum(at**2,1)))*alat*tpi
!    rcut = rcut-rcut/50.0d0
! q
! q = (k0 +q) - k0
! xq_coul(:) = xk(:,ikq) - xk_kpoints(:,ik0)
    barcoul(:,:) = (0.0d0,0.0d0)
! -q = k0 - (k0 + q)
    xq_coul(:) = xk_kpoints(:,ik0) - xk(:,ikq)
    do ig = 1, ngmsex
         qg = sqrt((g(1,ig)  + xq_coul(1))**2.d0 + (g(2,ig) + xq_coul(2))**2.d0   &
                + (g(3,ig )  + xq_coul(3))**2.d0)

         qg2 =   (g(1,ig)    + xq_coul(1))**2.d0 + (g(2,ig) + xq_coul(2))**2.d0   &
                + ((g(3,ig)) + xq_coul(3))**2.d0 

        if (qg < eps8) limit=.true.
        if(.not.limit) then
            spal = 1.0d0 - cos (rcut * tpiba * qg)
            barcoul (ig, ig) = e2 * fpi / (tpiba2*qg2) * dcmplx(spal, 0.0d0)
            limit=.false.
        else
            write(6,'("Taking Limit.")')
            barcoul(ig,ig)= (fpi*e2*(rcut**2))/2 
            limit=.false.
        endif
    enddo

!   WRITE(6, '("FFT_COULOMB")')
ALLOCATE ( barcoulr    (nrsex,  nrsex)  )
    barcoulr(:,:) = (0.0d0, 0.0d0)
    do ig = 1, ngmsex
       aux(:) = czero
       do igp = 1, ngmsex
          aux(nlsex(igp)) = barcoul(ig,igp)
       enddo
       call cfft3d (aux, nr1sex, nr2sex, nr3sex, nr1sex, nr2sex, nr3sex, +1)
       do irp = 1, nrsex
          barcoulr(ig, irp) = aux(irp) / omega
       enddo
    enddo

DEALLOCATE(barcoul)
    do irp = 1, nrsex
       aux = czero
       do ig = 1, ngmsex
          aux(nlsex(ig)) = conjg( barcoulr(ig,irp) )
       enddo
       call cfft3d (aux, nr1sex, nr2sex, nr3sex, nr1sex, nr2sex, nr3sex, +1)
       barcoulr(1:nrsex,irp) = conjg ( aux )
    enddo
!might want to keep it as Gv rather than vG
       sigma_ex = sigma_ex + wq(iq)* (0.0d0,1.0d0) / tpi *  greenf_nar * barcoulr
DEALLOCATE(barcoulr)
DEALLOCATE(greenf_nar)
ENDDO ! on q

    ALLOCATE ( sigma_g_ex  (ngmsex, ngmsex) ) 
    sigma_g_ex(:,:) = (0.0d0,0.0d0)

!take matrix element here:
!The results did not differ at all from taking matrix elements in G space for silicon at Gamma. 
!Still useful to have for later purposes.(26 lines)
!if (nksq.gt.1) rewind (unit = iunigk)
!CALL davcio (evq, lrwfc, iuwfc, 1, -1)
!     if (nksq.gt.1) then
!          read (iunigk, err = 100, iostat = ios) npw, igk
!     endif
!     if (lgamma)  npwq = npw
!do ibnd = 1, nbnd
!   exch_element = DCMPLX(0.0d0, 0.0d0)
!   aux = czero
!   dvoxel = (omega/nrsex)**2
!   do ig = 1, npwq
!      if((igkq(ig).le.ngmsex).and.((igkq(ig)).gt.0)) then
!          aux(nlsex(ig)) = evq(igkq(ig), ibnd)  
!      endif
!   enddo
!   call cfft3d (aux, nr1sex, nr2sex, nr3sex, nr1sex, nr2sex, nr3sex, +1)
!   aux = (aux/sqrt(omega))
!   do irp = 1, nrsex
!      do ir = 1, nrsex
!         exch_element = exch_element + (conjg(aux(ir))*sigma_ex(ir,irp)*aux(irp))*dvoxel
!      enddo
!   enddo
!   write (6, '("Exchange Element:")')
!   write (6,*) exch_element*RYTOEV
!   write (6,*) Omega
!enddo
!real space exchange element.

!WRITE(6, '("FFT_SIGMA")')
    CALL sigma_r2g_ex(sigma_ex, sigma_g_ex)
!Storing full sigma_ex operator in G  space...
    CALL davcio(sigma_g_ex, lrsex, iunsex, 1, 1)
    CALL stop_clock('sigma_exch')

    DEALLOCATE(sigma_g_ex)
    DEALLOCATE (sigma_ex)
END SUBROUTINE sigma_exch
