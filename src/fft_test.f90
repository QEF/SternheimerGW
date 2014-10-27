SUBROUTINE fft_test(ik0)

  USE kinds,         ONLY : DP
  USE gwsigma,       ONLY : ngmsco, sigma, sigma_g, nrsco, nlsco, fft6_g2r, ecutsco, ngmsig
  USE wvfct,         ONLY : nbnd, npw, npwx, igk, g2kin, et
  USE eqv,           ONLY : evq, eprec
  USE gvect,         ONLY : g, ngm, ecutwfc, nl
  USE cell_base,     ONLY : tpiba2, tpiba
  USE qpoint,        ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE disp,          ONLY : nqs, nq1, nq2, nq3, wq, x_q 
  USE control_gw,    ONLY : lgamma, eta
  USE units_gw,      ONLY : iuncoul, iungreen, iunsigma, lrsigma, lrcoul, lrgrn, iuwfc, lrwfc
  USE freq_gw,       ONLY : fpol, fiu, nfs, nfsmax,&
                            nwcoul, nwgreen, nwalloc, nwsigma, wtmp, wcoul,&
                            wgreen, wsigma, wsigmamin, wsigmamax,&
                            deltaw, wcoulmax, ind_w0mw, ind_w0pw

  IMPLICIT NONE
 
  INTEGER :: ig, igp, npe, irr, icounter, ir, irp
  INTEGER :: igstart, igstop, igpert
  INTEGER :: iq, ipol
  INTEGER :: ikmq, ik0, ik
  INTEGER :: rec0, ios
  INTEGER :: counter

!G arrays:
  COMPLEX(DP), ALLOCATABLE :: greenf_g(:,:), greenfp(:,:), greenfm(:,:)
  INTEGER  :: iwim, iw 
  INTEGER  :: iw0, iwp, iw0mw, iw0pw
  REAL(DP) :: w_ryd(nwcoul)
  
!For running PWSCF need some variables 
  LOGICAL :: do_band, do_iq, setup_pw, exst

  ALLOCATE (igkq  ( npwx))    
  ALLOCATE (g2kin ( npwx))    

  ALLOCATE ( greenf_g       (ngmsco, ngmsco)          )
  ALLOCATE ( greenfp        (nrsco, nrsco)            )
  ALLOCATE ( greenfm        (nrsco, nrsco)            )

  igkq(:) = 0 

  CALL allocate_fft()
  CALL ggen()

  write(6,*)
  write(6,*)nlsco
  write(6,*)

   do iq = 1, nqs
      CALL prepare_kmq(do_band, do_iq, setup_pw, iq, ik0)
      CALL gk_sort(xq(1), ngm, g, (4*ecutsco)/tpiba2, npwq, igkq, g2kin)

      counter = 0
      do ig = 1, npwx
         if((igkq(ig).lt.ngmsco).and.((igkq(ig)).gt.0)) then
        !if((g2kin(ig).le.(4*ecutsco/tpiba2)).and.((igkq(ig)).gt.0)) then
             counter = counter + 1
             igkq(counter) = igkq(ig)
         endif
      enddo
      
      write(6,'("igkq")') 
      write(6,*) igkq
      write(6,*) 
      write(6,'("nlsco(igkq)")') 

      write(6,*) counter
      do ig = 1, counter, 6
        !write(6,*) nlsco(igkq(ig)) 
         write(6,'(6i)') nlsco(igkq(ig:ig+5)) 
      enddo

      write(6,*) 

     !do iw = 1, nwgreen
     !DO iw0 = 1, nwsigma

      iw = 1
      iw0 = 1
      iw0mw = ind_w0mw (iw0,iw)
      iw0pw = ind_w0pw (iw0,iw)

     !rec0 = (iw0mw-1) * 1 * nqs + (ik0-1) * nqs + (iq-1) + 1
     !CALL davcio( greenf_g, lrgrn, iungreen, rec0, -1 )
     !CALL fft6_g2r (counter, nrsco, nlsco, greenf_g, greenfm, 1, 1)
     !write(6,*) greenfm
     !ENDDO !on iw0  
   enddo

END SUBROUTINE fft_test
