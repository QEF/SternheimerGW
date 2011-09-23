SUBROUTINE gw_product(ik0) 

! G TIMES W PRODUCT
  USE kinds,         ONLY : DP
  USE gvect,         ONLY : ngm, nrxx, g, nr1, nr2, nr3, nrx1, nrx2, nrx3, nl
  USE gsmooth,       ONLY : nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nls, ngms
  USE ions_base,     ONLY : nat
  USE lsda_mod,      ONLY : nspin
  USE constants,     ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE disp,          ONLY : nqs, nq1, nq2, nq3, wq
  USE control_gw,    ONLY : lgamma, eta
  USE klist,         ONLY : wk, xk
  USE io_files,      ONLY : prefix, iunigk
  USE wvfct,         ONLY : nbnd, npw, npwx, igk, g2kin, et
  USE cell_base,     ONLY : omega, tpiba2
  USE eqv,           ONLY : evq, eprec
  USE freq_gw,       ONLY : fpol, fiu, nfs, nfsmax,&
                         nwcoul, nwgreen, nwalloc, nwsigma, wtmp, wcoul,&
                         wgreen, wsigma, wsigmamin, wsigmamax,&
                         deltaw, wcoulmax, ind_w0mw, ind_w0pw
  USE units_gw,      ONLY : iuncoul, iungreen, iunsigma, lrsigma, lrcoul, lrgrn, iuwfc, lrwfc
  USE gwsigma,       ONLY : ngmsig, sigma, sigma_g, nrsig
  USE qpoint,        ONLY : xq, npwq, igkq, nksq, ikks, ikqs



  IMPLICIT NONE

!For running PWSCF need some variables 
  LOGICAL :: do_band, do_iq, setup_pw, exst
 
!Pade arrays
  COMPLEX(DP), ALLOCATABLE :: z(:), u(:), a(:)

!W arrays 
  COMPLEX(DP), ALLOCATABLE :: scrcoul_g (:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: scrcoul_pade_g (:,:)
  COMPLEX(DP), ALLOCATABLE :: scrcoul(:,:)

!G arrays
!Analytic component

  COMPLEX(DP), ALLOCATABLE :: greenf_g(:,:), greenfp(:,:), greenfm(:,:)

!Non-analytic component
  COMPLEX(DP), ALLOCATABLE ::  sigma_ex(:,:), greenf_na (:,:), greenf_nar(:,:), barcoul(:,:), barcoulr(:,:), gr_N(:,:)
  REAL(DP) :: qg2, xq0s(3), qg, xxq(3)
  COMPLEX(DP) :: cprefac 


  !various frequency counters 
  INTEGER :: iwim, iw 
  INTEGER  :: iw0, iwp, iw0mw, iw0pw
  REAL(DP) w_ryd(nwcoul)
  REAL(DP) :: rcut, spal
  COMPLEX(DP) :: ci
  INTEGER :: ig, igp, npe, irr, icounter
  INTEGER :: igstart, igstop, igpert
  INTEGER :: iq, ipol, ibnd
  INTEGER :: ikmq, ik0, ik
  INTEGER :: rec0, ios

  !CHECK FOR NAN's
  REAL(DP) :: ar, ai


! HL Need to think about how all of this is going to be parallelized. 
! #ifdef __PARA
!       scrcoul_g = czero
!       if (me.eq.1.and.mypool.eq.1) then
! #endif
! #ifdef __PARA
!       endif
!       use poolreduce to broadcast the results to every pool
!       call poolreduce ( 2 * ngms * ngms * nwim, scrcoul_g)
! #endif
! Real Space represenation of sigma the nrsig, nlsig() values/arrays need
! to be generated in ggensigma.f90.

     ALLOCATE ( sigma_ex (nrsig, nrsig) )
     ALLOCATE ( sigma_g(ngmsig, ngmsig, nwsigma ))
     ALLOCATE ( sigma (nrsig, nrsig, nwsigma, 1))

   !W arrays
    !HL since in coulomb I only need to calculate the 11-12 frequencies
    !and then do the analytic continuation I should only allocate an array of the size nfs to the
    !frequency dimension.
     ALLOCATE ( scrcoul_g (ngmsig, ngmsig, nfs, 1) )
     ALLOCATE ( scrcoul_pade_g (ngmsig, ngmsig) )

!Probably should define W(r,r') in real space.
     ALLOCATE ( scrcoul (nrsig, nrsig) )

! G arrays
     ALLOCATE ( greenf_g   (ngmsig, ngmsig)      )
     ALLOCATE ( greenfp    (nrsig, nrsig)        )
     ALLOCATE ( greenfm    (nrsig, nrsig)        )
     ALLOCATE ( greenf_na  (nrsig,nrsig)         )
     ALLOCATE ( greenf_nar (nrsig,nrsig)         )
     ALLOCATE ( gr_N       (nrsig,nrsig)         ) 

! for G^{NA} 
     ALLOCATE ( barcoul    (nrsig,nrsig) )
     ALLOCATE ( barcoulr   (nrsig,nrsig) )

!Also locally allocating arrays which will contain the pade approximants.
     ALLOCATE  (z(nfs), a(nfs))

! Array for coulomb frequencies.
     w_ryd(:) = wcoul(:)/RYTOEV

! Array for proper FFT correspondence.
! CALCULATE SIGMA^{C}(r,r',w) = \int G(r,r',w + w')(W(r,r',w') - v(r,r')) dw'

! WRITE(6,'(" wq ")')
! WRITE(6,*) wq(:) 

   WRITE(6," ")
   WRITE(6,'(4x,"Direct product GW for k0(",i3," ) = (",3f12.7," )")') ik0, (xk(ipol, ik0), ipol=1,3)
   WRITE(6," ")

!now sum over {q} the products G(k0-q)W(q).
!HL ci is i! 

     ci = (0.0d0, 1.d0)
     sigma(:,:,:,:) = (0.0d0, 0.0d0)

CALL start_clock('gwproduct')

DO iq = 1, nqs 
 WRITE(6, '("Summing q = ", i4 )') iq


! In parallel case this is zeroed  scrcoul_g = (0.0d0, 0.0d0)
! read Pade coefficients of screened coulomb interaction (W-v) and broadcast
! read ( iuncoul, rec = iq, iostat = ios) scrcoul_g
! #ifdef __PARA
!       scrcoul_g = czero
!       if (me.eq.1.and.mypool.eq.1) then
! #endif
! HL the length of record of the coulomb file is: 
! lrcoul = 2 * ngmsig * ngmsig * nfs.
! What is currently written to it comes from coulomb.f90 and consists of:
! the co-efficients of the pade approximant. 
!Generate proper shift indices for FFTs on different points in the B.Z:      

     CALL davcio(scrcoul_g, lrcoul, iuncoul, iq, -1 )

!#ifdef __PARA
!      endif
!      use poolreduce to broadcast the results to every pool
!      call poolreduce ( 2 * ngms * ngms * nwim, scrcoul_g)
!#endif

!Start integration over iw +/- wcoul. 
!All G-vectors from coulomb should have gamma ordering

     DO iw = 1, nwcoul
        do ig = 1, ngmsig
          do igp = 1, ngmsig

            do iwim = 1, nfs
               z(iwim) = dcmplx( 0.d0, fiu(iwim))
               a(iwim) = scrcoul_g (ig,igp,iwim,1)
               ! write(6,*)a(iwim)
            enddo

            do iwim = 1, nfs
               ar = real(a(iwim))
               ai = aimag(a(iwim))

               if ( ( ar .ne. ar ) .or. ( ai .ne. ai ) ) then
                     !write(6,*) (z(i),i=1,N)
                     !write(6,*) (u(i),i=1,N)
                     !write(6,*) (a(i),i=1,N)
                      a(:) = (0.0d0, 0.0d0)
                      write (6,'("pade-coeffs nan ", 3i4)')ig, igp, iq 
               endif
            enddo
   !Use Vidberg Serene to generate value of screened coulomb at all iw -wcoul to +wcoul.
             call pade_eval ( nfs, z, a, dcmplx( w_ryd(iw), eta), scrcoul_pade_g (ig,igp))

   !write(6,*) scrcoul_pade_g(ig,igp) 
          enddo
        enddo

        ! WRITE(6,*) w_ryd(iw), scrcoul_pade_g (:,:) 
        ! 2. fft the scrcoulomb to real space.
        ! write(6,'("starting fft6")')

        call fft6_g2r ( scrcoul_pade_g, scrcoul)

     !HL need to think about the weight of each of these q-points. Currently hard wire to 1/16
     !ci = (0.d0, 1.d0)

        !cprefac = (deltaw/RYTOEV) * (1.0d0/216.0d0) * (0.0d0, 1.0d0)/ tpi
        !cprefac = (deltaw/RYTOEV) * wq(iq) *(1.0d0/216.0d0) * (0.0d0, 1.0d0)/ tpi

        cprefac = (deltaw/RYTOEV) * wq(iq) * (0.0d0, 1.0d0)/ tpi

        do iw0 = 1, nwsigma
            iw0mw = ind_w0mw (iw0,iw)
            iw0pw = ind_w0pw (iw0,iw)
!
!generate green's function in real space for iw0mw (greenfm)
!PARALLEL RELATED STUFF  
!       #ifdef __PARA
!             greenf_g = czero
!             if (me.eq.1.and.mypool.eq.1) then
!       #endif
!           rec0 = (iw0pw-1) * nk0 * nq + (ik0-1) * nq + (iq-1) + 1

            rec0 = (iw0mw-1) * 1 * nqs + (ik0-1) * nqs + (iq-1) + 1
            CALL davcio( greenf_g, lrgrn, iungreen, rec0, -1 )

!#ifdef __PARA
!          endif
!          use poolreduce to broadcast the results to every pool
!          call poolreduce ( 2 * ngms * ngms, greenf_g )
!#endif

           CALL fft6_g2r ( greenf_g, greenfm )

! generate green's function in real space for iw0pw (greenfp).
! HL more parallel
!         #ifdef __PARA
!                 greenf_g = czero
!                 if (me.eq.1.and.mypool.eq.1) then
!         #endif
!         rec0 = (iw0pw-1) * nk0 * nq + (ik0-1) * nq + (iq-1) + 1

          rec0 = (iw0pw-1) * 1 * nqs + (ik0-1) * nqs + (iq-1) + 1
          CALL davcio( greenf_g, lrgrn, iungreen, rec0, -1 )

!         read ( iungreen, rec = rec0, iostat = ios) greenf_g
!#ifdef __PARA
!         endif
!         use poolreduce to broadcast the results to every pool
!         call poolreduce ( 2 * ngms * ngms, greenf_g )
!#endif
!         greenf_g is ngms*ngms, greenf is nrs*nrs

          call fft6_g2r ( greenf_g, greenfp )

!@HL Might want to think about a slightly more pleasant way to do this integral.
! Rex (the man himself) suggested it should be possible to do the convolution on the imaginary axis 
! And Then analytically continue the self-energy to the real axis. 

      sigma (:,:,iw0,1) = sigma (:,:,iw0,1) + cprefac * (greenfp(:,:) + greenfm(:,:)) * scrcoul(:,:)
        enddo !on iw0  
     ENDDO ! on frequency convolution over w'
     CALL clean_pw_gw(iq)

ENDDO ! end loop on {k0-q} and {q} 

!EXCHANGE PART OF THE SELF-ENERGY 
!Here I need to figure out the easy way of regenerating necessary wave functions.
!There are also loads of double fourier transforms so the code corresponds with the convention in the paper.

sigma_ex(:,:) = (0.0d0, 0.0d0)

DO iq = 1, nqs
 WRITE(6, '("Exchange = ", i4 )') iq

! NON-ANALYTIC PART OF THE GREEN'S FUNCTION
! HL Need the wave functions at k0mq = xk0(:,ik0) - xq(:,iq) again. 
! These will have all already been calculated

      CALL prepare_kmq(do_band, do_iq, setup_pw, iq, ik0)
      CALL run_pwscf(do_band)
      CALL initialize_gw()

!HL: Again we need to make sure the gamma wave function is what's 
!    read in when we do the calculation at q = 0.

     if(lgamma) then
        CALL davcio (evq, lrwfc, iuwfc, 1, -1)
     else
        CALL davcio (evq, lrwfc, iuwfc, 2, -1)
     endif
     
      greenf_na = (0.0d0, 0.0d0)

      if (lgamma) write(6, '("lgamma=.true.")')
      if (nksq.gt.1) rewind (unit = iunigk)


      if (nksq.gt.1) then
           read (iunigk, err = 100, iostat = ios) npw, igk
 100        call errore ('green_linsys', 'reading igk', abs (ios) )
        endif


        if (.not.lgamma.and.nksq.gt.1) then
           read (iunigk, err = 200, iostat = ios) npwq, igkq
 200        call errore ('green_linsys', 'reading igkq', abs (ios) )
        endif

        if (lgamma) npwq = npw 
 
!@HL again parallelism across Gs... 
!do ig = igstart, igstop

WRITE(6,'("NBND", i4)')nbnd
     do ig = 1, ngmsig

!      do ig = 1, npwq
!       do igp = 1, ngmsig
!HL     Same deal here. igkq has max dim npwx, but can contain G-vectors up to ecutwfc.
!        do igp = 1, npwq
!          do ibnd = 1, nbnd
!!         no spin factor here! in <nk|Sigma|nk> only states of the same spin couple [cf HL86]
!           if(ig.and.igkq(igp).lt.ngmsig.and.igkq(igp).gt.0) then
!              greenf_na (ig,igkq(igp)) = greenf_na (ig,igkq(igp)) + &
!                          tpi * (0.0d0, 1.0d0) * (evq(ig,ibnd)) * conjg(evq(igp,ibnd))
!            endif
!          enddo
!        enddo
!      enddo

        do igp = 1, npwq
          do ibnd = 1, nbnd
!         no spin factor here! in <nk|Sigma|nk> only states of the same spin couple [cf HL86]
          !if(igkq(ig).and.igkq(igp).lt.ngmsig.and.igkq(igp).and.igkq(ig).gt.0) then
            if( ((igkq(igp).lt.ngmsig) .and. (igkq(ig)).lt.ngmsig) .and. ((igkq(igp).gt.0) & 
               .and. (igkq(ig)).gt.0)) then
              greenf_na (igkq(ig),igkq(igp)) = greenf_na (igkq(ig),igkq(igp)) + &
                          tpi * (0.0d0, 1.0d0) * (evq(ig,ibnd)) * conjg(evq(igp,ibnd))
            endif
          enddo
        enddo
      enddo


!#ifdef __PARA
!   use poolreduce to bring together the results from each pool
!   call poolreduce ( 2 * nrs * nrs, greenf_na)
!#endif
! In SGWI fft6_g2r is used sporadically. I will use it consistently here (much tidier)
! also should write fft6_r2g again much (tidier). Might require a look at NumRec FFT
! refresher course. Anyway now we have greenf_na, and greenf_nar and one's in G space
! one is in (r)eal space.  

      call fft6_g2r(greenf_na, greenf_nar)

! Now Sigma^{ex}(r,r') = - \sum_{v}\psi^{*}_{v}(r)\psi_{v}(r')v(r,r'). 
! Requires spencer/alavi and then another 6-d fourier transform.
! COULOMB EXCHANGE TERM

      rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))
      xq0s = (/ 0.01 , 0.00, 0.00 /) ! this should be set from input
      barcoul(:,:) = (0.0d0,0.0d0)
      barcoulr(:,:) = (0.0d0, 0.0d0)
      do ig = 1, ngmsig
        qg2 = (g(1,ig )+xq(1))**2.d0 + (g(2,ig )+xq(2))**2.d0 + (g(3,ig )+ xq(3))**2.d0 
        if (qg2 < eps8) qg2 =  (g(1,ig )+xq0s(1))**2.d0 + (g(2,ig )+xq0s(2))**2.d0 & 
                                                          + (g(3,ig )+ xq0s(3))**2.d0
        spal = 1.0d0 - cos ( rcut * sqrt(tpiba2) * sqrt(qg2) )
!       Looks like we are only taking diagonal elements G=G'. Which makes sense for the bare coulomb
        barcoul (ig,ig) = dcmplx(e2 * fpi / (tpiba2*qg2), 0.0d0) * dcmplx(spal, 0.0d0)
      enddo


! #ifdef __PARA
!      use poolreduce to bring together the results from each pool
!      call poolreduce ( 2 * nrs * nrs, barcoul)
! #endif


     call fft6_g2r(barcoul, barcoulr)

!HL sigma_ex(:,:) = sigma_ex(:,:) + wq (iq) * ci / twopi * greenf_na(:,:) * barcoul(:,:)
!@HL Here sigma_ex is in real space.

!     sigma_ex(:,:) = sigma_ex(:,:) + wq(iq)*(1.0d0/216.0d0) * ci / tpi * greenf_nar(:,:) * barcoulr(:,:)

!If wq has been normalized to one...

     sigma_ex(:,:) = sigma_ex(:,:) + wq(iq)* ci / tpi * greenf_nar(:,:) * barcoulr(:,:)
     CALL clean_pw_gw(iq)
ENDDO ! on q

   WRITE(6," ")
   WRITE(6,'(4x,"FORMING FULL Sigma")')
!
!  \Sigma = \Sigma^c + \Sigma^ex
!
   do iw = 1, nwsigma
        sigma(:,:,iw,1) = sigma(:,:,iw,1) + sigma_ex(:,:)

       !sigma(:,:,iw,1) = sigma(:,:,iw,1) 
   enddo

   WRITE(6,'(4x,"Fourier Transform of Sigma")')

   CALL sigma_r2g(sigma,sigma_g) 

!  Now write Sigma in G space to file. 
!  WRITE ( iunsigma, rec = ik0, iostat = ios) sigma_g

   WRITE(6,'(4x,"Writing Sigma to File")')

   CALL davcio (sigma_g, lrsigma, iunsigma, ik0, 1)
   CALL stop_clock('gwproduct') 

!HL should de-allocate lots of things here. 
   DEALLOCATE (barcoul, barcoulr)
   DEALLOCATE (greenf_g)
   DEALLOCATE (greenfp, greenfm)
   DEALLOCATE (greenf_na, greenf_nar)
   DEALLOCATE (scrcoul_g, scrcoul_pade_g, scrcoul)
   DEALLOCATE (sigma_ex, sigma_g, sigma)
   DEALLOCATE (z, a)
 
   RETURN
END SUBROUTINE gw_product

