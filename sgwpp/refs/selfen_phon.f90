 !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE selfen_phon
  !-----------------------------------------------------------------------
  !
  !  compute the imaginary part of the phonon self energy due to electron-
  !  phonon interaction in the Migdal approximation. This corresponds to 
  !  the phonon linewidth (half width). The phonon frequency is taken into
  !  account in the energy selection rule.
  !
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation
  !
  !  RM 24/02/2014
  !  removed the calculation of transport lambda
  !  that part was wrong since vkq corresponded to dmef(:,:,:,ikq) 
  !  at the last q point calculated and not the current q 
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : at, bg
  USE io_global, ONLY : stdout
  USE phcom,     ONLY : lgamma, nmodes
  USE epwcom,    ONLY : nbndsub, lrepmatf, iunepmatf, fsthick, &
                        eptemp, ngaussw, degaussw, iuetf,     &
                        wmin, wmax, nw, nbndskip, a2f, epf_mem, etf_mem, &
                        nsmear, delta_smear, eig_read, eps_acustic, &
                        rand_k, nkf1, nkf2, nkf3, efermi_read, fermi_energy
  USE pwcom,     ONLY : nelec, ef, isk, nbnd
  USE el_phon,   ONLY : epf17, ibndmax, ibndmin, etf, &
                        etfq, wkf, xqf, wqf, nksf, nxqf,   &
                        nksqf, wf, nkstotf, xkf, xqf, lambda_all
  USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, pi, ci, cone, czero, zero
#ifdef __PARA
  USE mp,        ONLY : mp_barrier,mp_sum
  USE mp_global, ONLY : me_pool,inter_pool_comm,my_pool_id
#endif
  !
  implicit none
  !
  integer :: ik, ikk, ikq, ibnd, jbnd, imode, nrec, iq, fermicount, ismear
  complex(kind=DP) epf(ibndmax-ibndmin+1, ibndmax-ibndmin+1)
  real(kind=DP) :: g2, ekk, ekq, wq, ef0, wgkk, wgkq, weight, dosef, w0g1, w0g2,    &
                   degaussw0, eptemp0, lambda_tot
  !
  real(kind=DP), external :: efermig, dos_ef, w0gauss, wgauss
  real(kind=DP) :: gamma(nmodes)
  !
  WRITE(6,'(/5x,a)') repeat('=',67)
  WRITE(6,'(5x,"Phonon (Imaginary) Self-Energy in the Migdal Approximation")') 
  WRITE(6,'(5x,a/)') repeat('=',67)
  !
  IF ( fsthick .lt. 1.d3 ) &
     WRITE(stdout, '(/5x,a,f10.6,a)' ) &
     'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
  WRITE(stdout, '(/5x,a,f10.6,a)' ) &
     'Golden Rule strictly enforced with T = ',eptemp(1) * ryd2ev, ' eV'
  !
  IF ( .not. ALLOCATED (lambda_all) )  ALLOCATE( lambda_all  (nmodes, nxqf, nsmear) )
  lambda_all(:,:,:)   = zero
  !
  ! here we loop on smearing values JN - this may be an option later on 
  !
  DO ismear = 1, nsmear
     !
     degaussw0 = (ismear-1) * delta_smear + degaussw
     eptemp0 = (ismear-1) * delta_smear + eptemp(1)
     !
     ! Fermi level and corresponding DOS
     !
     IF ( efermi_read ) THEN
        !
        ef0 = fermi_energy
        !
     ELSE
        !
        ef0 = efermig(etf,nbndsub,nksf,nelec,wkf,degaussw0,ngaussw,0,isk)
        ! if some bands are skipped (nbndskip.neq.0), nelec has already been recalculated 
        ! in ephwann_shuffle
        !
     ENDIF
     !
     dosef = dos_ef (ngaussw, degaussw0, ef0, etf, wkf, nksf, nbndsub)
     !   N(Ef) in the equation for lambda is the DOS per spin
     dosef = dosef / two
     !
     WRITE(6, 100) degaussw0 * ryd2ev, ngaussw
     WRITE(6, 101) dosef / ryd2ev, ef0 * ryd2ev
     !
     ! loop over all q points of the fine mesh (this is in the k-para case 
     ! it should always be k-para for selfen_phon)
     ! 
     DO iq = 1, nxqf
        !
        CALL start_clock('PH SELF-ENERGY')
        !
        fermicount = 0
        gamma(:) = zero
        !
        DO ik = 1, nksqf
           !
           IF (lgamma) THEN
              ikk = ik
              ikq = ik
           ELSE
              ikk = 2 * ik - 1
              ikq = ikk + 1
           ENDIF
           !
           ! we read the hamiltonian eigenvalues (those at k+q depend on q!) 
           !
           IF (etf_mem) THEN
              etf (ibndmin:ibndmax, ikk) = etfq (ibndmin:ibndmax, ikk, iq)
              etf (ibndmin:ibndmax, ikq) = etfq (ibndmin:ibndmax, ikq, iq)
           ELSE
              nrec = (iq-1) * nksf + ikk
              CALL davcio ( etf (ibndmin:ibndmax, ikk), ibndmax-ibndmin+1, iuetf, nrec, - 1)
              nrec = (iq-1) * nksf + ikq
              CALL davcio ( etf (ibndmin:ibndmax, ikq), ibndmax-ibndmin+1, iuetf, nrec, - 1)
           ENDIF
           !
           ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
           IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .AND. &
                ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) THEN
              !
              fermicount = fermicount + 1
              !
              DO imode = 1, nmodes
                 !
                 ! the phonon frequency
                 wq = wf (imode, iq)
                 !
                 !  we read the e-p matrix from disk / memory
                 !
                 IF (etf_mem) then
                    epf(:,:) = epf17 ( ik, iq, :, :, imode)
                 ELSE
                    nrec = (iq-1) * nmodes * nksqf + (imode-1) * nksqf + ik
                    CALL dasmio ( epf, ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, -1)
                 ENDIF
                 !
                 DO ibnd = 1, ibndmax-ibndmin+1
                    !
                    !  the fermi occupation for k
                    ekk = etf (ibndmin-1+ibnd, ikk) - ef0
                    wgkk = wgauss( -ekk/eptemp0, -99)
                    !
                    DO jbnd = 1, ibndmax-ibndmin+1
                       !
                       !  the fermi occupation for k+q
                       ekq = etf (ibndmin-1+jbnd, ikq) - ef0
                       wgkq = wgauss( -ekq/eptemp0, -99)  
                       !
                       ! here we take into account the zero-point sqrt(hbar/2M\omega)
                       ! with hbar = 1 and M already contained in the eigenmodes
                       ! g2 is Ry^2, wkf must already account for the spin factor
                       !
                       ! NON FUNZIONA SCAMBIANDO i,j
                       ! the coupling from Gamma acoustic phonons is negligible
                       IF ( wq .gt. eps_acustic ) THEN
                          g2 = abs(epf (jbnd, ibnd))**two / ( two * wq )
                       ELSE
                          g2 = 0.d0
                       ENDIF
                       !
                       ! = k-point weight * [f(E_k) - f(E_k+q)]/ [E_k+q - E_k -w_q +id]
                       ! This is the imaginary part of the phonon self-energy, sans the matrix elements
                       !
                       !weight = wkf (ikk) * (wgkk - wgkq) * &
                       !   aimag ( cone / ( ekq - ekk - wq - ci * degaussw0 ) ) 
                       !
                       ! the below expression is positive-definite, but also an approximation
                       ! which neglects some fine features
                       !
                       w0g1 = w0gauss ( ekk / degaussw0, 0) / degaussw0
                       w0g2 = w0gauss ( ekq / degaussw0, 0) / degaussw0
                       weight = pi * wq * wkf (ikk) * w0g1 * w0g2
                       !
                       gamma   ( imode ) = gamma   ( imode ) + weight * g2
                       !
                    ENDDO ! jbnd
                    !
                 ENDDO   ! ibnd
                 !
              ENDDO ! loop on q-modes
              !
           ENDIF ! endif fsthick
           !
           CALL stop_clock('PH SELF-ENERGY')
           !
        ENDDO ! loop on k
#ifdef __PARA
        !
        ! collect contributions from all pools (sum over k-points)
        ! this finishes the integral over the BZ  (k)
        !
        CALL mp_sum(gamma,inter_pool_comm) 
        CALL mp_sum(fermicount, inter_pool_comm)
        CALL mp_barrier()
        !
#endif
        !
        WRITE(6,'(/5x,"iq = ",i5," coord.: ", 3f9.5, " wt: ", f9.5)') iq, xqf(:,iq), wqf(iq)
        WRITE(6,'(5x,a)') repeat('-',67)
        !
        lambda_tot = zero
        DO imode = 1, nmodes
           ! 
           wq = wf (imode, iq)
           IF ( wq .gt. eps_acustic ) THEN 
              lambda_all( imode, iq, ismear )  = gamma( imode ) / pi / wq**two / dosef
           ENDIF
           lambda_tot = lambda_tot + lambda_all( imode, iq, ismear )
           !
           WRITE(6, 102) imode, lambda_all(imode,iq,ismear), ryd2mev*gamma(imode), ryd2mev*wq
  CALL flush(6)
        ENDDO
        !
        WRITE(6,103) lambda_tot
        WRITE(6,'(5x,a/)') repeat('-',67)
        !
        ! test ONLY
#ifdef __PARA
        IF (me_pool == 0) &
#endif
        !     
        WRITE( stdout, '(/5x,a,i8,a,i8/)' ) &
             'Number of (k,k+q) pairs on the Fermi surface: ',fermicount, ' out of ', nkstotf/2
        !
     ENDDO ! loop on q
     !
  ENDDO !smears
  !
  ! generate the Eliashberg spectral function
  !
  IF (a2f) call eliashberg_a2f
  !
100 FORMAT(5x,'Gaussian Broadening: ',f10.6,' eV, ngauss=',i4)
101 FORMAT(5x,'DOS =',f10.6,' states/spin/eV/Unit Cell at Ef=',f10.6,' eV')
102 FORMAT(5x,'lambda( ',i3,' )=',f12.4,'   gamma=',f12.4,' meV','   omega=',f12.4,' meV')
103 FORMAT(5x,'lambda( tot )=',f12.4)
  !
  RETURN
  !
  END SUBROUTINE selfen_phon
  !
  !-----------------------------------------------------------------------
  SUBROUTINE selfen_phon_fly (iq )
  !-----------------------------------------------------------------------
  !
  !  compute the imaginary part of the phonon self energy due to electron-
  !  phonon interaction in the Migdal approximation. This corresponds to 
  !  the phonon linewidth (half width). The phonon frequency is taken into
  !  account in the energy selection rule.
  !
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation.  This routine is similar to the one above
  !  but it is ONLY called from within ephwann_shuffle and calculates 
  !  the selfenergy for one phonon at a time.  Much smaller footprint on
  !  the disk
  !
  !  RM 24/02/2014
  !  redefined the size of coskkq, vkk, vkq within the fermi windwow
  !  cleaned up the subroutine
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : at, bg
  USE io_global, ONLY : stdout
  USE io_files,  ONLY : find_free_unit, prefix, tmp_dir
  use phcom,     ONLY : lgamma, nmodes
  use epwcom,    ONLY : nbndsub, lrepmatf, iunepmatf, fsthick, &
                        eptemp, ngaussw, degaussw, iuetf,     &
                        wmin, wmax, nw, nbndskip, epf_mem, etf_mem, &
                        nsmear, delta_smear, eig_read, eps_acustic, &
                        rand_k, nkf1, nkf2, nkf3, efermi_read, fermi_energy
  use pwcom,     ONLY : nelec, ef, isk, nbnd
  use el_phon,   ONLY : epf17, ibndmax, ibndmin, etf, &
                        etfq, wkf, xqf, wqf, nksf, nxqf,   &
                        nksqf, wf, nkstotf, xkf, xqf, &
                        lambda_all, lambda_v_all, nrr_k, &
                        dmef, ndegen_k, irvec
  USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero, pi, ci, cone, czero
#ifdef __PARA
  use mp,        ONLY : mp_barrier,mp_sum
  use mp_global, ONLY : me_pool,inter_pool_comm,my_pool_id,npool
#endif
  !
  implicit none
  !
  integer :: ik, ikk, ikq, ibnd, jbnd, imode, nrec, iq, nq, fermicount, ismear
  complex(kind=DP) epf (ibndmax-ibndmin+1, ibndmax-ibndmin+1)
  real(kind=DP) :: g2, ekk, ekq, wq, ef0, wgkk, wgkq, &
                   weight, dosef, w0g1, w0g2, degaussw0, eptemp0, lambda_tot
  !
  real(kind=DP), external :: efermig, dos_ef, w0gauss, wgauss
  real(kind=DP) :: gamma(nmodes),gamma_v(nmodes)
  real(kind=DP) :: gamma_all(nmodes,nxqf,nsmear), gamma_v_all(nmodes,nxqf,nsmear)
  real(kind=DP) :: coskkq(ibndmax-ibndmin+1, ibndmax-ibndmin+1)
  real(kind=DP) :: DDOT, vkk(3,ibndmax-ibndmin+1), vkq(3,ibndmax-ibndmin+1)
  !
  !
  IF ( iq .eq. 1 ) THEN 
     WRITE(6,'(/5x,a)') repeat('=',67)
     WRITE(6,'(5x,"Phonon (Imaginary) Self-Energy in the Migdal Approximation (on the fly)")') 
     WRITE(6,'(5x,a/)') repeat('=',67)
     !
     IF ( fsthick.lt.1.d3 ) &
          WRITE(stdout, '(/5x,a,f10.6,a)' ) &
          'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
     WRITE(stdout, '(/5x,a,f10.6,a)' ) &
          'Golden Rule strictly enforced with T = ',eptemp(1) * ryd2ev, ' eV'
     !
     IF ( .not. ALLOCATED (lambda_all) )    ALLOCATE( lambda_all  (nmodes, nxqf, nsmear) )
     IF ( .not. ALLOCATED (lambda_v_all) )  ALLOCATE( lambda_v_all(nmodes, nxqf, nsmear) )
     lambda_all(:,:,:)   = zero
     lambda_v_all(:,:,:) = zero
     !
  ENDIF
  !
  DO ismear = 1, nsmear
     !
     degaussw0 = (ismear-1) * delta_smear + degaussw
     eptemp0   = (ismear-1) * delta_smear + eptemp(1)
     !
     ! Fermi level and corresponding DOS
     !
     IF ( efermi_read ) THEN
        !
        ef0 = fermi_energy
        !
     ELSE
        !
        ef0 = efermig(etf,nbndsub,nksf,nelec,wkf,degaussw0,ngaussw,0,isk)
        ! if some bands are skipped (nbndskip.neq.0), nelec has already been recalculated 
        ! in ephwann_shuffle
        !
     ENDIF
     !
     dosef = dos_ef (ngaussw, degaussw0, ef0, etf, wkf, nksf, nbndsub)
     !   N(Ef) in the equation for lambda is the DOS per spin
     dosef = dosef / two
     !
     IF ( iq .eq. 1 ) THEN 
        WRITE (6, 100) degaussw0 * ryd2ev, ngaussw
        WRITE (6, 101) dosef / ryd2ev, ef0 * ryd2ev
        !WRITE (6, 101) dosef / ryd2ev, ef  * ryd2ev
     ENDIF
     !
     CALL start_clock('PH SELF-ENERGY')
     !
     fermicount = 0
     gamma(:)   = zero
     gamma_v(:) = zero
     !
     DO ik = 1, nksqf
        !
        IF (lgamma) THEN
           ikk = ik
           ikq = ik
        ELSE
           ikk = 2 * ik - 1
           ikq = ikk + 1
        ENDIF
        ! 
        coskkq = 0.d0
        DO ibnd = 1, ibndmax-ibndmin+1
           DO jbnd = 1, ibndmax-ibndmin+1
              ! coskkq = (vk dot vkq) / |vk|^2  appears in Grimvall 8.20
              ! this is different from :   coskkq = (vk dot vkq) / |vk||vkq|
              ! In principle the only coskkq contributing to lambda_tr are both near the
              ! Fermi surface and the magnitudes will not differ greatly between vk and vkq
              ! we may implement the approximation to the angle between k and k+q vectors also 
              ! listed in Grimvall
              !
              ! v_(k,i) = 1/m <ki|p|ki> = 2 * dmef (:, i,i,k)
              ! 1/m  = 2 in Rydberg atomic units
              !
              vkk(:, ibnd ) = 2.0 * REAL (dmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk ) )
              vkq(:, jbnd ) = 2.0 * REAL (dmef (:, ibndmin-1+jbnd, ibndmin-1+jbnd, ikq ) )
              if ( abs ( vkk(1,ibnd)**2 + vkk(2,ibnd)**2 + vkk(3,ibnd)**2) .gt. 1.d-4) &
                   coskkq(ibnd, jbnd ) = DDOT(3, vkk(:,ibnd ), 1, vkq(:,jbnd),1)  / &
                   DDOT(3, vkk(:,ibnd), 1, vkk(:,ibnd),1)
           ENDDO
        ENDDO
        !
        ! we read the hamiltonian eigenvalues (those at k+q depend on q!) 
        !
        ! when we see references to iq for file readinq, it is always = 1 for on the fly calculations
        IF (etf_mem) then
           etf (ibndmin:ibndmax, ikk) = etfq (ibndmin:ibndmax, ikk,  1)
           etf (ibndmin:ibndmax, ikq) = etfq (ibndmin:ibndmax, ikq,  1)
        ELSE
           nrec = ikk
           CALL davcio ( etf (ibndmin:ibndmax, ikk), ibndmax-ibndmin+1, iuetf, nrec, - 1)
           nrec = ikq
           CALL davcio ( etf (ibndmin:ibndmax, ikq), ibndmax-ibndmin+1, iuetf, nrec, - 1)
        ENDIF
        !
        ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
        IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .AND. &
             ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) THEN
           !
           fermicount = fermicount + 1
           !
           DO imode = 1, nmodes
              !
              ! the phonon frequency
              wq = wf (imode, iq)
              !
              !  we read the e-p matrix from disk / memory
              !
              IF (etf_mem) then
                 epf(:,:) = epf17 ( ik,  1, :, :, imode)
              ELSE
                 nrec = (imode-1) * nksqf + ik
                 CALL dasmio ( epf, ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, -1)
              ENDIF
              !
              DO ibnd = 1, ibndmax-ibndmin+1
                 !
                 !  the fermi occupation for k
                 ekk = etf (ibndmin-1+ibnd, ikk) - ef0
                 wgkk = wgauss( -ekk/eptemp0, -99)
                 !
                 DO jbnd = 1, ibndmax-ibndmin+1
                    !
                    !  the fermi occupation for k+q
                    ekq = etf (ibndmin-1+jbnd, ikq) - ef0
                    wgkq = wgauss( -ekq/eptemp0, -99)  
                    !
                    ! here we take into account the zero-point sqrt(hbar/2M\omega)
                    ! with hbar = 1 and M already contained in the eigenmodes
                    ! g2 is Ry^2, wkf must already account for the spin factor
                    !
                    ! NON FUNZIONA SCAMBIANDO i,j
                    ! the coupling from Gamma acoustic phonons is negligible
                    IF ( wq .gt. eps_acustic ) THEN
                       g2 = abs(epf (jbnd, ibnd))**two / ( two * wq )
                    ELSE
                       g2 = 0.d0
                    ENDIF
                    !
                    ! = k-point weight * [f(E_k) - f(E_k+q)]/ [E_k+q - E_k -w_q +id]
                    ! This is the imaginary part of the phonon self-energy, sans the matrix elements
                    !
                    !weight = wkf (ikk) * (wgkk - wgkq) * &
                    !     aimag ( cone / ( ekq - ekk - wq - ci * degaussw0 ) ) 
                    !
                    ! the below expression is positive-definite, but also an approximation
                    ! which neglects some fine features
                    !
                    w0g1 = w0gauss ( ekk / degaussw0, 0) / degaussw0
                    w0g2 = w0gauss ( ekq / degaussw0, 0) / degaussw0
                    weight = pi * wq * wkf (ikk) * w0g1 * w0g2
                    !
                    gamma   (imode) = gamma   (imode) + weight * g2 
                    gamma_v (imode) = gamma_v (imode) + weight * g2 * (1-coskkq(ibnd, jbnd) ) 
                    !
                 ENDDO ! jbnd
                 !
              ENDDO   ! ibnd
              !
           ENDDO ! loop on q-modes
           !
        ENDIF ! endif fsthick
        !
        CALL stop_clock('PH SELF-ENERGY')
        !
     ENDDO ! loop on k
     !
#ifdef __PARA
     !
     ! collect contributions from all pools (sum over k-points)
     ! this finishes the integral over the BZ  (k)
     !
     CALL mp_sum(gamma,inter_pool_comm) 
     CALL mp_sum(gamma_v,inter_pool_comm) 
     CALL mp_sum(fermicount, inter_pool_comm)
     CALL mp_barrier()
     !
#endif
     !
     WRITE(6,'(/5x,"ismear = ",i5," iq = ",i5," coord.: ", 3f9.5, " wt: ", f9.5)') ismear, iq, xqf(:,iq), wqf(iq)
     WRITE(6,'(5x,a)') repeat('-',67)
     !
     lambda_tot = 0
     !
     DO imode = 1, nmodes
        ! 
        wq = wf (imode, iq)
        IF ( wq .gt. eps_acustic ) THEN 
           lambda_all  ( imode, iq, ismear ) = gamma  ( imode ) / pi / wq**two / dosef
           lambda_v_all( imode, iq, ismear ) = gamma_v( imode ) / pi / wq**two / dosef
        ENDIF
        gamma_all  ( imode, iq, ismear ) = gamma  ( imode )
        gamma_v_all( imode, iq, ismear ) = gamma_v( imode )
        lambda_tot = lambda_tot + lambda_all( imode, iq, ismear )
        !
        WRITE(6, 102) imode, lambda_all(imode,iq,ismear), ryd2mev*gamma_all(imode,iq,ismear), ryd2mev*wq
        !
     ENDDO
     !
     WRITE(6, 103) lambda_tot
     WRITE(6,'(5x,a/)') repeat('-',67)
     ! 
     ! test ONLY
#ifdef __PARA
     IF (me_pool == 0) &
#endif
     !     
     WRITE( stdout, '(/5x,a,i8,a,i8/)' ) &
           'Number of (k,k+q) pairs on the Fermi surface: ',fermicount, ' out of ', nkstotf/2
     !
  ENDDO !smears
  !
100 FORMAT(5x,'Gaussian Broadening: ',f10.6,' eV, ngauss=',i4)
101 FORMAT(5x,'DOS =',f10.6,' states/spin/eV/Unit Cell at Ef=',f10.6,' eV')
102 FORMAT(5x,'lambda( ',i3,' )=',f12.4,'   gamma=',f12.4,' meV','   omega=',f12.4,' meV')
103 FORMAT(5x,'lambda( tot )=',f12.4)
  !
  RETURN
  !
END SUBROUTINE selfen_phon_fly
!

