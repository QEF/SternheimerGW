! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!-----------------------------------------------------------------------
PROGRAM gw
  !-----------------------------------------------------------------------
  ! ... This is the main driver of the GW code.
  ! ... It reads all the quantities calculated by pwscf, it
  ! ... checks if some recover file is present and determines
  ! ... which calculation needs to be done. Finally, it makes
  ! ... a loop over the q points. At a generic q, if necessary, it
  ! ... recalculates the band structure by calling pwscf again.
  ! NC = norm conserving pseudopotentials
  ! US = ultrasoft pseudopotentials
  ! PAW = projector augmented-wave
  ! [1] LDA, [2] [1]+GGA, [3] [2]+LSDA/sGGA, [4] [3]+Spin-orbit/nonmagnetic,
  ! [5] [4]+Spin-orbit/magnetic

  USE io_global,        ONLY : stdout, ionode_id, ionode
  USE wvfct,            ONLY : nbnd,npwx
  USE disp,             ONLY : nqs, num_k_pts, xk_kpoints, w_of_q_start
  USE output,           ONLY : fildrho
  USE check_stop,       ONLY : check_stop_init
  USE gw_restart,       ONLY : gw_writefile, destroy_status_run
  USE save_gw,          ONLY : clean_input_variables

  USE mp_global,        ONLY: mp_startup, nimage, npool, intra_image_comm, inter_image_comm, &
                              nproc_pool, mpime, nproc, my_pool_id, me_pool, &
                              mp_global_end, inter_pool_comm

  USE mp,                 ONLY: mp_barrier, mp_bcast, mp_sum, mp_end
  USE parallel_include,   ONLY: mpi_comm_world 

  USE path_io_routines,   ONLY : io_path_start
  USE control_gw,         ONLY : done_bands, reduce_io, recover, tmp_dir_gw, &
                               ext_restart, bands_computed, bands_computed, nbnd_occ, lgamma,&
                               do_coulomb, do_sigma_c, do_sigma_exx, do_green, do_sigma_matel,&
                               do_q0_only, multishift, do_sigma_extra, solve_direct, tinvert

  USE input_parameters, ONLY : pseudo_dir
  USE io_files,         ONLY : prefix, tmp_dir
  USE control_flags,    ONLY : restart
  USE qpoint,           ONLY : xq
  USE save_gw,          ONLY : tmp_dir_save
  USE environment,      ONLY: environment_start
  USE freq_gw,          ONLY : nfs, nwsigma
  USE units_gw,         ONLY : iuncoul, iungreen, lrgrn, lrcoul, iunsigma, lrsigma, lrsex, iunsex,&
                               iunresid, lrresid, iunalphabeta, lralphabeta, iunsigext, lrsigext
  USE basis,            ONLY : starting_wfc, starting_pot, startingconfig
  USE gwsigma,          ONLY : nr1sex, nr2sex, nr3sex, nrsex, nlsex, ecutsex, &
                               nr1sco, nr2sco, nr3sco, nrsco, nlsco, ecutsco, &
                               ngmsig, ngmsex, ecutsig, ngmsco, ngmgrn, ngmpol, nbnd_sig
  USE gvect,            ONLY : nl, g
  USE kinds,            ONLY : DP
  USE gwsymm,           ONLY : ngmunique, ig_unique, use_symm, sym_friend, sym_ig

  
  IMPLICIT NONE
  INTEGER :: iq, ik
  INTEGER :: ios
  LOGICAL :: do_band, do_iq, setup_pw, exst
  CHARACTER (LEN=9)   :: code = 'GW'
  CHARACTER (LEN=256) :: auxdyn
 
!Variables to split Planewaves perturbations across processors.
  INTEGER :: igstart, igstop, ngpool, ngr, igs, ig
  COMPLEX(DP), ALLOCATABLE :: scrcoul_g(:,:,:,:)

!Initialize MPI, clocks, print initial messages
!/Modules/mp.f90

#ifdef __PARA
   CALL mp_startup ( )
   IF (nimage>1) CALL io_path_start()
#endif
   
   !/Modules/environment.f90 prints out all parallel information and opening message.
    CALL environment_start ( code )

    WRITE(stdout, '(/5x, "Reading variables")') 
    CALL gwq_readin()
    WRITE(stdout, '(/5x, "Finished reading variables")') 

!   HL
!   Check stop init Modules/check_stop.f90
!   This module contains functions to check if the code should
!   be smoothly stopped.

    CALL check_stop_init()

!   This routine checks the initial status of the GW run, initializes the qmesh, and prepares
!   the control of the dispersion calculation. 
    CALL check_initial_status(auxdyn)

!   Generate frequency grid for GW convolution. 
    CALL freqbins()
    CALL clean_pw( .FALSE. )
    CALL allocate_fft()
    CALL ggen()

!   Generating Exchange and Correlation grid.
    CALL sig_fft_g(nr1sco, nr2sco, nr3sco, nrsco, ecutsig, 1)
    CALL sig_fft_g(nr1sex, nr2sex, nr3sex, nrsex, ecutsex, 2)
    CALL clean_pw( .FALSE. )

    ALLOCATE ( scrcoul_g( ngmpol, ngmpol, nfs, 1) )
    ALLOCATE ( ig_unique( ngmpol) )
    ALLOCATE ( sym_ig(ngmpol))
    ALLOCATE ( sym_friend(ngmpol))

    iuncoul = 28
    lrcoul = 2 * ngmpol * ngmpol * nfs

    iungreen = 31
    lrgrn  = 2 * ngmgrn * ngmgrn

if(multishift) then
    iunresid = 34
    lrresid  = 2*npwx

! Could probably keep the alphabeta coefficients in memory.
    iunalphabeta = 35
    lralphabeta  = 4
endif

IF (ionode) THEN
       iuncoul = 28
       lrcoul = 2 * ngmpol * ngmpol * nfs
       CALL diropn (iuncoul, 'coul', lrcoul, exst)
!   Green's function file
       iungreen = 31
       lrgrn  = 2 * ngmgrn * ngmgrn
       CALL diropn (iungreen, 'green', lrgrn, exst)
!   Sigma file
       iunsigma = 32
       lrsigma = 2 * ngmsco * ngmsco * nwsigma
       CALL diropn(iunsigma, 'sigma', lrsigma, exst)
!   Should sigma_ex need to be written to file:
       iunsex = 33
       lrsex = 2 * ngmsex * ngmsex
       CALL diropn(iunsex, 'sigma_ex', lrsex, exst)
!   Should sigma_extra need to be written to file:
       iunsigext = 36
       lrsigext = 2 * ngmsco * ngmsco
       CALL diropn(iunsigext, 'sig_ext', lrsigext, exst)
ENDIF

IF(do_coulomb) THEN
     DO iq = w_of_q_start, nqs
        scrcoul_g(:,:,:,:) = (0.0d0, 0.0d0)
!Prepare k, k+q grids, run nscf calculation, find small group of q.
        CALL prepare_q(do_band, do_iq, setup_pw, iq)
        CALL run_pwscf(do_band)
        CALL initialize_gw()


!Determine the unique G vectors in the small group of q if symmetry is being used.
!If not then all the vectors up to ngmpol (correlation cutoff) are treated as unique and calculated.
        if(use_symm) then
           WRITE(6,'("")')
           WRITE(6,'("SYMMETRIZING COULOMB Perturbations")')
           WRITE(6,'("")')
           CALL stern_symm()
        else
           ngmunique = ngmpol
           do ig = 1, ngmpol
              ig_unique(ig) = ig
           enddo
        endif
!Distribute unique G-vectors between processors:
#ifdef __PARA
      npool = nproc / nproc_pool
      if (npool.gt.1) then
      ! number of g-vec per pool and reminder
        ngpool = ngmunique / npool
        ngr = ngmunique - ngpool * npool
      ! the remainder goes to the first ngr pools
        if ( my_pool_id < ngr ) ngpool = ngpool + 1
        igs = ngpool * my_pool_id + 1
        if ( my_pool_id >= ngr ) igs = igs + ngr
      ! the index of the first and the last g vec in this pool
        igstart = igs
        igstop = igs - 1 + ngpool
        write (stdout,'(/4x,"Max n. of PW perturbations per pool = ",i5)') igstop-igstart+1
      else
#endif
       igstart = 1
       igstop = ngmunique
#ifdef __PARA
      endif
#endif
!CALCULATE W(G,G';iw):
       if((igstop-igstart+1).ne.0) then
           CALL coulomb(iq, igstart, igstop, scrcoul_g)
       endif

!COLLECT G-VECTORS:
       CALL mp_barrier(inter_pool_comm)
       CALL mp_sum ( scrcoul_g, inter_pool_comm )
       CALL mp_barrier(inter_pool_comm)

!Write W_{q}(G,G';iw) to file:
       IF (ionode) THEN
           CALL unfold_w(scrcoul_g,iq)
!If solve direct now need to invert epsilon:
           IF(solve_direct.and.tinvert) CALL invert_epsilon(scrcoul_g, iq)
           CALL davcio(scrcoul_g, lrcoul, iuncoul, iq, +1, ios)
       ENDIF

       CALL mp_barrier(inter_pool_comm)
       CALL clean_pw_gw(iq)
       CALL mp_barrier(inter_pool_comm)

       if(do_q0_only) GOTO 124
   END DO
   WRITE(stdout, '("Finished Calculating Screened Coulomb")') 
ENDIF
124 CONTINUE
!Free up memory from scrcoul_g 
   call mp_barrier(inter_pool_comm)
   DEALLOCATE( scrcoul_g )
   DEALLOCATE( ig_unique )

   DO ik = 1, 1
   !DO ik = 1, num_k_pts
       if(do_green.and.multishift) CALL diropn(iunresid, 'resid', lrresid, exst)
       if(do_green.and.multishift) CALL diropn(iunalphabeta, 'alphbet', lralphabeta, exst)
       xq(:) = xk_kpoints(:, ik)
       WRITE(6,'(4x,"Sigma_k", 3f12.7)') xk_kpoints(:,ik)
       do_iq=.TRUE.
       lgamma = ( xq(1) == 0.D0 .AND. xq(2) == 0.D0 .AND. xq(3) == 0.D0 )
       setup_pw = .TRUE.
       do_band  = .TRUE.
! Generates small group of k and then forms IBZ_{k}.
       CALL run_pwscf_green(do_band)
       CALL initialize_gw()
! CALCULATE G(r,r'; w) 
! WRITE(stdout, '(/5x, "GREEN LINEAR SYSTEM SOLVER")')
       if(do_green) write(6,'("Do green_linsys")')

       if(do_green.and.(.not.multishift)) then
            CALL green_linsys(ik)
            call mp_barrier()
            if(do_sigma_c) CALL sigma_c(ik)
       endif

! CALCULATE Sigma_corr(r,r';w) = i\int G(r,r'; w + w')(W(r,r';w') - v(r,r')) dw'
       if(do_green.and.multishift) then 
            CALL green_linsys_shift(ik)
            call mp_barrier()
            if(do_sigma_c) CALL sigma_c(ik)
       endif

       if(do_green.and.multishift) then 
          CLOSE(UNIT = iunresid, STATUS = 'DELETE')
          CLOSE(UNIT = iunalphabeta, STATUS = 'DELETE')
       endif

       if(ionode) then 
!CALCULATE Sigma_ex(r,r') = iG(r,r')v(r,r')
         if(do_sigma_exx)   CALL sigma_exch(ik)
       endif

       if(ionode) WRITE(6, '("Finished CALCULATING SIGMA")') 
       CALL mp_barrier(inter_pool_comm)
       CALL clean_pw_gw(ik)
       CALL mp_barrier(inter_pool_comm)
   ENDDO

   DO ik = 1, 1
         if(do_sigma_matel)  then
        !Calculates QP Corrections for bands 1:nbnd_sig.
            nbnd = nbnd_sig 
            xq(:) = xk_kpoints(:, ik)
            do_iq=.TRUE.
            lgamma = ( xq(1) == 0.D0 .AND. xq(2) == 0.D0 .AND. xq(3) == 0.D0 )
            setup_pw = .TRUE.
            do_band  = .TRUE.
            CALL run_pwscf_green(do_band)
            CALL initialize_gw()
            WRITE(6,'(4x,"Sigma exchange")')
            !if(ionode.and.do_sigma_exx) CALL sigma_exchg(ik)
            CALL mp_barrier(inter_pool_comm)
            CALL sigma_matel(ik) 
            CALL mp_barrier(inter_pool_comm)
            CALL clean_pw_gw(ik)
         endif
   ENDDO

   call mp_barrier(inter_pool_comm)
   CALL gw_writefile('init',0)
   CALL clean_input_variables()
   CALL collect_grid_files()
   CALL destroy_status_run()
   IF (bands_computed) CALL print_clock_pw()
   CALL stop_gw( .TRUE. )
END PROGRAM gw
