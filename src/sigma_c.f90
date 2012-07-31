SUBROUTINE sigma_c(ik0) 
! G TIMES W PRODUCT
  USE io_global,     ONLY : stdout, ionode_id, ionode
  USE kinds,         ONLY : DP
  USE lsda_mod,      ONLY : nspin
  USE constants,     ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE disp,          ONLY : nqs, nq1, nq2, nq3, wq, x_q, xk_kpoints
  USE control_gw,    ONLY : lgamma, eta, godbyneeds, padecont, cohsex
  USE klist,         ONLY : wk, xk
  USE io_files,      ONLY : prefix, iunigk, prefix, tmp_dir
  USE wvfct,         ONLY : nbnd, npw, npwx, igk, g2kin, et
  USE eqv,           ONLY : evq, eprec
  USE freq_gw,       ONLY : fpol, fiu, nfs, nfsmax,&
                            nwcoul, nwgreen, nwalloc, nwsigma, wtmp, wcoul,&
                            wgreen, wsigma, wsigmamin, wsigmamax,&
                            deltaw, wcoulmax, ind_w0mw, ind_w0pw
  USE units_gw,      ONLY : iuncoul, iungreen, iunsigma, lrsigma, lrcoul, lrgrn, iuwfc, lrwfc
  USE qpoint,        ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE gwsigma,       ONLY : ngmsco, sigma, sigma_g, nrsco, nlsco, fft6_g2r, ecutsco, ngmsig,&
                            nr1sco, nr2sco, nr3sco, ngmgrn, ngmpol
  USE gvect,         ONLY : g, ngm, ecutwfc, nl
  USE cell_base,     ONLY : tpiba2, tpiba, omega, alat
  USE symm_base,     ONLY : nsym, s, time_reversal, t_rev, ftau, invs
  USE mp_global,     ONLY : inter_pool_comm, intra_pool_comm, mp_global_end, mpime, npool, &
                            nproc_pool, me_pool, my_pool_id, nproc
  USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
  USE fft_scalar,    ONLY : cfft3d

  IMPLICIT NONE

  COMPLEX(DP)         :: ci, czero
  COMPLEX(DP)         :: aux (nrsco)
!For running PWSCF need some variables 
  LOGICAL             :: pade_catch
!Pade arrays
  COMPLEX(DP), ALLOCATABLE :: z(:), u(:), a(:)
!W arrays 
  COMPLEX(DP), ALLOCATABLE :: scrcoul_g (:,:,:)
  COMPLEX(DP), ALLOCATABLE :: scrcoul_g_R (:,:,:)
  COMPLEX(DP), ALLOCATABLE :: scrcoul_pade_g (:,:)
  COMPLEX(DP), ALLOCATABLE :: scrcoul(:,:)
!G arrays:
  COMPLEX(DP), ALLOCATABLE :: greenf_g(:,:), greenfp(:,:), greenfm(:,:)
!Integration Variable 
  COMPLEX(DP) :: cprefac
!FREQUENCY GRIDS/COUNTERS
  INTEGER  :: iwim, iw, ikq 
  INTEGER  :: iw0, iw0mw, iw0pw
  REAL(DP) :: w_ryd(nwcoul)
!COUNTERS
  INTEGER :: ig, igp, irr, icounter, ir, irp
  INTEGER :: iqstart, iqstop, iqs, nkr
  INTEGER :: iq, ipol, iqrec
  INTEGER :: ikmq, ik0, ik, nkpool
  INTEGER :: rec0, ios
  INTEGER :: counter, ierr
!SYMMETRY
  REAL(DP)              :: xq_ibk(3), xq_ibz(3)
  INTEGER               :: isym
  INTEGER, ALLOCATABLE  :: gmapsym(:,:)
  COMPLEX(DP)           :: eigv(ngm,48)

!CHECK FOR NAN's
  REAL(DP)     :: ar, ai

!File related:
  character(len=256) :: tempfile, filename
!Complete file name
  integer*8 :: unf_recl

#define DIRECT_IO_FACTOR 8 
! #ifdef __PARA
!       scrcoul_g = czero
!       if (me.eq.1.and.mypool.eq.1) then
! #endif
! #ifdef __PARA
!       endif
!       use poolreduce to broadcast the results to every pool
!       call poolreduce ( 2 * ngms * ngms * nwim, scrcoul_g)
! #endif
! iG(W-v)

!I think these should all (necessarily) stay as ngmsig
   ALLOCATE ( scrcoul_g      (ngmpol, ngmpol, nfs)     )
   ALLOCATE ( scrcoul_g_R    (ngmpol, ngmpol, nfs)     )
   ALLOCATE ( scrcoul_pade_g (ngmpol, ngmpol)          )
   ALLOCATE ( greenf_g       (ngmgrn, ngmgrn)          )

!These go on the big grid...
   ALLOCATE ( scrcoul        (nrsco, nrsco)            )
   ALLOCATE ( greenfp        (nrsco, nrsco)            )
   ALLOCATE ( greenfm        (nrsco, nrsco)            )
   ALLOCATE ( sigma          (nrsco, nrsco, nwsigma)   )
   ALLOCATE ( gmapsym        (ngm, 48)                 )
!  Pade Approximants.
   ALLOCATE  (z(nfs), a(nfs))

! Array for coulomb frequencies.
   w_ryd(:) = wcoul(:)/RYTOEV

   if(allocated(sigma)) then
     WRITE(6,'(4x,"Sigma allocated")')
   else
     WRITE(6,'(4x,"Sigma too large!")')
     CALL mp_global_end()
     STOP
   endif


   WRITE(6," ")
   WRITE(6,'(4x,"Direct product GW for k0(",i3," ) = (",3f12.7," )")') ik0, (xk(ipol, ik0), ipol=1,3)
   WRITE(6," ")
   WRITE(6,'(4x, "ngmsco, ", i4, " nwsigma, ", i4)') ngmsco, nwsigma
   WRITE(6,'(4x, "nrsco, ", i4, " nfs, ", i4)') nrsco, nfs

   ci = (0.0d0, 1.d0)
   czero = (0.0d0, 0.0d0)
   sigma(:,:,:) = (0.0d0, 0.0d0)

   CALL start_clock('sigmac')
!Here we generate the G-vector rotation array.
   CALL gmap_sym(nsym, s, ftau, gmapsym, eigv, invs)

!Set appropriate weights for points in the brillouin zone.
!Weights of all the k-points are in odd positions in list.
!nksq is number of k points not including k+q.
   do iq = 1, nksq
      if(lgamma) then
         write(6, '(" lgamma ")')
         wq(iq) = 0.5d0*wk(iq) 
      else
         wq(iq) = 0.5d0*wk(2*iq-1) 
      endif
   enddo

!Every processor needs access to the files: _gw0si.coul1 and _gw0si.green1
call mp_barrier(inter_pool_comm)

#ifdef __PARA
if(.not.ionode) then
!OPEN coulomb file (only written to by head node).
   filename = trim(prefix)//"."//"coul1"
   tempfile = trim(tmp_dir) // trim(filename)
   unf_recl = DIRECT_IO_FACTOR * int(lrcoul, kind=kind(unf_recl))
   open(iuncoul, file = trim(adjustl(tempfile)), iostat = ios, &
   form = 'unformatted', status = 'OLD', access = 'direct', recl = unf_recl)

!OPEN green file (only written to by head node as well).
   filename = trim(prefix)//"."//"green1"
   tempfile = trim(tmp_dir) // trim(filename)
   unf_recl = DIRECT_IO_FACTOR * int(lrgrn, kind=kind(unf_recl))
   open(iungreen, file = trim(adjustl(tempfile)), iostat = ios, &
   form = 'unformatted', status = 'OLD', access = 'direct', recl = unf_recl)
endif
#endif
#ifdef __PARA
      npool = nproc / nproc_pool
      write(stdout,'("npool", i4, i5)') npool, nksq
      if (npool.gt.1) then
      ! number of g-vec per pool and reminder
        nkpool = nksq / npool
        nkr = nksq - nkpool * npool

      ! the remainder goes to the first nkr pools
        if ( my_pool_id < nkr ) nkpool = nkpool + 1
        iqs = nkpool * my_pool_id + 1
        if ( my_pool_id >= nkr ) iqs = iqs + nkr

      ! the index of the first and the last g vec in this pool
        iqstart = iqs
        iqstop  = iqs - 1 + nkpool
        write (stdout,'(/4x,"Max n. of Kpoint in Sigma_C per pool = ",i5)') iqstop-iqstart+1
      else
#endif
       iqstart = 1
       iqstop = nksq
#ifdef __PARA
      endif
#endif
!ONLY PROCESSORS WITH K points to process: 
IF(iqstop-iqstart+1.ne.0) THEN
!  WRITE(1000+mpime, '("mpime ", i4, "  iqstart, iqstop: ", 2i5)')mpime, iqstart, iqstop
   DO iq = iqstart, iqstop
      if (lgamma) then
          ikq = iq
      else
!k+q is in even positions of list (k,k+q)
          ikq = 2*iq
      endif
!q point for convolution \sum_{q \in IBZ_{k}} G_{k+q} W_{-q}
!-q = k0 - (k0 + q)
      xq_ibk(:) = xk_kpoints(:,ik0) - xk(:,ikq)

!Find which symmetry operation rotates xq_ibk back to
!The irreducible brillouin zone and which q \in IBZ it corresponds to.
     call find_q_ibz(xq_ibk, s, iqrec, isym)
     write(6, *)  
     write(6, '("xq_IBK point")')
     write(6, '(3f11.7)') xq_ibk
     write(6, '("equivalent xq_IBZ point, symop, iqrec")')
     write(6, '(3f11.7, 3i4)') x_q(:,iqrec), isym, iqrec, iuncoul
     write(6,*)

!lrcoul = 2 * ngmsig * ngmsig * nfs.
!CALL davcio(scrcoul_g, lrcoul, iuncoul, iq, -1 )
!Read the equivalent q_point in the IBZ of crystal:

!If we get really facsistic about memory could:
!ALLOCATE(scrcoul_g)
     CALL davcio(scrcoul_g, lrcoul, iuncoul, iqrec, -1)

!Rotate G_vectors for FFT.
!In EPW FG checked that gmapsym(gmapsym(ig,isym),invs(isym)) = ig
!I have checked that here as well and it works.
!write(6,'(11i4)')gmapsym(gmapsym(:ngmsco,isym),invs(isym))
!Take a look at rotated G-vect map:
!for systems with ftau != 0 should also multiply each matrix element
!by phase, e^{iG'\tau}.
!Another one of these nested loops. 
!Two strategies or alternatives:
!1) Could pad scrcoul_g_R so that all the vectors still fall in side of it up to ngmsco
!   then trim them off later.
!2) Modify gmapsym so that it only keeps vectors up to ngmsco...

     do ig = 1, ngmpol
        do igp = 1, ngmpol
           if((gmapsym(ig,isym).lt.ngmpol).and.(gmapsym(igp,isym).lt.ngmpol)) then
               do iwim = 1, nfs
                  !Is the symmetry stuff killing me now?
                  scrcoul_g_R(gmapsym(ig,isym), gmapsym(igp,isym), iwim) = scrcoul_g(ig,igp,iwim)
               enddo
           endif
        enddo
     enddo


!deallocate(scrcoul_g)
!Start integration over iw +/- wcoul. 
    WRITE(6,'("Starting Frequency Integration")')
    DO iw = 1, nwcoul
        scrcoul_pade_g(:,:) = (0.0d0, 0.0d0)
        do ig = 1, ngmpol
           do igp = 1, ngmpol
             do iwim = 1, nfs
                 z(iwim) = dcmplx( 0.d0, fiu(iwim))
             !normal ordering.
             !a(iwim) = scrcoul_g (ig,igp,iwim)
             !For rotated g-vects.
                 a(iwim) = scrcoul_g_R (ig,igp,iwim)
             enddo
             pade_catch=.false.
             do iwim = 1, nfs
                 ar = real(a(iwim))
                 ai = aimag(a(iwim))
                 if ( ( ar .ne. ar ) .or. ( ai .ne. ai ) ) then
                      a(:) = (0.0d0, 0.0d0)
                      pade_catch = .true.
                 endif
             enddo
             if(padecont) then
                call pade_eval ( nfs, z, a, dcmplx( w_ryd(iw), eta), scrcoul_pade_g (ig,igp))
             else if (cohsex) then
                  scrcoul_pade_g(ig,igp) = a(1)
             else if(godbyneeds) then
                     scrcoul_pade_g(ig,igp)=(a(2)/(dcmplx(w_ryd(iw),0.0d0)-a(1)+ci*eta ))-((a(2))/(dcmplx(w_ryd(iw), 0.0d0) + a(1) - ci*eta))
             else 
                  WRITE(6,'("No screening model chosen!")')
                  STOP
                  call mp_global_end()
             endif
           enddo
        enddo
!W(G,G';w)
         czero = (0.0d0, 0.0d0)
         scrcoul(:,:) = czero
         do ig = 1, ngmpol
            aux(:) = czero
            do igp = 1, ngmpol
               aux(nlsco(igp)) = scrcoul_pade_g(ig,igp)
            enddo
            call cfft3d (aux, nr1sco, nr2sco, nr3sco, nr1sco, nr2sco, nr3sco, +1)
            do irp = 1, nrsco
               scrcoul(ig, irp) = aux(irp) / omega
            enddo
         enddo
         do irp = 1, nrsco
            aux = czero
            do ig = 1, ngmpol
               aux(nlsco(ig)) = conjg(scrcoul(ig,irp))
            enddo
            call cfft3d (aux, nr1sco, nr2sco, nr3sco, nr1sco, nr2sco, nr3sco, +1)
            scrcoul(1:nrsco,irp) = conjg ( aux )
         enddo
!Now have W(r,r';omega')
         cprefac = (deltaw/RYTOEV) * wq(iq) * (0.0d0, 1.0d0)/ tpi
         DO iw0 = 1, nwsigma
            iw0mw = ind_w0mw (iw0,iw)
            iw0pw = ind_w0pw (iw0,iw)
!rec0 = (iw0mw-1) * 1 * nqs + (ik0-1) * nqs + (iq-1) + 1
            rec0 = (iw0mw-1) * 1 * nksq + (iq-1) + 1
            CALL davcio( greenf_g, lrgrn, iungreen, rec0, -1 )
            greenfm(:,:) = czero
            do ig = 1, ngmgrn
               aux(:) = czero
               do igp = 1, ngmgrn
                  aux(nlsco(igp)) = greenf_g(ig,igp)
               enddo
               call cfft3d (aux, nr1sco, nr2sco, nr3sco, nr1sco, nr2sco, nr3sco, +1)
               do irp = 1, nrsco
                  greenfm(ig, irp) = aux(irp)/omega
               enddo
            enddo
            do irp = 1, nrsco
               aux = czero
               do ig = 1, ngmgrn
                  aux(nlsco(ig)) = conjg(greenfm(ig,irp))
               enddo
               call cfft3d (aux, nr1sco, nr2sco, nr3sco, nr1sco, nr2sco, nr3sco, +1)
               greenfm(1:nrsco,irp) = conjg ( aux )
            enddo
!Now have G(r,r';omega-omega')
!rec0 = (iw0pw-1) * 1 * nqs + (ik0-1) * nqs + (iq-1) + 1
            rec0 = (iw0pw-1) * 1 * nksq + (iq-1) + 1
            CALL davcio(greenf_g, lrgrn, iungreen, rec0, -1)
!Inlining FFT:
            greenfp(:,:) = czero
            do ig = 1, ngmgrn
               aux(:) = czero
               do igp = 1, ngmgrn
                  aux(nlsco(igp)) = greenf_g(ig,igp)
               enddo
              call cfft3d (aux, nr1sco, nr2sco, nr3sco, nr1sco, nr2sco, nr3sco, +1)
               do irp = 1, nrsco
                  greenfp(ig, irp) = aux(irp) / omega
               enddo
            enddo
            do irp = 1, nrsco
               aux = czero
               do ig = 1, ngmgrn
                  aux(nlsco(ig)) = conjg(greenfp(ig,irp))
               enddo
               call cfft3d (aux, nr1sco, nr2sco, nr3sco, nr1sco, nr2sco, nr3sco, +1)
               greenfp(1:nrsco,irp) = conjg ( aux )
            enddo
!Now have G(\r,\r';\omega+\omega')
               sigma (:,:,iw0) = sigma (:,:,iw0) + cprefac * (greenfp(:,:) + greenfm(:,:)) * scrcoul(:,:)
        ENDDO !on iw0  
      ENDDO ! on frequency convolution over w'

!now could add extra ananlytic piece: -\inf to wsigmamin-wcoul
!and from wsigmamin+wcoul to +\inf
!G(r,r';w'-w0) =  FFT[\delta(\G,G')/(w-w0-w1)]
!W(r,r';w')
!\Sigma(r, r'; w0) = \Sigma(r,r'; w0) + \int_{\-inf}^{wcoul) G(r,r'; w0 - w')W(r,\r';w')dw'
!                                     + \int_{wcoul}^{\inf) G(r,r'; w0 - w')W(r,\r';w')dw'
    ENDDO ! end loop iqstart, iqstop 
ENDIF

    DEALLOCATE ( gmapsym          )
    DEALLOCATE ( greenfm, greenfp )
    DEALLOCATE ( greenf_g         )
    DEALLOCATE ( scrcoul          )
    DEALLOCATE ( scrcoul_pade_g   )
    DEALLOCATE ( scrcoul_g        )
    DEALLOCATE ( z,a )

#ifdef __PARA
    CALL mp_barrier(inter_pool_comm)
#endif

#ifdef __PARA
!sum up sigma across processors.
    CALL mp_sum(sigma, inter_pool_comm)
#endif __PARA

CALL mp_barrier(inter_pool_comm)

IF (ionode) then
 !I use to say there was no chance of doing this
 !on the big grid, but now we do it on a bigger grid
 !and later using windowing we will do it on the full density grid.
  ALLOCATE ( sigma_g (ngmsco, ngmsco, nwsigma))
  if(allocated(sigma_g)) then
     WRITE(6,'(4x,"Sigma_g allocated")')
  else
     WRITE(6,'(4x,"Sigma_g too large!")')
     CALL mp_global_end()
     STOP
  endif

  WRITE(6,'(4x,"Sigma in G-Space")')
    sigma_g = (0.0d0,0.0d0)
    do iw = 1, nwsigma
      do ir = 1, nrsco
        aux = (0.0d0, 0.0d0)
        do irp = 1, nrsco
           aux(irp) = sigma(ir,irp,iw)
        enddo
        call cfft3d (aux, nr1sco, nr2sco, nr3sco, nr1sco, nr2sco, nr3sco, -1)
        do igp = 1, ngmsco
           sigma (ir, igp, iw) = aux(nlsco(igp))
        enddo
      enddo
      do igp = 1, ngmsco
        aux = czero
        do ir = 1, nrsco
          aux(ir) = conjg ( sigma(ir,igp,iw) )
        enddo
        call cfft3d (aux, nr1sco, nr2sco, nr3sco, nr1sco, nr2sco, nr3sco, -1)
        do ig = 1, ngmsco
           sigma (ig,igp,iw) = conjg ( aux( nlsco( ig )) ) * omega
        enddo
      enddo
!CALL davcio (sigma_g, lrsigma, iunsigma, iw, 1)
    enddo

    do ig = ngmsco + 1, nrsco
       do igp = ngmsco + 1, nrsco
          do iw = 1, nwsigma
             sigma (ig, igp, iw) = (0.0d0, 0.0d0)
          enddo 
       enddo
    enddo

!sigma_g = sigma(1:ngmsco,1:ngmsco,:)
    do ig = 1, ngmsco
     do igp = 1, ngmsco
        do iw = 1, nwsigma
           sigma_g(ig,igp,iw)  = sigma(ig,igp,iw)
        enddo 
     enddo
    enddo
!Now write Sigma in G space to file. 
    WRITE(6,'(4x,"Writing Sigma to File")')
!HL Original:
!Just storing in first record no matter what k-point.
    CALL davcio (sigma_g, lrsigma, iunsigma, 1, 1)

!or could store sigma same way green's fxn is stored...
    CALL stop_clock('sigmac') 
    DEALLOCATE ( sigma_g  )
ENDIF !ionode

    call mp_barrier(inter_pool_comm)
    DEALLOCATE ( sigma    )
    RETURN
END SUBROUTINE sigma_c

