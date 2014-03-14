SUBROUTINE sigma_c_serial(ik0, ikq, green, sigma, iw0)
! G TIMES W PRODUCT
  USE io_global,     ONLY : stdout, ionode_id, ionode
  USE kinds,         ONLY : DP
  USE lsda_mod,      ONLY : nspin
  USE constants,     ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE disp,          ONLY : nqs, nq1, nq2, nq3, wq, x_q, xk_kpoints
  USE control_gw,    ONLY : lgamma, eta, godbyneeds, padecont, cohsex, modielec
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
  USE gwsigma,       ONLY : ngmsco, nrsco, nlsco, fft6_g2r, ecutsco, ngmsig,&
                            nr1sco, nr2sco, nr3sco, ngmgrn, ngmpol
  USE gvect,         ONLY : g, ngm, ecutwfc, nl
  USE cell_base,     ONLY : tpiba2, tpiba, omega, alat, at
  USE symm_base,     ONLY : nsym, s, time_reversal, t_rev, ftau, invs, nrot
  USE mp_global,     ONLY : inter_pool_comm, intra_pool_comm, mp_global_end, mpime, npool, &
                            nproc_pool, me_pool, my_pool_id, nproc
  USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
  USE fft_scalar,    ONLY : cfft3d
  USE modes,         ONLY : nsymq, invsymq, gi, gimq, irgq, irotmq, minus_q
  USE wavefunctions_module, ONLY : evc
  USE control_flags,        ONLY : noinv

  IMPLICIT NONE
  COMPLEX(DP)         :: ci, czero
  COMPLEX(DP)         :: phase
  COMPLEX(DP)         :: aux (nrsco)
!For running PWSCF need some variables 
  LOGICAL             :: pade_catch
  LOGICAL             :: found_q
  LOGICAL             :: limit, limq, inv_q, found
!Pade arrays
  COMPLEX(DP), ALLOCATABLE :: z(:), u(:), a(:)
!W arrays 
  COMPLEX(DP), ALLOCATABLE :: scrcoul_g (:,:,:)
  COMPLEX(DP), ALLOCATABLE :: scrcoul_g_R (:,:,:)
  COMPLEX(DP), ALLOCATABLE :: scrcoul_pade_g (:,:)
  COMPLEX(DP), ALLOCATABLE :: scrcoul(:,:)
  COMPLEX(DP)              :: sigma(nrsco, nrsco)
!v array
!COMPLEX(DP), ALLOCATABLE ::  barcoul(:,:), barcoulr(:,:), barcoul_R(:,:)
  REAL(DP) :: qg2, qg, qxy, qz
!G arrays:
  COMPLEX(DP), ALLOCATABLE :: greenf_g(:,:), greenfr(:,:)
  COMPLEX(DP), INTENT(IN) :: green(ngmgrn, ngmgrn, nwgreen)
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
  INTEGER :: inversym, screening
!SYMMETRY
  REAL(DP)              :: xq_ibk(3), xq_ibz(3)
  INTEGER               :: isym, jsym
  INTEGER, ALLOCATABLE  :: gmapsym(:,:)
!For G^NA
  INTEGER     :: igkq_ig(npwx) 
  INTEGER     :: igkq_tmp(npwx) 
  INTEGER     :: ss(3,3)
  COMPLEX(DP), ALLOCATABLE  ::  eigv(:,:)
!q-vector of coulomb potential xq_coul := k_{0} - xk(ik)
  REAL(DP) :: xq_coul(3)
  REAL(DP) :: rcut, spal, zcut
  INTEGER  :: ibnd
!CHECK FOR NAN's
  REAL(DP)     :: ar, ai
!For dirac delta fxn.
  REAL(DP)     :: dirac, x, support

!File related:
  character(len=256) :: tempfile, filename
!Complete file name
  integer*8 :: unf_recl

#define DIRECT_IO_FACTOR 8 
!iG(W-v)
   ALLOCATE ( scrcoul_g      (ngmpol, ngmpol, nfs)     )
   ALLOCATE ( scrcoul_g_R    (ngmpol, ngmpol, nfs)     )
   ALLOCATE ( scrcoul_pade_g (ngmpol, ngmpol)          )
   ALLOCATE ( greenf_g       (ngmgrn, ngmgrn)          )


!These go on the big grid...
   ALLOCATE ( scrcoul        (nrsco, nrsco)            )
   ALLOCATE ( greenfr        (nrsco, nrsco)            )

!Technically only need gmapsym up to ngmpol or ngmgrn...
   ALLOCATE ( gmapsym  (ngm, nrot)    )
   ALLOCATE ( eigv     (ngm, nrot)    )
   ALLOCATE ( z(nfs), a(nfs), u(nfs)  )

   w_ryd(:) = wcoul(:)/RYTOEV

!   WRITE(6," ")
!   WRITE(6,'(4x,"Direct product GW for k0(",i3," ) = (",3f12.7," )")') ik0, (xk(ipol, ik0), ipol=1,3)
!   WRITE(6,'(4x, "ngmsco, ", i4, " nwsigma, ", i4)') ngmsco, nwsigma
!   WRITE(6,'(4x, "nrsco, ", i4, " nfs, ", i4)') nrsco, nfs

   ci    = (0.0d0, 1.d0)
   czero = (0.0d0, 0.0d0)
   limit = .false.

   greenf_g(:,:) = czero

   CALL gmap_sym(nrot, s, ftau, gmapsym, eigv, invs)

!  WRITE(6,'("nsym, nsymq, nsymbrav ", 3i4)'), nsym, nsymq, nrot 
!Set appropriate weights for points in the brillouin zone.
!Weights of all the k-points are in odd positions in list.
!nksq is number of k points not including k+q.
   wq(:) = 0.0d0
   do iq = 1, nksq
      if(lgamma) then
         wq(iq) = 0.5d0*wk(iq) 
      else
         wq(iq) = 0.5d0*wk(2*iq-1) 
      endif
   enddo

!ONLY PROCESSORS WITH K points to process: 
!WRITE(1000+mpime, '("mpime ", i4, "  iqstart, iqstop: ", 2i5)')mpime, iqstart, iqstop
!if (nksq.gt.1) rewind (unit = iunigk)
!   if (lgamma) then
!       ikq = iq
!   else
!k+q is in even positions of list (k,k+q)
!       ikq = 2*iq
!   endif
!q point for convolution \sum_{q \in IBZ_{k}} G_{k+q} W_{-q}
!  q = (k0 + q) - k0
!  xq_ibk(:) = xk(:,ikq) - xk_kpoints(:, ik0)
!  by actually searching for -q and calculating w of -q  
!  we ensure there is no confusion in which symmetry routine is required.
!  SYMMFIX
!  -q =  k0 - (k0 + q)

   xq_ibk(:) = xk_kpoints(:, ik0) - xk(:, ikq)

!Find which symmetry operation rotates xq_ibk back to
!The irreducible brillouin zone and which q \in IBZ it corresponds to.
!q is stored in the list x_q as positive q but all the calculations have
!been done at -q therefore we are just going to calculate \Sum G_{k+q}W_{-q}
   inv_q=.false.
   call find_q_ibz(xq_ibk, s, iqrec, isym, found_q, inv_q)

   if(lgamma) npwq=npw 

!   write(6, *)  
!   write(6, '("xq_IBK point")')
!   write(6, '(3f11.7)') xq_ibk
!   write(6, '("equivalent xq_IBZ point, symop, iqrec")')
!   write(6, '(3f11.7, 2i4)') x_q(:,iqrec), isym, iqrec
!   write(6,*)

   if(iw0.eq.1) then
       write(1000+mpime, *)  
       write(1000+mpime, '("xq_IBK point")')
       write(1000+mpime, '(3f11.7)') xq_ibk
       write(1000+mpime, '("equivalent xq_IBZ point, symop, iqrec")')
       write(1000+mpime, '(3f11.7, 2i4)') x_q(:, iqrec), isym, iqrec
    endif

!Dielectric Function should be written to file at this point
!So we read that in, rotate it, and then apply the coulomb operator.
     scrcoul_g(:,:,:)   = dcmplx(0.0d0, 0.0d0)
     scrcoul_g_R(:,:,:) = dcmplx(0.0d0, 0.0d0)
     if(.not.modielec) CALL davcio(scrcoul_g, lrcoul, iuncoul, iqrec, -1)

!Rotate G_vectors for FFT.
     do igp = 1, ngmpol
        do ig = 1, ngmpol
           if((gmapsym(ig,isym).le.ngmpol).and.(gmapsym(igp,isym).le.ngmpol) &
               .and.(gmapsym(ig,isym).gt.0).and.(gmapsym(ig,isym).gt.0)) then
               do iwim = 1, nfs
!Rotating dielectric matrix with phase factor:
               phase = eigv(ig,isym)*conjg(eigv(igp,isym))
               scrcoul_g_R(ig, igp, iwim) = scrcoul_g(gmapsym(ig,isym), gmapsym(igp,isym),iwim)*phase
               enddo
           endif
        enddo
     enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Generate bare coulomb:                      !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!sax cutoff:                                      !!
!!!rcut = 0.50d0*minval(sqrt(sum(at**2,1)))*alat*tpi!!
!!!rcut = rcut-rcut/50.0d0                          !!
!!!                                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))
if(.not.modielec) then
DO iw = 1, nfs
   DO ig = 1, ngmpol
!SPHERICAL SCREENING
       qg2 = (g(1,ig) + xq_ibk(1))**2 + (g(2,ig) + xq_ibk(2))**2 + (g(3,ig)+xq_ibk(3))**2
       limq = (qg2.lt.eps8) 
       IF(.not.limq) then
           do igp = 1, ngmpol
              scrcoul_g_R(ig, igp, iw) = scrcoul_g_R(ig,igp,iw)*dcmplx(e2*fpi/(tpiba2*qg2), 0.0d0)
           enddo
       ENDIF
       qg = sqrt(qg2)
       spal = 1.0d0 - cos(rcut*sqrt(tpiba2)*qg)
!Normal case using truncated coulomb potential.
       if(.not.limq) then
          do igp = 1, ngmpol
              scrcoul_g_R(ig, igp, iw) = scrcoul_g_R(ig,igp,iw)*dcmplx(spal, 0.0d0)
          enddo
       else
!should only occur case iq->0, ig = 0 use vcut (q(0) = (4pi*e2*Rcut^{2})/2
!             write(6,'("Taking Limit.")')
!             write(6,*) (rcut)
!             write(6,*) (fpi*e2*(rcut**2))/2.0d0
!             write(6,*) ig, iw
!             write(6,*) g(:, ig)
!             write(6,*) xq_ibk(:)
!             write(6,*) qg2
!for omega=0,q-->0, G=0 the real part of the head of the dielectric matrix should be real
!we enforce that here:
             if(iw.eq.1) then
                scrcoul_g_R(ig, igp, iw) = real(scrcoul_g_R(ig,igp,iw))
             endif
             do igp = 1, ngmpol
                scrcoul_g_R(ig, igp, iw) = scrcoul_g_R(ig,igp,iw)*dcmplx((fpi*e2*(rcut**2))/2.0d0, 0.0d0)
             enddo
       endif
   ENDDO!ig
ENDDO!nfs
endif

!zeroing wings of W again!
if(ikq.eq.1) then
!    WRITE(6,'(4x,"Zeroing wings")')
    do iw = 1, nfs
     do ig = 2, ngmpol
        scrcoul_g_R(ig,1,iw)  = dcmplx(0.0d0, 0.d0)
     enddo
     do igp = 2, ngmpol
          scrcoul_g_R(1,igp,iw) = dcmplx(0.0d0, 0.d0)
     enddo
    enddo
endif

    if(.not.modielec) then 
        if(godbyneeds) then
        do ig = 1, ngmpol
          do igp = 1, ngmpol 
!         for godby-needs plasmon pole the algebra is done assuming real frequency*i...
!         that is: the calculation is done at i*wp but we pass a real number as the freq. 
               do iw = 1, nfs
                 z(iw) = dcmplx(aimag(fiu(iw)), 0.0d0)
                 u(iw) = scrcoul_g_R(ig, igp, iw)
               enddo
               CALL godby_needs_coeffs(nfs, z, u, a)
               do iw = 1, nfs 
!         just overwrite scrcoul_g_R with godby-needs coefficients.
                  scrcoul_g_R (ig, igp, iw) = a(iw)
               enddo
          enddo
        enddo
        else if (padecont) then
          do igp = 1, ngmpol
           do ig = 1, ngmpol
!     Pade input points on the imaginary axis
              do iw = 1, nfs
                 z(iw) = fiu(iw)
                 u(iw) = scrcoul_g_R (ig, igp, iw)
              enddo
              call pade_coeff ( nfs, z, u, a)
!     Overwrite scrcoul with Pade coefficients to be passed to pade_eval.
             do iw = 1, nfs 
                scrcoul_g_R (ig, igp, iw) = a(iw)
             enddo
           enddo !enddo on ig
        enddo  !enddo on igp
        else if(.not.padecont.and..not.godbyneeds) then
                 WRITE(6,'("No screening model chosen!")')
        endif
    endif
!Start integration over iw +/- wcoul. 
    do iw = 1, nwcoul
        scrcoul_pade_g(:,:) = (0.0d0, 0.0d0)
        if(.not.modielec) then
          do ig = 1, ngmpol
               do igp = 1, ngmpol
                  do iwim = 1, nfs
                      z(iwim) = fiu(iwim)
                      a(iwim) = scrcoul_g_R (ig,igp,iwim)
                enddo
               pade_catch=.false.
                if(padecont) then
                  call pade_eval ( nfs, z, a, scrcoul_g_R(ig,igp,1), dcmplx(w_ryd(iw), eta), scrcoul_pade_g (ig,igp))
               else if(godbyneeds) then
                 scrcoul_pade_g(ig,igp) = a(2)/(dcmplx(w_ryd(iw)**2,0.0d0)-(a(1)-(0.0d0,1.0d0)*eta)**2)
               else 
                    WRITE(6,'("No screening model chosen!")')
                    STOP
                    call mp_global_end()
               endif
            enddo
          enddo
        else if(modielec)  then
           do ig = 1, ngmpol
              call mod_diel(ig, xq_ibk, w_ryd(iw), scrcoul_pade_g(ig,ig), 1)
              qg2 = (g(1,ig) + xq_ibk(1))**2 + (g(2,ig) + xq_ibk(2))**2 + (g(3,ig)+xq_ibk(3))**2
              limq = (qg2.lt.eps8) 
              IF(.not.limq) then
                 scrcoul_pade_g(ig, ig) = scrcoul_pade_g(ig,ig)*dcmplx(e2*fpi/(tpiba2*qg2), 0.0d0)
              ENDIF
              qg = sqrt(qg2)
              spal = 1.0d0 - cos(rcut*sqrt(tpiba2)*qg)
!Normal case using truncated coulomb potential.
              if(.not.limq) then
                    scrcoul_pade_g(ig, ig) = scrcoul_pade_g(ig,ig)*dcmplx(spal, 0.0d0)
              else
                    scrcoul_pade_g(ig, ig) = scrcoul_pade_g(ig,ig)*dcmplx((fpi*e2*(rcut**2))/2.0d0, 0.0d0)
              endif
           enddo
        endif

!W_{q=0}(G,G';w)
!if(ikq.eq.1) then
!    WRITE(6,'(4x,"Zeroing wings")')
!     do ig = 2, ngmpol
!        scrcoul_pade_g(ig,1)  = dcmplx(0.0d0, 0.d0)
!     enddo
!     do igp = 2, ngmpol
!          scrcoul_pade_g(1,igp) = dcmplx(0.0d0, 0.d0)
!     enddo
!endif

        czero = (0.0d0, 0.0d0)
        scrcoul(:,:) = czero
        call fft6(scrcoul_pade_g(1,1), scrcoul(1,1), 1)

!Now have W(r,r';omega')
!simpson quadrature: int_w1^wN f(w)dw = deltaw*
![ 1/3 ( f1 + fN ) + 4/3 sum_even f_even + 2/3 sum_odd f_odd]
        if(lgamma) then
            iq = ikq
        else
            iq = ikq/2
            print*, iq, wq(iq)
        endif

        cprefac = (deltaw/RYTOEV) * wq(iq) * (0.0d0, 1.0d0)/ tpi

        if ( iw/2*2.eq.iw ) then
            cprefac = cprefac * 4.d0/3.d0
        else
            cprefac = cprefac * 2.d0/3.d0
        endif

!Presumably the fxn should be zero at the cutoff but might as well be consistent...
        if (iw.eq.1) then
             cprefac = (1.0d0/3.0d0)*(deltaw/RYTOEV) * wq(iq) * (0.0d0, 1.0d0)/ tpi
        else if (iw.eq.nwcoul) then
             cprefac = (1.0d0/3.0d0)*(deltaw/RYTOEV) * wq(iq) * (0.0d0, 1.0d0)/ tpi
        endif

        iw0mw = ind_w0mw (iw0,iw)
        iw0pw = ind_w0pw (iw0,iw)


        rec0 = (iw0mw-1) + 1
        greenf_g(:,:) = green(:,:, rec0)
        greenfr(:,:) = czero
        call fft6(greenf_g(1,1), greenfr(1,1), 1)

        sigma (:,:) = sigma (:,:) + cprefac * greenfr(:,:)*scrcoul(:,:)

!Now have G(r,r';omega-omega')
        rec0 = (iw0pw-1) + 1
        greenf_g(:,:) = green(:,:,rec0)
        greenfr(:,:) = czero
        call fft6(greenf_g(1,1), greenfr(1,1),1)

        sigma (:,:) = sigma (:,:) + cprefac * greenfr(:,:)*scrcoul(:,:)

    enddo

    DEALLOCATE ( gmapsym          )
    DEALLOCATE ( greenfr          )
    DEALLOCATE ( greenf_g         )
    DEALLOCATE ( scrcoul          )
    DEALLOCATE ( scrcoul_pade_g   )
    DEALLOCATE ( scrcoul_g, scrcoul_g_R )
    DEALLOCATE ( z,a,u )
    RETURN
END SUBROUTINE sigma_c_serial
