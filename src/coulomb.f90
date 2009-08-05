  !
  !----------------------------------------------------------------
  subroutine coulomb ( vr, xxq, nw, w, scrcoul, igstart, igstop)
  !----------------------------------------------------------------
  ! 
  ! the screened Coulomb interaction W_q (G,G',w) for a given q
  ! 
  ! we include the spin degeneracy in the kpoint weights
  !
  !----------------------------------------------------------------
  !
  use parameters
  use constants
  use gspace
  use kspace
#ifdef __PARA
  USE para
  USE mp_global,  ONLY : nproc, mpime, nproc_pool, my_pool_id, me_pool
  USE mp, ONLY:  mp_barrier
#endif
  implicit none
  !
  real(dbl) :: xxq(3), ui, uj, uk, qg2, qgg, xktmp(3), g0(3), alpha_mix
  integer :: iq, count, i, j, k, ik, ipol, ikk, ikq, ig0, igmg0, nw
  ! igmg0 = i of (G - G0)
  !
  integer ::  maxscf, maxbcgsolve, maxscf_
  integer :: ngpool, igs, ngr, igstart, igstop
  integer :: shift(nks)
  integer, parameter :: ng0vec = 27
  complex(dbl) :: scrcoul(nrs, nrs, nw)
  integer :: ig, iw, igp, ierr, ibnd, ios, recl, unf_recl, ikf, fold(nq)
  real(dbl) :: g2kin(ngm), et(nbnd_occ, nks), w(nw), w_ryd(nw)
  complex(dbl) :: vr(nr), psi(ngm, nbnd_occ), dvbare(nr), dvscf(nr)
  complex(dbl) :: dvscfp (nr), dvscfm (nr)
  ! these are for +w and -w
  complex(dbl) :: z(nw), u(nw), a(nw)
  ! these are for the Pade continuation to the real axis
  complex(kind=DP) :: aux (ngm, nbnd_occ), evq (ngm, nbnd_occ)
  character (len=3) :: nd_nmbr0
  character (len=256) :: barfile, dwfpfile, dwfmfile
  !
  logical :: convt
  ! return .true. from solve_linter if convergence has been achieved

  !
  recl = 2 * nbnd_occ * ngm  ! 2 stands for complex
  unf_recl = DIRECT_IO_FACTOR * recl
  barfile  = './silicon'//'.bar'
  dwfpfile = './silicon'//'.dwfp'
  dwfmfile = './silicon'//'.dwfm'
  !
#ifdef __PARA
  call set_ndnmbr ( mypool, me_pool, nprocp, npool, nd_nmbr0)
  barfile =  trim(barfile)//'.'//nd_nmbr0
  dwfpfile = trim(dwfpfile)//'.'//nd_nmbr0
  dwfmfile = trim(dwfmfile)//'.'//nd_nmbr0
#endif
  !
  open ( iubar, file = barfile, iostat = ios, form = 'unformatted', &
       status = 'unknown', access = 'direct', recl = unf_recl)
  open ( iudwfp, file = dwfpfile, iostat = ios, form = 'unformatted', &
       status = 'unknown', access = 'direct', recl = unf_recl)
  open ( iudwfm, file = dwfmfile, iostat = ios, form = 'unformatted', &
       status = 'unknown', access = 'direct', recl = unf_recl)
!@
  !
  ! IMPORTANT: if we decrease maxbcgsolve down to 1 or 2, the
  ! final result can become wrong (as compared to maxbcgsolve = Inf)
  !
  maxscf = 1000
  maxbcgsolve = 5
  alpha_mix = 0.5
!@

  ! initialize
  scrcoul = 0.d0
  w_ryd = w/ryd2ev

  !
  !  generate uniform {k} and {k+q} grids 
  !
  !  The {k} grid is taken to coincide with the {q} grid generated
  !  in gwhs.f90, the {k+q} grid is obtained by folding
  !  In this way we calculate the occupied states only once in gwhs.f90
  !
  do ik = 1, nksq
    xk (:, ik) = xq (:, ik) 
    ! include spin degeneracy
    wk = two / float ( nksq )
  enddo
  !
  !  find folding index and G-vector map
  !
  do ik = 1, nksq
    call ktokpmq (xk(:,ik), xxq, +1, fold(ik) )
    !
    !  the folding G-vector
    !
    g0 = xk(:, fold(ik)) - ( xk(:,ik) + xxq )
    !
    shift(ik) = 0
    do ig0 = 1, ng0vec
      if ( ( abs(g0vec(1,ig0)-g0(1)).lt.1.d-8 ) .and. &
           ( abs(g0vec(2,ig0)-g0(2)).lt.1.d-8 ) .and. &
           ( abs(g0vec(3,ig0)-g0(3)).lt.1.d-8 ) ) then
         shift(ik) = ig0
      endif
    enddo
    if (shift(ik).eq.0) call error ('coulomb','folding vector not found',0)
    !
  enddo
  !
  !  double the grid and add k+q in the even positions
  !
  do ik = nksq, 1, -1
    ikk = 2 * ik - 1
    ikq = 2 * ik
    xk(:,ikk) = xk(:,ik)
    xk(:,ikq) = xk(:,ik)+xxq
    wk(ikk) = wk(ik)
    wk(ikq) = 0.d0
  enddo
  !
  do ik = 1, nksq
    !
    ikk = 2 * ik - 1
    ikq = 2 * ik
    !
    !  the folded k+q
    !
    ikf = fold(ik)
    !
    !  read occupied wavefunctions for the folded k+q point
    !  c_{k+q} (G) = c_{k+q+G0} (G-G0)
    !
    read ( iunwfc, rec = 2 * ikf - 1, iostat = ios) aux
    !
    !  set the phase factor of evq for the folding
    !
    !  WARNING: here we loose some information since
    !  the sphere G-G0 has half the surface outside the
    !  cutoff. It is important to make sure that the cutoff
    !  used is enough to guarantee that the lost half-surface
    !  does not create problems. In the EPW code I solved
    !  the issue by translating the G-sphere in the fft
    !  mapping. It's painful, so I will do it only in extremis. 
    !  FG Aug 2008 
    !
    do ibnd = 1, nbnd_occ
      do ig = 1, ngm
        igmg0 = gmap (ig, shift(ik))
        if (igmg0.eq.0) then 
           evq (ig,ibnd) = czero
        else
           evq (ig,ibnd) = aux ( igmg0, ibnd)
        endif
      enddo
    enddo
    !
    !  and write it again in the right place
    !
    write ( iunwfc, rec = ikq, iostat = ios) evq
    !
    ! DEBUG: here I checked that by running
    ! call eigenstates ( xk(:,ikq), vr, g2kin, evq, eval_occ(:,ikf) ) 
    ! the evq (wfs at k+q) are the same as those obained above
    ! (within a phase factor and gauge in degenerate cases -
    ! I actually checked sum_ibnd |evq(:,ibnd)|^2 )
    !
    ! the eigenvalues
    ! 
    et(:,ikk) = eval_occ(:,ik)
    et(:,ikq) = eval_occ(:,ikf) 
    !
  enddo
  !
  !  loop on bare perturbations ig and fixed q point
  !
  !  [note: no barriers inside this loop as the number of perturbation
  !  is different across processors]
  !
  do ig = igstart, igstop 
    !
!   write(6,'(4x,"ig = ",i5)') ig
    qg2 = (g(1,ig)+xxq(1))**2.d0 + (g(2,ig)+xxq(2))**2.d0 + (g(3,ig)+xxq(3))**2.d0
    !
!@    do iw = 1, nw
   do iw = 1, 2
      !
      dvbare = czero
      dvscf  = czero
      if (qg2 > 1.d-8) then
        !
        write(6,'(4x,"Screened Coulomb: q =",3f7.3,"  G =",3f7.3,"  w(eV) =",3f7.3)') &
          xxq,g(:,ig), w(iw)
        !
        write(6,'(4x,3x,"iw = ",i5)') iw
        !
        dvbare ( nl ( ig ) ) = cone ! this is in SIZE-order of G-vectors 
        !                           ! (look at fft_test.f90)
        call cfft3 ( dvbare, nr1, nr2, nr3,  1)

!      !
!      ! DEBUG - tested - OK when I use 555-111. NOT OK when I use 555-000
!      !
!      ! here I only perform 1 iteration. In output we have
!      ! eps instead of eps^-1 (eps = delta - v_c * chi_0 )
!      !
!      call solve_linter_dyn_nonSCF ( dvbare, dvscf, xxq, et, vr, w_ryd(iw) )
!      !
!      ! go to G-space, this is  -v_c * chi_0
!      call cfft3 ( dvscf , nr1, nr2, nr3, -1)
!      ! now add the delta function
!      dvscf ( nl (ig) ) = dvscf ( nl (ig) ) + cone
!      ! transform dvbare to G space
!      call cfft3 ( dvbare, nr1, nr2, nr3, -1)
!      write(6,'(f9.5,5x,2f9.5,5x,2f9.5)') w(iw), dvbare ( nl(1) ), dvscf ( nl(1) )
!      !                                                            ^^^^^^^^^^^^^^^
!      !                                                            eps(0,0,xq0,w)
!      ! END DEBUG
!      !

        !
        ! solve self-consistently the linear system for this perturbation 
        ! _dyn is for the dynamical case (w/=0)
        !
        call solve_linter_dyn ( dvbare, dvscf, xxq, et, vr, w_ryd(iw), maxscf, maxbcgsolve, alpha_mix, .true., convt )
        !
        if (.not.convt) then
          !
          ! force full CG minimization at each SCF step
          !
          write(800,'(4x,"Screened Coulomb: q =",3f7.3,"  G'' =",3f7.3,"  w(eV) =",3f7.3)') &
             xxq,g(:,ig), w(iw)
          write(800,'(4x,"Convergence of mixed minimization not achieved: forcing full SCF minimization")')
          maxscf_ = 10000
          call solve_linter_dyn ( dvbare, dvscf, xxq, et, vr, w_ryd(iw), maxscf_, maxbcgsolve, alpha_mix, .true., convt )
          if (.not.convt) call error ('solve_linter_dyn','scf convergence not achieved',1)
        endif
        !
        ! transform dvscf to G space 
        !
        call cfft3 ( dvscf , nr1, nr2, nr3, -1)
        !
        write(6,'(4x,2f9.5)') dvscf ( nl(1) )
        !                     ^^^^^^^^^^^^^^^
        !                      eps^-1(0,0,q)
#ifdef __PARA
!       write(1000+mypool,'(4x,2f9.5)') dvscf ( nl(1) )
#endif
        !
      endif
      !
      ! symmetrized inverse dielectric matrix (finite limits for q->0)
      !
      do igp = 1, ngms
        dvscf ( nl(igp) ) = dvscf ( nl(igp) ) * &
           sqrt ( (g(1,igp)+xxq(1))**2.d0 + (g(2,igp)+xxq(2))**2.d0 + (g(3,igp)+xxq(3))**2.d0 ) / &
           sqrt ( (g(1,ig )+xxq(1))**2.d0 + (g(2,ig )+xxq(2))**2.d0 + (g(3,ig )+xxq(3))**2.d0 )
!       write(6,'("INVEPS",2x,2i4,2f20.10,6(2x,f9.1))') ig, igp,    &
!           real(dvscf ( nl(igp) )), aimag(dvscf ( nl(igp) )),      &
!           g(1,ig), g(2,ig), g(3,ig), g(1,igp), g(2,igp), g(3,igp)
!       write(500+iw,'(2i4,2f20.10)') ig, igp, real(dvscf ( nl(igp) )), aimag(dvscf ( nl(igp) ))
      enddo
      !
      ! keep only the G-vectors 1:ngms for the screened Coulomb 
      !
      do igp = 1, ngms
        qgg = sqrt( (g(1,ig )+xxq(1))**2.d0 + (g(2,ig )+xxq(2))**2.d0 + (g(3,ig )+xxq(3))**2.d0 ) &
            * sqrt( (g(1,igp)+xxq(1))**2.d0 + (g(2,igp)+xxq(2))**2.d0 + (g(3,igp)+xxq(3))**2.d0 )
        scrcoul (ig,igp,iw) = dvscf ( nl(igp) ) * dcmplx ( e2 * fpi / (tpiba2*qgg), zero )
      enddo
      !
      ! scrcoul (-,igp,-) is now in SIZE order of G-vectors
      !
    enddo
    !
    ! Here goes the Pade continuation to real axis
    !
    do igp = 1, ngms
      !
      ! Pade input points on the imaginary axis
      !
      do iw = 1, nw
        z(iw) = dcmplx( 0.d0, w_ryd(iw))
        u(iw) = scrcoul (ig,igp,iw)
      enddo 
      ! 
      ! Pade coefficients
      !
      call pade_coeff ( nw, z, u, a)
      !
      ! Pade output points on the real axis (at a distance eta)
      ! (I use the same grid for simplicity - this can be changed)
      !
      do iw = 1, nw
        call pade_eval ( nw, z, a, dcmplx( w_ryd(iw), eta), scrcoul (ig,igp,iw))
      enddo
      !
    enddo
    !
  enddo
  ! 
  close ( iubar, status = 'delete' ) 
  close ( iudwfp, status = 'delete' ) 
  close ( iudwfm, status = 'delete' ) 
  !
  return
  end subroutine coulomb
  !----------------------------------------------------------------
  !
