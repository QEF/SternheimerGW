  !
  !----------------------------------------------------------------
  subroutine coulomb_q0G0 ( vr, xq0, nw, w, scrcoul )
  !----------------------------------------------------------------
  ! 
  ! the screened Coulomb interaction W_q (G,G',w) for a given q
  ! 
  ! we include the spin degeneracy in the kpoint weights
  !
  ! this subroutine is specifically for the case q=0, G=0
  ! We use a small q (called xq0) and calculate the KS states 
  ! on the grid {k}+xq0
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
  real(dbl) :: xq0(3), ui, uj, uk, qg2, qgg, xktmp(3), g0(3), alpha_mix
  real(dbl) :: qg, rcut, spal
  integer :: iq, count, i, j, k, ik, ipol, ikk, ikq, ig0, igmg0, nw
  ! igmg0 = i of (G - G0)
  !
  integer :: maxscf, maxbcgsolve
  integer :: shift(nks)
  integer, parameter :: ng0vec = 27
  complex(dbl) :: scrcoul(nrs, nrs, nw)
  integer :: ig, iw, igp, ierr, ibnd, ios, recl, unf_recl, ikf, fold(nq)
  real(dbl) :: g2kin(ngm), et(nbnd_occ, nks), w(nw), w_ryd(nw)
  integer :: nwim
  ! number of frequencies on the imaginary axis
  real(dbl), allocatable :: wim(:), wim_ryd(:)
  ! frequencies on the imaginary axis
  complex(dbl), allocatable :: z(:), u(:), a(:)
  ! these are for the Pade continuation to the real axis
  complex(dbl) :: vr(nr), psi(ngm, nbnd_occ), dvbare(nr), dvscf(nr)
  complex(dbl) :: dvscfp (nr), dvscfm (nr) 
  ! these are for +w and -w
  complex(kind=DP) :: aux (ngm, nbnd_occ), evq (ngm, nbnd_occ)
  real(dbl) :: kplusg(3)
  character (len=3) :: nd_nmbr0
  character (len=256) :: barfile, dwfpfile, dwfmfile
  !
  logical :: convt
  ! return .true. is convergence has been achieved in solve_linter_dyn
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
  !
!@
  ! this simply means that we perform the brute-force minimization
  ! by fully relaxing dpsi at every dvscf update
  maxscf = 1000
  maxbcgsolve = 5
  alpha_mix = 0.5
!@

  !
  ! read from file the frequencies on the imaginary axis (eV)
  ! (every node reads in parallel - it's a tiny file) 
  !
  open ( 45, file = "./imfreq.dat", form = 'formatted', status = 'unknown')
  read (45, *) 
  read (45, *) nwim
  if (.not.allocated(wim)) allocate ( wim(nwim) )
  if (.not.allocated(wim_ryd)) allocate ( wim_ryd(nwim) )
  if (.not.allocated(z)) allocate ( z(nwim) )
  if (.not.allocated(u)) allocate ( u(nwim) )
  if (.not.allocated(a)) allocate ( a(nwim) )
  do iw = 1, nwim
    read (45,*) wim(iw) 
    if (wim(iw).lt.0.d0) call error ('coulomb','imaginary frequencies must be positive',1)
  enddo
  if (nwim.gt.20) call error ('coulomb','too many imaginary frequencies',nwim)
  close(45)
  !
 
  ! initialize
  scrcoul = 0.d0
  w_ryd = w/ryd2ev
  wim_ryd = wim/ryd2ev

  !
  !  generate uniform {k} and {k+q} grids 
  !
  !  The {k} grid is taken to coincide with the {q} grid generated
  !  in gwhs.f90, the {k+q} grid is obtained as {k}+xq0
  !  where xq0 is the small q vector taken to represent the limit q->0
  !  In this case we need to calculate the occupied states for the 
  !  {k}+xq0 grid
  !
  do ik = 1, nksq
    xk (:, ik) = xq (:, ik) 
    ! include spin degeneracy
    wk = two / float ( nksq )
  enddo
  !
  !  double the grid and add k+q in the even positions
  !
  do ik = nksq, 1, -1
    ikk = 2 * ik - 1
    ikq = 2 * ik
    xk(:,ikk) = xk(:,ik)
    xk(:,ikq) = xk(:,ikk)+xq0
    wk(ikk) = wk(ik)
    wk(ikq) = 0.d0
  enddo
  !
  !  generate occupied wavefunctions for the k+xq0 points
  !
  write(6,'(/4x,a)') repeat('-',67)
  write(6,'(4x,"{k}+xq0 grid: occupied eigenvalues (eV)")')
  write(6,'(4x,a/)') repeat('-',67)
  !
  recl = 2 * nbnd_occ * ngm  ! 2 stands for complex
  !
  do ik = 1, nksq
    !
    ikk = 2 * ik - 1
    ikq = 2 * ik
    !
    ! the eigenvalues on the {k} grid
    ! 
    et(:,ikk) = eval_occ(:,ik)
    !
    ! the k-dependent kinetic energy in Ry
    ! [Eq. (14) of Ihm,Zunger,Cohen J Phys C 12, 4409 (1979)]
    !
    do ig = 1, ngm
      kplusg = xk(:,ikq) + g(:,ig)
      g2kin ( ig ) = tpiba2 * dot_product ( kplusg, kplusg )
    enddo
    !
    call eigenstates2 ( xk(:,ikq), vr, g2kin, evq, et(:,ikq) )
    !
    !  direct write to file - take into account the k/k+q alternation
    !
    write ( iunwfc, rec = ikq, iostat = ios) evq
    !
    write ( 6, '(4x,"k+xq0 (",i3," )",10(3x,f7.3))') ik, et(:,ikq)*ryd2ev
    !
  enddo
  write(6,'(4x,a/)') repeat('-',67)
  !
  !  only G=0 perturbation here 
  !
  ig = 1
  qg2 = (g(1,ig)+xq0(1))**2.d0 + (g(2,ig)+xq0(2))**2.d0 + (g(3,ig)+xq0(3))**2.d0
  !
! write(6,'(4x,"ig = ",i5)') ig
  !
  do iw = 1, nwim
!@  do iw = 1, 2
    !
!   write(6,'(4x,3x,"iw = ",i5)') iw
    write(6,'(4x,"Screened Coulomb: q =",3f7.3,"  G =",3f7.3,"  w(eV) =",3f7.3)') &
      xq0,g(:,ig), wim(iw)
    !
    dvbare = czero
    dvscf  = czero
    !
    dvbare ( nl ( ig ) ) = cone 
    call cfft3 ( dvbare, nr1, nr2, nr3,  1)
    !
!      !
!      ! DEBUG - tested - OK when I use 555-111. NOT OK when I use 555-000
!      !
!      ! here I only perform 1 iteration. In output we have
!      ! eps instead of eps^-1 (eps = delta - v_c * chi_0 )
!      !
!      call solve_linter_dyn_nonSCF ( dvbare, dvscf, xq0, et, vr, w_ryd(iw) )
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
    call solve_linter_dyn ( dvbare, dvscf, xq0, et, vr, wim_ryd(iw), maxscf, maxbcgsolve, alpha_mix, .true., convt )
    !
    ! transform dvscf to G space 
    !
    call cfft3 ( dvscf , nr1, nr2, nr3, -1)
    !
    write(6,'(4x,"INVEPS: ",3f9.5)') wim(iw), dvscf ( nl(1) )
    !                                         ^^^^^^^^^^^^^^^
    !                                          eps^-1(0,0,q)
#ifdef __PARA
!   write(1000+mypool,'(4x,2f9.5)') dvscf ( nl(1) )
#endif
    !
    ! symmetrized inverse dielectric matrix (finite limits for q->0)
    !
    do igp = 1, ngms
      dvscf ( nl(igp) ) = dvscf ( nl(igp) ) * &
         sqrt ( (g(1,igp)+xq0(1))**2.d0 + (g(2,igp)+xq0(2))**2.d0 + (g(3,igp)+xq0(3))**2.d0 ) / &
         sqrt ( (g(1,ig )+xq0(1))**2.d0 + (g(2,ig )+xq0(2))**2.d0 + (g(3,ig )+xq0(3))**2.d0 )
!     write(6,'("INVEPS",2x,2i4,2f20.10,6(2x,f9.1))') ig, igp,    &
!         real(dvscf ( nl(igp) )), aimag(dvscf ( nl(igp) )),      &
!         g(1,ig), g(2,ig), g(3,ig), g(1,igp), g(2,igp), g(3,igp)
!     write(500+iw,'(2i4,2f20.10)') ig, igp, real(dvscf ( nl(igp) )), aimag(dvscf ( nl(igp) ))
    enddo
    ! 
    ! keep only the G-vectors 1:ngms for the screened Coulomb 
    !
    do igp = 1, ngms
       qgg = sqrt( (g(1,ig )+xq0(1))**2.d0 + (g(2,ig )+xq0(2))**2.d0 + (g(3,ig )+xq0(3))**2.d0 ) &
           * sqrt( (g(1,igp)+xq0(1))**2.d0 + (g(2,igp)+xq0(2))**2.d0 + (g(3,igp)+xq0(3))**2.d0 )
       scrcoul (ig,igp,iw) = dvscf ( nl(igp) ) * dcmplx ( e2 * fpi / (tpiba2*qgg), zero )
     scrcoul (ig,igp,iw) = dvscf ( nl(igp) ) 
    enddo
    !
    ! Spencer/Alavi truncation of the bare coulomb interaction
    ! [PRB 77,193110 (2008]
    !
    ! WARNING: I tested epsilon, but never tested W
    !
    rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))
    qg = sqrt( (g(1,ig )+xq0(1))**2.d0 + (g(2,ig )+xq0(2))**2.d0 + (g(3,ig )+xq0(3))**2.d0 )
    spal = one - cos ( rcut * tpiba * qg )
    !
    do igp = 1, ngms
      scrcoul (ig,igp,iw) = scrcoul (ig,igp,iw) * dcmplx ( spal, zero)
    enddo
    !
  enddo
  !
  ! Here goes the Pade continuation to real axis
  !
  do igp = 1, ngms
    !
    ! Pade input points on the imaginary axis
    !
    do iw = 1, nwim
      z(iw) = dcmplx( 0.d0, wim_ryd(iw))
      u(iw) = scrcoul (ig,igp,iw)
    enddo
    !
    ! Pade coefficients
    !
    call pade_coeff ( nwim, z, u, a)
    !
    ! Pade output points on the real axis (at a distance eta)
    !
    do iw = 1, nw
      call pade_eval ( nwim, z, a, dcmplx( w_ryd(iw), eta), scrcoul (ig,igp,iw))
    enddo
    !
  enddo
  !
! do iw = 1, nw
!   write(6,'(4x,"PADE INVEPS: ",3f9.5)') w(iw), scrcoul (1,1,iw) 
! enddo
  !
  close ( iubar, status = 'delete' ) 
  close ( iudwfp, status = 'delete' ) 
  close ( iudwfm, status = 'delete' ) 
  !
  return
  end subroutine coulomb_q0G0
  !----------------------------------------------------------------
  !
