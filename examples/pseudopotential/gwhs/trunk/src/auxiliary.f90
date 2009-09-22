  !
  !-----------------------------------------------------------------------
  subroutine h_psi ( psi, hpsi, g2kin, vr)
  !-----------------------------------------------------------------------
  !
  ! hpsi = H_k | psi >
  !
  use gspace
  use constants
  use parameters
  implicit none
  !
  real(kind=DP) :: psi (ngm), hpsi (ngm), g2kin(ngm)
  complex(kind=DP) :: vr (nr), psic (nr), rhoca(nr1,nr2,nr3)
  integer :: ig, ir, i, j, k
  !
  !
  ! FFT test
  !
  allocate ( rhoca (nr1, nr2, nr3) )
  !
  do i = 1, nr1
   do j = 1, nr2
    do k = 1, nr3
      ir = i + (j-1) * nr1 + (k-1) * nr1 * nr2
      rhoca (i, j, k) = float(i) 
      vr ( ir ) = rhoca (i, j, k) 
      write (16,'(3i5,f9.5)') i,j,k,real(rhoca (i, j, k))
    enddo
   end do
  end do
  !
  ! go to G space
  call cfft3 ( vr, nr1, nr2, nr3, -1)
  psic = czero
  do ig = 1, ngm
    psic ( ig ) = vr ( nl ( ig ) )  ! psic is in order of increasing |G|^2
  enddo
  !
  ! now back to R space
  do ig = 1, ngm
    vr ( nl ( ig ) ) = psic (ig)
  enddo
  call cfft3 ( vr, nr1, nr2, nr3,  1)
  !
  do i = 1, nr1
   do j = 1, nr2
    do k = 1, nr3
      ir = i + (j-1) * nr1 + (k-1) * nr1 * nr2
      rhoca (i, j, k) = vr (ir)
      write (6,'(3i5,f9.5)') i,j,k,real(rhoca (i, j, k))
    enddo
   end do
  end do
  !
  stop
  !
  return
  end subroutine h_psi
  !
  !
  !----------------------------------------------------------------
  program epm
  !----------------------------------------------------------------
  ! 
  ! test-case for the GW-DFPT proposal 
  ! empirical pseudopotential method implementation
  ! for nonpolar tetrahedral semiconductors 
  !
  ! Fri May 30, 4.26 PM: 
  !  great! it gives exactly the same Si bandstructure as in
  !  Cohen&Bergstresser!
  !
  ! NOTE: I use only one G-grid. The reason is that, while V(G-G')
  ! has in principle components 2G_max, in the EPS scheme the  
  ! largest value happens for G2=11, therefore we can simply use
  ! the smooth grid...
  ! To perform a true scaling test we need to include the double
  ! grid (and count the operations on the smooth and on the dense
  ! grids).
  !
  ! NOTE2: the k-dependent Hamiltonian is real and symmetric 
  !  in reciprocal space
  !  (but not in real space because of the momentum component)
  !  therefore we can find eigenstates with c(G) all real
  !  (the eigenvectors of a rs matrix form a orthogonal matrix)
  !  of course psi(r) does not need to be real
  !
  ! Note that the Hamiltonian in G-space is real for the
  ! diamond structure
  !
  !----------------------------------------------------------------
  !
  use gspace
  use parameters
  use constants
  implicit none
  !
  ! variables
  !
  real(dbl) :: gcutm, arg
  real(dbl) :: xk (3,nk), deltag(3), kplusg(3), xq(3,nq)
  integer :: ig, ik, ig1, ig2, ideltag, notfound, ierr, ibnd, iq
  logical :: allowed, equiv
  real(dbl), allocatable :: t(:,:), ss(:), vs(:), h(:,:,:)
  complex(dbl), allocatable :: vr(:)
  real(kind=DP), parameter :: eps8 = 1.0D-8
  real(kind=DP), allocatable :: eval(:), u(:,:), fv1(:), fv2(:), hk(:,:)
  real(kind=DP), allocatable :: psi(:,:), e(:)
  !
  ! temporary variables - debug
  integer :: i, j, k
  real(kind=DP), allocatable :: psincr(:), hpsincr(:)
  complex(kind=DP), allocatable :: psincc(:)
  !

  !
  !----------------------------------------------------------------
  ! DEFINE THE CRYSTAL LATTICE
  !----------------------------------------------------------------
  !
  !  direct lattice vectors (cart. coord. in units of a_0)
  !  [Yu and Cardona, pag 23]
  !
  at( :, 1) = (/ 0.0,  0.5,  0.5 /)  ! a1
  at( :, 2) = (/ 0.5,  0.0,  0.5 /)  ! a2
  at( :, 3) = (/ 0.5,  0.5,  0.0 /)  ! a3
  !     
  !  reciprocal lattice vectors (cart. coord. in units 2 pi/a_0)
  !
  bg( :, 1) = (/ -1.0,  1.0,  1.0 /)  ! b1
  bg( :, 2) = (/  1.0, -1.0,  1.0 /)  ! b2
  bg( :, 3) = (/  1.0,  1.0, -1.0 /)  ! b3
  !
  !  atomic coordinates (cart. coord. in units of a_0)
  !
  tau( :, 1) = (/  0.125,  0.125,  0.125 /)  
  tau( :, 2) = (/ -0.125, -0.125, -0.125 /)  
  !
  !----------------------------------------------------------------
  ! GENERATE THE G-VECTORS
  !----------------------------------------------------------------
  !
  ! the G^2 cutoff in units of 2pi/a_0
  ! Note that in Ry units the kinetic energy is G^2, not G^2/2
  ! (note for the Hamiltonian we need to double the size, 2Gmax, hence the factor 4)
  !
  gcutm = four * ecutwfc / tpiba2      
! gcutm = ecutwfc / tpiba2      
  !
  ! set the fft grid
  !
  ! estimate nr1 and check if it is an allowed value for FFT
  !
  nr1 = 1 + int (2 * sqrt (gcutm) * sqrt( at(1,1)**2 + at(2,1)**2 + at(3,1)**2 ) ) 
  nr2 = 1 + int (2 * sqrt (gcutm) * sqrt( at(1,2)**2 + at(2,2)**2 + at(3,2)**2 ) ) 
  nr3 = 1 + int (2 * sqrt (gcutm) * sqrt( at(1,3)**2 + at(2,3)**2 + at(3,3)**2 ) ) 
  !
  do while (.not.allowed(nr1)) 
    nr1 = nr1 + 1
  enddo
  do while (.not.allowed(nr2)) 
    nr2 = nr2 + 1
  enddo
  do while (.not.allowed(nr3)) 
    nr3 = nr3 + 1
  enddo
  !
  call ggen( gcutm)
!@@
goto 10   ! use this to skip Haydok and do Coulomb only
!@@
  !
  !----------------------------------------------------------------
  ! GENERATE THE k-POINTS
  !----------------------------------------------------------------
  !
  ! brutally defined by hand - here I should read from file
  !
  do ik = 1, nk
    read (5,*) xk(1,ik), xk(2,ik), xk(3,ik) 
  enddo
  !

  !----------------------------------------------------------------
  ! CONSTRUCT THE HAMILTONIAN
  !----------------------------------------------------------------
  !
  allocate ( t(ngm, nk), ss(ngm), vs(ngm) )
  !
  ! the structure factor
  !
  do ig = 1, ngm
    arg = twopi * ( g(1,ig) * tau( 1, 1) + g(2,ig) * tau( 2, 1) + g(3,ig) * tau( 3, 1) )
    ss (ig) = cos ( arg )
  enddo
  !
  ! the empirical pseudopotential
  !
  vs = zero
  ! integer comparison - careful with other structures
  do ig = 1, ngm
    if     ( int ( gl(igtongl(ig)) ) .eq.  3 ) then
      vs (ig) =  v3
    elseif ( int ( gl(igtongl(ig)) ) .eq.  8 ) then
      vs (ig) =  v8
    elseif ( int ( gl(igtongl(ig)) ) .eq. 11 ) then
      vs (ig) =  v11
    endif
  enddo
  !
  ! the k-dependent kinetic energy in Ry
  ! [Eq. (14) of Ihm,Zunger,Cohen J Phys C 12, 4409 (1979)]
  !
  do ik = 1, nk
    do ig = 1, ngm
      kplusg = xk(:,ik) + g(:,ig)
      t ( ig, ik) = tpiba2 * dot_product ( kplusg, kplusg )
    enddo
  enddo
  !
  allocate ( h (ngm, ngm, nk) )
  do ig = 1, ngm
    do ik = 1, nk
        h ( ig, ig, ik) = t ( ig, ik)
    enddo
  enddo
  !
  do ig1 = 1, ngm
   do ig2 = ig1+1, ngm
     !
     ! define ideltag
     ! 
     deltag = g(:,ig1) - g(:,ig2)
     ideltag = 1 
     do while ( .not. equiv ( deltag, g(:,ideltag) ) .and. ideltag.lt.ngm )
       ideltag = ideltag + 1
     enddo   
     !
     ! the free-electron dispersions look ok (vs=0)
     ! when compared to Yu/Cardona Fig. 2.8
     ! 
     do ik = 1, nk
       if (ideltag.ne.ngm) &
          h ( ig1, ig2, ik) = ss (ideltag) * vs (ideltag) 
          h ( ig2, ig1, ik) = h ( ig1, ig2, ik) 
     enddo
     !
   enddo
  enddo
  !
  ! the empirical pseudopotential in real space 
  ! for further use in h_psi
  !
  allocate ( vr(nr) )
  vr = czero
  do ig = 1, ngm
    vr ( nl ( ig ) ) = dcmplx ( ss (ig) * vs (ig), zero )
  enddo
  call cfft3 ( vr, nr1, nr2, nr3,  1)
  !
  deallocate ( ss, vs)
  !
! !----------------------------------------------------------------
! ! BRUTE-FORCE DIAGONALIZATION
! !----------------------------------------------------------------
! !
! allocate ( eval(ngm), hk(ngm,ngm), u(ngm,ngm), fv1(ngm), fv2(ngm) )
! !
! do ik = 1, nk
!   !
!   hk = h(:,:,ik)
!   !
!   ! 0 = no eigenvectors
!   call rs ( ngm, ngm, hk, eval, 0, u, fv1, fv2, ierr)
!   !
!   write(10,'(i3,10(3x,f13.5))') ik, eval(1:nbnd)*ryd2ev
!   !
! enddo
! !
! deallocate ( eval, hk, u, fv1, fv2 )
  !
  !----------------------------------------------------------------
  ! CONJUGATE GRADIENTS AND GREEN'S FUNCTION FROM RECURSION METHOD
  !----------------------------------------------------------------
  !
  ! I checked that the band-structure from CG is identical to the
  ! one obtained by a full diagonalization -> FFT and H|psi> are ok.
  !
  allocate ( e (nbnd), psi (ngm, nbnd) )
  !
  write(6,'(/4x,a)') repeat('-',67)
  write(6,'(4x,"Eigenvalues (eV) from preconditioned conjugate gradient")') 
  write(6,'(4x,a/)') repeat('-',67)
  !
  allocate ( eval(ngm0), hk(ngm0,ngm0), u(ngm0,ngm0), fv1(ngm0), fv2(ngm0) )
  !
! do ik = 1, nk  !@@@ only for debug                       
  do ik = 1, 2 ! also tried 6,7 - fine
    !
    !----------------------------------------------------------------
    ! CONJUGATE GRADIENTS DIAGONALIZATION
    !----------------------------------------------------------------
    !
    ! trial wavefunctions for CG algorithm, cutoff ecut0 (modules.f90)
    !
    hk = h(1:ngm0,1:ngm0,ik)
    call rs ( ngm0, ngm0, hk, eval, 1, u, fv1, fv2, ierr)
    psi = zero
    do ibnd = 1, nbnd
     do ig = 1, ngm0
        psi(ig,ibnd) = u(ig,ibnd)
     enddo
    enddo
    !
    call cgdiag (nbnd, psi, e, t(1,ik), vr )
    !
    write(6,'(4x,i3,8(2x,f6.3))') ik, e(1:8)*ryd2ev
    write(6,'(4x,a/)') repeat('-',67)
    !
    !----------------------------------------------------------------
    ! GREEN'S FUNCTION G(r,r',k,w) IN THE PSINC BASIS 
    !----------------------------------------------------------------
    !
    ! I put this in the CG loop since I want to compare the GF
    ! calculated by expanding on empty states with the one obtained
    ! by the recursion method
    !
!@    do ig1 = 1, ngm
!@      do ig2 = 1, ngm
    do ig1 = 1, 5
      do ig2 = 30, 31
       !
!!     call green ( ig1, ig2, t(1,ik), vr )
       call green ( ig1, ig2, t(1,ik), vr, psi, e ) ! @@@ debug
       !
     enddo
    enddo
    !
  enddo
!@@
!stop ! keep this stop if only Green's function is needed
10 continue
!@@
  !
  !----------------------------------------------------------------
  ! SCREENED COULOMB INTERACTION
  !----------------------------------------------------------------
  !
  !  q point for test purposes - really need a loop here
  !
  xq(:, 1) = (/  0.001, 0.0, 0.0/)
  xq(:, 2) = (/  0.010, 0.0, 0.0/)
  xq(:, 3) = (/  0.050, 0.0, 0.0/)
  xq(:, 4) = (/  0.100, 0.0, 0.0/)
  xq(:, 5) = (/  0.150, 0.0, 0.0/)
  xq(:, 6) = (/  0.200, 0.0, 0.0/) 
  xq(:, 7) = (/  0.300, 0.0, 0.0/) 
  xq(:, 8) = (/  0.500, 0.0, 0.0/) 
  xq(:, 9) = (/  0.750, 0.0, 0.0/) 
  xq(:,10) = (/  1.000, 0.0, 0.0/) 
  xq(:,11) = (/  1.500, 0.0, 0.0/) 
  xq(:,12) = (/  2.000, 0.0, 0.0/) 
  xq(:,13) = (/  6.000, 0.0, 0.0/) 
  xq(:,14) = (/ 10.000, 0.0, 0.0/) 
  xq(:,15) = (/ 2.0, 0.0, 0.0/) 
  xq(:,16) = (/ 1.0, 1.0, 1.0/) 
  !
! do iq = 1, nq
! do iq = 1, 14
  do iq = 2, 2 
     call coulomb ( xq(:,iq) )
!    call dielec_mat ( xq(:,iq), nwcoul )
  enddo
  !
  !----------------------------------------------------------------
  ! SELF-ENERGY
  !----------------------------------------------------------------
  !
  end
  !





  !
  !----------------------------------------------------------------
  subroutine dielec_mat ( xxq, nw )
  !----------------------------------------------------------------
  ! 
  ! the static dielectric function eps(G=0,G'=0,q,w=0)
  ! 
  !----------------------------------------------------------------
  !
  use parameters
  use constants
  use gspace
  use kspace
  implicit none
  !
  integer :: iq, count, i, j, k, ipol, ikk, ikq, nw
  integer :: ig, iw, igp, ierr, ibnd, ios, recl, unf_recl, notconv
  integer :: ideltag, ig1, ig2, ig3, ik, jbnd
  real(dbl) :: xxq(3), ui, uj, uk, qg2, xxk(3), et(nbnd, nks), epsil, vqg, w
  real(dbl) :: g2kin(ngm), arg, deltag(3), avg_iter, precondition (ngm)
  real(DP), allocatable :: enval(:), u(:,:), fv1(:), fv2(:), hk(:,:), ss(:), vs(:)
  complex(dbl) :: vr(nr), psi(ngm, nbnd), dvbare(nr), dvscf(nr), ZDOTC, chi0
  complex(kind=DP) :: evc (ngm, nbnd), evq (ngm, nbnd), a(nksq,nbnd_occ,nbnd)
  logical :: equiv
  !
  allocate ( enval(ngm0), u(ngm0,ngm0), fv1(ngm0), fv2(ngm0), & 
    ss(ngm), vs(ngm), hk (ngm0, ngm0) )
  !
  allocate ( xk (3,nks), wk(nks) )
  !

  recl = 2 * nbnd * ngm  ! 2 stands for complex
  unf_recl = DIRECT_IO_FACTOR * recl
  open ( iunwfc, file = "./silicon.wfc", iostat = ios, form = 'unformatted', &
       status = 'unknown', access = 'direct', recl = unf_recl)
  !
  ! the empirical pseudopotential
  !
  vs = zero
  ! integer comparison - careful with other structures
  do ig = 1, ngm
    if     ( int ( gl(igtongl(ig)) ) .eq.  3 ) then
      vs (ig) =  v3
    elseif ( int ( gl(igtongl(ig)) ) .eq.  8 ) then
      vs (ig) =  v8
    elseif ( int ( gl(igtongl(ig)) ) .eq. 11 ) then
      vs (ig) =  v11
    endif
  enddo
  !
  ! the structure factor
  !
  do ig = 1, ngm
    arg = twopi * ( g(1,ig) * tau( 1, 1) + g(2,ig) * tau( 2, 1) + g(3,ig) * tau( 3, 1) )
    ss (ig) = cos ( arg )
  enddo
  !
  ! the empirical pseudopotential in real space
  ! for further use in h_psi
  !
  vr = czero
  do ig = 1, ngm
    vr ( nl ( ig ) ) = dcmplx ( ss (ig) * vs (ig), zero )
  enddo
  call cfft3 ( vr, nr1, nr2, nr3,  1)
  !

  !
  !  generate uniform {k} and {k+q} grids - no symmetry-reduction for now
  !  (MP method)
  !
  count = 0
  do i = 1, nq1
    ui = (q1 + 2.d0 * i - nq1 - 1.d0) / (2.d0 * nq1)
    do j = 1, nq2
      uj = (q2 + 2.d0 * j - nq2 - 1.d0) / (2.d0 * nq2)
      do k = 1, nq3
        uk = (q3 + 2.d0 * k - nq3 - 1.d0) / (2.d0 * nq3)
        count = count + 1
        xk (:, count) = ui * bg(:,1) + uj * bg(:,2) + uk * bg(:,3)
      enddo
    enddo
  enddo
  !
  ! include spin degeneracy
  wk = 2.d0 / float ( count )
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
! write(6,'(/4x,a)') repeat('-',67)
! write(6,'(4x,"Uniform k-point grid for the screened coulomb interaction"/)') 
! do ik = 1, nks
!    write ( 6, '(4x,"k(",i4,") = (",3f12.7,"), wk =",f12.7)') &
!        ik, (xk (ipol, ik) , ipol = 1, 3) , wk (ik)
! enddo
! write(6,'(4x,a/)') repeat('-',67)
  !
  do ik = 1, nks
    !
    if (ik/100*100-ik.eq.0) write(6,*) ik,nks
    !
    !  find wavefunctions for this kpoint
    !
    xxk = xk(:, ik)
    !
    ! the k-dependent kinetic energy in Ry
    !
    do ig = 1, ngm
      g2kin(ig) = ( (xxk(1) + g(1,ig))**2.d0 + &
                    (xxk(2) + g(2,ig))**2.d0 + &
                    (xxk(3) + g(3,ig))**2.d0 ) * tpiba2
    enddo
    !
    ! starting guess from direct diagonalization - low cutoff
    !
    do ig = 1, ngm0
      hk ( ig, ig) = g2kin ( ig )
    enddo
    !
    do ig1 = 1, ngm0
     do ig2 = ig1+1, ngm0
       !
       ! define ideltag
       ! 
       deltag = g(:,ig1) - g(:,ig2)
       ideltag = 1 
       do while ( .not. equiv ( deltag, g(:,ideltag) ) .and. ideltag.lt.ngm )
         ideltag = ideltag + 1
       enddo   
       !
       if (ideltag.ne.ngm) then
         hk ( ig1, ig2) = ss (ideltag) * vs (ideltag) 
         hk ( ig2, ig1) = hk ( ig1, ig2) 
       endif
       !
     enddo
    enddo
    !
    call rs ( ngm0, ngm0, hk, enval, 1, u, fv1, fv2, ierr)
    psi = czero
    do ibnd = 1, nbnd
     do ig = 1, ngm0
        psi(ig,ibnd) = dcmplx ( u(ig,ibnd), zero)
     enddo
    enddo
    !
    ! conjugate gradients diagonalization
    !
    precondition = max( 1.d0, g2kin )
    !
    call ccgdiagg (ngm, ngm, nbnd, psi, et(:,ik), precondition, eps, &
       maxter, .true., notconv, avg_iter, g2kin, vr)
    !
    !
!   write(6,'(4x,i3,8(2x,f6.3))') ik, et( 1: 8,ik)*ryd2ev
!   write(6,'(4x,3x,8(2x,f6.3))')     et( 9:16,ik)*ryd2ev
!   write(6,'(4x,3x,8(2x,f6.3))')     et(17:24,ik)*ryd2ev
!   write(6,'(4x,a/)') repeat('-',67)
    !
    !  direct write to file
    !
    write ( iunwfc, rec = ik, iostat = ios) psi 
    !
    ! to read: read ( iunwfc, rec = ik, iostat = ios) psi
    !
  enddo
  ! 
  do ik = 1, nksq
     !
     ikk = 2 * ik - 1
     ikq = ikk + 1
     !
     ! reads unperturbed wavefuctions u(k) and u(k+q)
     !
     read ( iunwfc, rec = ikk, iostat = ios) evc
     read ( iunwfc, rec = ikq, iostat = ios) evq
     !
     ! calculate the matrix elements <c,k+q|v,k> (G=0,G'=0)
     !
     do ibnd = 1, nbnd_occ
       do jbnd = nbnd_occ+1, nbnd
          a(ik,ibnd,jbnd) =  ZDOTC(ngm, evq(:,jbnd), 1, evc(:,ibnd), 1)
       enddo
     enddo
     !
  enddo
  !
  ! the dielectric matrix (G=0,G'=0,w=0)
  !
  ! the factor 2 takes into accout the (c,v) and (v,c) orderings of the mat.
  ! elements. the crystal volume & spin factor is 2/(nksq*omega), with 2/nks
  ! being given by wk(ikk)
  !
  qg2 = xxq(1)**2.d0 + xxq(2)**2.d0 + xxq(3)**2.d0
  vqg = fpi * e2 / omega / (tpiba2 * qg2) 
  !
  do iw = 1, nw
    !
    w = fmin + ( wmax - wmin ) * float(iw-1) / float(nw-1) / ryd2ev
    !
    chi0 = czero
    do ik = 1, nksq
       !
       ikk = 2 * ik - 1
       ikq = ikk + 1
       !
       do ibnd = 1, nbnd_occ
         do jbnd = nbnd_occ+1, nbnd 
            chi0 = chi0 + 2.d0 * wk(ikk) * 0.5d0                           &
              *   conjg ( a(ik,ibnd,jbnd) ) * a(ik,ibnd,jbnd) *            &
                ( cone / ( et (ibnd,ikk) - et (jbnd,ikq) + ci * eta - w)   &
                + cone / ( et (ibnd,ikk) - et (jbnd,ikq) + ci * eta + w) ) 
         enddo
       enddo
    enddo
    !
!   epsil  = dreal ( one - vqg * chi0 )
    epsil  = dreal ( one / ( one - vqg * chi0 ) ) 
    !
!   write(6,'(2f12.5)') xxq(1), epsil
    !
    write(6, *)  w, epsil
    write(16,*)  w, epsil
    !
  enddo
  ! 
  deallocate ( xk, wk )
  deallocate ( enval, hk, u, fv1, fv2, ss, vs )
  !
  return
  end subroutine dielec_mat
  !----------------------------------------------------------------
  !
  !
  !-----------------------------------------------------------------------
  subroutine solve_linter ( dvbare, dvscf, xxq, et, vr, w)
  !-----------------------------------------------------------------------
  !
  ! this works for one perturbation at a time
  !
  !-----------------------------------------------------------------------
  !
  use parameters
  use constants
  use gspace
  use kspace
  implicit none
  !
  complex(kind=DP) :: dvbare(nr)
  ! the perturbation in real space
  complex(kind=DP) :: dvscf(nr), vr(nr)
  real(DP) :: et(nbnd_occ, nks), w
  !
  complex(kind=DP) :: vr_dyn (nr)
  ! local potential plus teh dynamicla part w + i * eta
  !
  integer :: ik, ikk, ikq, iter, ibnd, jbnd, ios, ig, ir
  complex(kind=DP) :: evc (ngm, nbnd_occ), evq (ngm, nbnd_occ)
  real(kind=DP) :: g2kin(ngm), dr2, wgt, qg2, xxq(3)
  !
  complex(kind=DP) :: dpsi(ngm,nbnd_occ), dvpsi(ngm,nbnd_occ), &
                      dvpsi0(ngm,nbnd_occ), dvscfin(nr), hpsi(ngm), &
                      hpsi2(ngm,nbnd_occ), ps2(nbnd_occ,nbnd_occ), &
                      dpsi0(ngm,nbnd_occ), drhoscf(nr), dvscfout(nr)
  complex(kind=DP) :: aux(nr), aux1(nr), aux2(nr)
  complex(kind=DP) :: ps(nbnd_occ), auxg(ngm)
  complex(kind=DP) :: ZDOTC
  real(DP) :: eprec(nbnd_occ), h_diag(ngm, nbnd_occ), anorm, meandvb
  logical :: conv_root, convt
  integer :: lter
  external ch_psi_all, ch_psi_all_eta

  !
  ! include the dynamnical part inside the local potential 
  ! (which is already complex)
  !
  vr_dyn = vr - w 
  !

  !
  !  loop over the iterations
  !
  iter = 0
  convt = .false.
  do while (iter.lt.nmax_iter .and. .not.convt)
     !
     iter = iter + 1
     drhoscf = czero
     !
     do ik = 1, nksq
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        ! reads unperturbed wavefuctions u(k) and u(k+q)
        !
        read ( iunwfc, rec = ikk, iostat = ios) evc
        read ( iunwfc, rec = ikq, iostat = ios) evq
        !
        ! compute the kinetic energy for k+q
        !
        do ig = 1, ngm
          g2kin(ig) = ( (xk(1,ikq) + g(1,ig))**2.d0 + &
                        (xk(2,ikq) + g(2,ig))**2.d0 + &
                        (xk(3,ikq) + g(3,ig))**2.d0 ) * tpiba2
        enddo
        !
        if (iter.eq.1) then
          !
          do ibnd = 1, nbnd_occ
            !
            !  dpsi and dvscfin are set to zero
            !
            dpsi = czero
            dvscfin = czero
            !
            ! dvbare*psi is calculated for this k point and all bands...
            ! (we compute the product in real space)
            !
            aux = czero
            do ig = 1, ngm
              aux ( nl ( ig ) ) = evc (ig, ibnd)
            enddo
            call cfft3 ( aux, nr1, nr2, nr3,  1)
            do ir = 1, nr
              aux (ir) = aux(ir) * dvbare (ir)
            enddo
            ! back to G-space (fft order of G-vectors)
            call cfft3 ( aux, nr1, nr2, nr3, -1)
            ! switch to magnitude-order of G-vectors
            do ig = 1, ngm
              dvpsi(ig, ibnd) = aux( nl(ig) ) 
            enddo
            !
          enddo
          !
          ! writes dvpsi for this k-point on iunit iubar
          !
          write ( iubar, rec = ik, iostat = ios) dvpsi
          !
        else
          !
          ! read  dvbare*psi for this k-point on iunit iubar
          !
          read ( iubar, rec = ik, iostat = ios) dvpsi
          !
          ! dvpsi =  dvbare*psi + dvscfin*psi 
          !
          do ibnd = 1, nbnd_occ
            !
            aux = czero
            do ig = 1, ngm
              aux ( nl ( ig ) ) = evc (ig, ibnd)
            enddo
            call cfft3 ( aux, nr1, nr2, nr3,  1)
            do ir = 1, nr
               aux (ir) = aux (ir) * dvscfin (ir)
            enddo
            call cfft3 ( aux, nr1, nr2, nr3, -1)
            do ig = 1, ngm
              dvpsi (ig, ibnd) = dvpsi (ig, ibnd) + aux ( nl(ig) )
            enddo
            !
          enddo
          !
          !  read dpsi for this k-point from iudwf
          !
          read ( iudwf, rec = ik, iostat = ios) dpsi
          !
        endif
        !
        !  ( 1 - P_occ^{k+q} ) * dvpsi
        !
        do ibnd = 1, nbnd_occ
           auxg = czero
           do jbnd = 1, nbnd_occ
              ps(jbnd) = - ZDOTC(ngm, evq(:,jbnd), 1, dvpsi(:,ibnd), 1)
              call ZAXPY (ngm, ps (jbnd), evq (:, jbnd), 1, auxg, 1)
           enddo
           call DAXPY (2 * ngm, one, auxg, 1, dvpsi (:, ibnd), 1)
        enddo
        !
        !  change the sign of the known term
        !
        call DSCAL (2 * ngm * nbnd_occ, - 1.d0, dvpsi, 1)
        !
        ! iterative solution of the linear system 
        ! (H-et)*dpsi=   - ( 1 - P_occ^{k+q} ) * (dvbare+dvscf)*psi
        !              [                  dvpsi                     ]
!
!       !
!       ! this preconditioner sounds complicated...
!       !
!       do ibnd = 1, nbnd_occ
!          do ig = 1, ngm
!             auxg (ig) = g2kin (ig) * evq (ig, ibnd)
!          enddo
!          eprec (ibnd) = 1.35d0 * ZDOTC (ngm, evq (1, ibnd), 1, auxg, 1)
!       enddo
!       do ibnd = 1, nbnd_occ
!          do ig = 1, ngm
!             h_diag(ig,ibnd)=1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd))
!          enddo
!       enddo
!
! this seems to work even better...
!
        h_diag=1.d0 
!
        !
        ! now obtain dpsi from cgsolve_all
        !
!@      call cgsolve_all (ch_psi_all, et(:,ikk), dvpsi, dpsi, h_diag, & 
!@           ngm, ngm, tr_cgsolve, ik, lter, conv_root, anorm, nbnd_occ, &
!@           g2kin, vr_dyn, evq )
        call bcgsolve_all (ch_psi_all_eta, et(:,ikk), dvpsi, dpsi, h_diag, &
             ngm, ngm, tr_cgsolve, ik, lter, conv_root, anorm, nbnd_occ, &
             g2kin, vr_dyn, evq )

!       if (.not.conv_root) &
!          write( 6, '(4x,"ik",i4," linter: one or more roots not converged ",e10.3)') &
!          ik , anorm
!       write(6,'("cgsolve_all:",2x,3i5,3x,e9.3)') iter, ik, lter, anorm
        !
!       !
!       ! DEBUG: calculate (H-et+alpha*Pv)*dpsi-dvpsi, this should be zero 
!       ! if dpsi is the correct solution - O.K. this is checked
!       !
!       call ZGEMM ('C', 'N', nbnd_occ , nbnd_occ, ngm, (1.d0, 0.d0) , evq, &
!         ngm, dpsi, ngm, (0.d0, 0.d0) , ps2, nbnd_occ)
!       call ZGEMM ('N', 'N', ngm, nbnd_occ, nbnd_occ, dcmplx(alpha_pv,0.d0), evq, &
!         ngm, ps2, nbnd_occ, czero, hpsi2, ngm)
!       do ibnd = 1, nbnd_occ
!         call h_psi_c ( dpsi(:,ibnd), hpsi, g2kin, vr_dyn)
!         hpsi = hpsi - et(ibnd,ikk) * dpsi(:,ibnd) - dvpsi(:,ibnd)
!         hpsi = hpsi +  hpsi2 (:, ibnd)
!         write(6,*) '--------------------'
!         do ig = 1, ngm
!           write(6,'(3i5,3(2x,2f15.5))') &
!             ik, ibnd, ig, 100*dpsi(ig,ibnd), 100*dvpsi(ig,ibnd), 100*hpsi(ig)
!         enddo
!       enddo
        !
        ! writes dpsi for this k point on iunit iudwf
        !
        write ( iudwf, rec = ik, iostat = ios) dpsi
        !
        ! contribution to drhoscf from this kpoint
        !
        wgt = 2.d0 * wk(ikk) / omega
        !
        do ibnd = 1, nbnd_occ
          !
          aux1 = czero
          do ig = 1, ngm
            aux1 ( nl ( ig ) ) = evc (ig, ibnd)
          enddo
          call cfft3 ( aux1, nr1, nr2, nr3,  1)
          !
          aux2 = czero
          do ig = 1, ngm
            aux2 ( nl ( ig ) ) = dpsi (ig, ibnd)
          enddo
          call cfft3 ( aux2, nr1, nr2, nr3,  1)
          !
          do ir = 1, nr
            drhoscf (ir) = drhoscf (ir) + wgt * conjg (aux1 (ir) ) * aux2 (ir)
          enddo
          !
        enddo
        !
     enddo 
     !
     ! here we have drhoscf for this iteration
     ! compute the corresponding Hartree potential (RPA)
     !
     dvscfout = czero
     call cfft3 ( drhoscf, nr1, nr2, nr3, -1)
     !
     ! here we enforce zero average variation of the charge density 
     ! if the bare perturbation does not have a constant term
     ! (otherwise the numerical error, coupled with a small denominator
     ! in the coulomb term, gives rise to a spurious dvscf response)
     !
     meandvb = sqrt ( (sum(dreal(dvbare)))**2.d0 + (sum(aimag(dvbare)))**2.d0 ) / float(nr)
     if (meandvb.lt.1.d-8) drhoscf ( nl(1) ) = 0.d0
     !
     do ig = 1, ngm
       qg2 = (g(1,ig)+xxq(1))**2 + (g(2,ig)+xxq(2))**2 + (g(3,ig)+xxq(3))**2
       if (qg2 > 1.d-8) &
         dvscfout ( nl(ig) ) =  e2 * fpi * drhoscf ( nl(ig) ) / (tpiba2 * qg2)

     enddo
     !
     call cfft3 ( dvscfout, nr1, nr2, nr3,  1)
     !
     ! we mix with the old potential
     !
     ! modif broyden mixing, complex potential
     !
     call mix_potential_c ( nr, dvscfout, dvscfin, alpha_mix, dr2, tr2_ph, iter, nmix_ph, convt)
     !
     ! convergence criterion
     !
     convt = dr2.lt.tr2_ph
     !
     write(6,'(4x, "scf iteration ",i3,": dr2 = ",e8.2)') iter, dr2
     !
  enddo
  !
  ! at this point dvscfin is the converged Hartree screening.
  ! the screened coulomb interaction corresponds to dv_bare + dv_hartree (RPA)
  !
  dvscf = dvscfin + dvbare
  !
  return
  end subroutine solve_linter
  !
  !-----------------------------------------------------------------------
  !
  !
  !-----------------------------------------------------------------------
  subroutine solve_linter_notSCF ( dvbare, dvscf, xxq, et, vr, w)
  !-----------------------------------------------------------------------
  !
  ! this works for one perturbation at a time
  !
  !-----------------------------------------------------------------------
  !
  use parameters
  use constants
  use gspace
  use kspace
  implicit none
  !
  complex(kind=DP) :: dvbare(nr)
  ! the perturbation in real space
  complex(kind=DP) :: dvscf(nr), vr(nr)
  real(DP) :: et(nbnd_occ, nks), w
  !
  complex(kind=DP) :: vr_dyn (nr)
  ! the local potential plus the dynamical frequency w 
  !
  integer :: ik, ikk, ikq, iter, ibnd, jbnd, ios, ig, ir
  complex(kind=DP) :: evc (ngm, nbnd_occ), evq (ngm, nbnd_occ)
  real(kind=DP) :: g2kin(ngm), dr2, wgt, qg2, xxq(3)
  !
  complex(kind=DP) :: dpsi(ngm,nbnd_occ), dvpsi(ngm,nbnd_occ), &
                      dvpsi0(ngm,nbnd_occ), dvscfin(nr), hpsi(ngm), &
                      hpsi2(ngm,nbnd_occ), ps2(nbnd_occ,nbnd_occ), &
                      dpsi0(ngm,nbnd_occ), drhoscf(nr), dvscfout(nr)
  complex(kind=DP) :: aux(nr), aux1(nr), aux2(nr)
  complex(kind=DP) :: ps(nbnd_occ), auxg(ngm)
  complex(kind=DP) :: ZDOTC
  real(DP) :: eprec(nbnd_occ), h_diag(ngm, nbnd_occ), anorm, meandvb
  logical :: conv_root, convt
  integer :: lter
  external ch_psi_all, ch_psi_all_eta
  
   !
   ! include the dynamical part inside the local potential 
   ! (which is already complex)
   ! the imaginary component is going to be inside ch_psi_all_eta.f90
   !
   vr_dyn = vr - w 
   ! 

   !
   !  loop over the iterations
   !
   iter = 0
   convt = .false.
!* do while (iter.lt.nmax_iter .and. .not.convt)
     !
     iter = iter + 1
     drhoscf = czero
     !
     do ik = 1, nksq
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        ! reads unperturbed wavefuctions u(k) and u(k+q)
        !
        read ( iunwfc, rec = ikk, iostat = ios) evc
        read ( iunwfc, rec = ikq, iostat = ios) evq
        !
        ! compute the kinetic energy for k+q
        !
        do ig = 1, ngm
          g2kin(ig) = ( (xk(1,ikq) + g(1,ig))**2.d0 + &
                        (xk(2,ikq) + g(2,ig))**2.d0 + &
                        (xk(3,ikq) + g(3,ig))**2.d0 ) * tpiba2
        enddo
        !
        if (iter.eq.1) then
          !
          do ibnd = 1, nbnd_occ
            !
            !  dpsi and dvscfin are set to zero
            !
            dpsi = czero
            dvscfin = czero
            !
            ! dvbare*psi is calculated for this k point and all bands...
            ! (we compute the product in real space)
            !
            aux = czero
            do ig = 1, ngm
              aux ( nl ( ig ) ) = evc (ig, ibnd)
            enddo
            call cfft3 ( aux, nr1, nr2, nr3,  1)
            do ir = 1, nr
              aux (ir) = aux(ir) * dvbare (ir)
            enddo
            ! back to G-space (fft order of G-vectors)
            call cfft3 ( aux, nr1, nr2, nr3, -1)
            ! switch to magnitude-order of G-vectors
            do ig = 1, ngm
              dvpsi(ig, ibnd) = aux( nl(ig) ) 
            enddo
            !
          enddo
          !
          ! writes dvpsi for this k-point on iunit iubar
          !
          write ( iubar, rec = ik, iostat = ios) dvpsi
          !
        else
          !
          ! read  dvbare*psi for this k-point on iunit iubar
          !
          read ( iubar, rec = ik, iostat = ios) dvpsi
          !
          ! dvpsi =  dvbare*psi + dvscfin*psi 
          !
          do ibnd = 1, nbnd_occ
            !
            aux = czero
            do ig = 1, ngm
              aux ( nl ( ig ) ) = evc (ig, ibnd)
            enddo
            call cfft3 ( aux, nr1, nr2, nr3,  1)
            do ir = 1, nr
               aux (ir) = aux (ir) * dvscfin (ir)
            enddo
            call cfft3 ( aux, nr1, nr2, nr3, -1)
            do ig = 1, ngm
              dvpsi (ig, ibnd) = dvpsi (ig, ibnd) + aux ( nl(ig) )
            enddo
            !
          enddo
          !
          !  read dpsi for this k-point from iudwf
          !
          read ( iudwf, rec = ik, iostat = ios) dpsi
          !
        endif
        !
        !  ( 1 - P_occ^{k+q} ) * dvpsi
        !
        do ibnd = 1, nbnd_occ
           auxg = czero
           do jbnd = 1, nbnd_occ
              ps(jbnd) = - ZDOTC(ngm, evq(:,jbnd), 1, dvpsi(:,ibnd), 1)
              call ZAXPY (ngm, ps (jbnd), evq (:, jbnd), 1, auxg, 1)
           enddo
           call DAXPY (2 * ngm, one, auxg, 1, dvpsi (:, ibnd), 1)
        enddo
        !
        !  change the sign of the known term
        !
        call DSCAL (2 * ngm * nbnd_occ, - 1.d0, dvpsi, 1)
        !
        ! iterative solution of the linear system 
        ! (H-et)*dpsi=   - ( 1 - P_occ^{k+q} ) * (dvbare+dvscf)*psi
        !              [                  dvpsi                     ]
!
!       !
!       ! this preconditioner sounds complicated...
!       !
!       do ibnd = 1, nbnd_occ
!          do ig = 1, ngm
!             auxg (ig) = g2kin (ig) * evq (ig, ibnd)
!          enddo
!          eprec (ibnd) = 1.35d0 * ZDOTC (ngm, evq (1, ibnd), 1, auxg, 1)
!       enddo
!       do ibnd = 1, nbnd_occ
!          do ig = 1, ngm
!             h_diag(ig,ibnd)=1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd))
!          enddo
!       enddo
!
! this seems to work even better...
!
        h_diag=1.d0 
!
        !
        ! now obtain dpsi from cgsolve_all
        !
!@        call cgsolve_all (ch_psi_all, et(:,ikk), dvpsi, dpsi, h_diag, & 
!@             ngm, ngm, tr_cgsolve, ik, lter, conv_root, anorm, nbnd_occ, &
!@             g2kin, vr_dyn, evq )
        call bcgsolve_all (ch_psi_all_eta, et(:,ikk), dvpsi, dpsi, h_diag, & 
             ngm, ngm, tr_cgsolve, ik, lter, conv_root, anorm, nbnd_occ, &
             g2kin, vr_dyn, evq )
!       if (.not.conv_root) &
!          write( 6, '(4x,"ik",i4," linter: one or more roots not converged ",e10.3)') &
!          ik , anorm
!       write(6,'("cgsolve_all:",2x,3i5,3x,e9.3)') iter, ik, lter, anorm
        !
!       !
!       ! DEBUG: calculate (H-et+alpha*Pv)*dpsi-dvpsi, this should be zero 
!       ! if dpsi is the correct solution - O.K. this is checked
!       !
!       call ZGEMM ('C', 'N', nbnd_occ , nbnd_occ, ngm, (1.d0, 0.d0) , evq, &
!         ngm, dpsi, ngm, (0.d0, 0.d0) , ps2, nbnd_occ)
!       call ZGEMM ('N', 'N', ngm, nbnd_occ, nbnd_occ, dcmplx(alpha_pv,0.d0), evq, &
!         ngm, ps2, nbnd_occ, czero, hpsi2, ngm)
!       do ibnd = 1, nbnd_occ
!         call h_psi_c ( dpsi(:,ibnd), hpsi, g2kin, vr_dyn)
!         hpsi = hpsi - et(ibnd,ikk) * dpsi(:,ibnd) - dvpsi(:,ibnd)
!         hpsi = hpsi +  hpsi2 (:, ibnd)
!         write(6,*) '--------------------'
!         do ig = 1, ngm
!           write(6,'(3i5,3(2x,2f15.5))') &
!             ik, ibnd, ig, 100*dpsi(ig,ibnd), 100*dvpsi(ig,ibnd), 100*hpsi(ig)
!         enddo
!       enddo
        !
        ! writes dpsi for this k point on iunit iudwf
        !
        write ( iudwf, rec = ik, iostat = ios) dpsi
        !
        ! contribution to drhoscf from this kpoint
        !
        wgt = 2.d0 * wk(ikk) / omega
        !
        do ibnd = 1, nbnd_occ
          !
          aux1 = czero
          do ig = 1, ngm
            aux1 ( nl ( ig ) ) = evc (ig, ibnd)
          enddo
          call cfft3 ( aux1, nr1, nr2, nr3,  1)
          !
          aux2 = czero
          do ig = 1, ngm
            aux2 ( nl ( ig ) ) = dpsi (ig, ibnd)
          enddo
          call cfft3 ( aux2, nr1, nr2, nr3,  1)
          !
          do ir = 1, nr
            drhoscf (ir) = drhoscf (ir) + wgt * conjg (aux1 (ir) ) * aux2 (ir)
          enddo
          !
        enddo
        !
     enddo 
     !
     ! here we have drhoscf for this iteration
     ! compute the corresponding Hartree potential (RPA)
     !
     dvscfout = czero
     call cfft3 ( drhoscf, nr1, nr2, nr3, -1)
     !
     ! here we enforce zero average variation of the charge density 
     ! if the bare perturbation does not have a constant term
     ! (otherwise the numerical error, coupled with a small denominator
     ! in the coulomb term, gives rise to a spurious dvscf response)
     !
     meandvb = sqrt ( (sum(dreal(dvbare)))**2.d0 + (sum(aimag(dvbare)))**2.d0 ) / float(nr)
     if (meandvb.lt.1.d-8) drhoscf ( nl(1) ) = 0.d0
     !
     do ig = 1, ngm
       qg2 = (g(1,ig)+xxq(1))**2 + (g(2,ig)+xxq(2))**2 + (g(3,ig)+xxq(3))**2
       if (qg2 > 1.d-8) &
         dvscfout ( nl(ig) ) =  e2 * fpi * drhoscf ( nl(ig) ) / (tpiba2 * qg2)

     enddo
     !
     call cfft3 ( dvscfout, nr1, nr2, nr3,  1)
     !
     ! we mix with the old potential
     !
     ! modif broyden mixing, complex potential
     !
!*   call mix_potential_c ( nr, dvscfout, dvscfin, alpha_mix, dr2, tr2_ph, iter, nmix_ph, convt)
     !
     ! convergence criterion
     !
     convt = dr2.lt.tr2_ph
     !
     write(6,'(4x, "scf iteration ",i3,": dr2 = ",e8.2)') iter, dr2
     !
!* enddo
   !
   ! at this point dvscfin is the converged Hartree screening.
   ! the screened coulomb interaction corresponds to dv_bare + dv_hartree (RPA)
   !
!* dvscf = dvscfin + dvbare
   !
   ! this is  - v_c * chi_0
   dvscf = - dvscfout 
   !
   return
   end subroutine solve_linter_notSCF
   !
   !-----------------------------------------------------------------------
   !
