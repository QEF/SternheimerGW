  !
  !----------------------------------------------------------------
  program gwhs
  !----------------------------------------------------------------
  ! 
  ! pilot code the GW-HS method
  ! empirical pseudopotential method implementation
  ! for nonpolar tetrahedral semiconductors 
  !
  ! NOTE: I use only one G-grid. The reason is that, while V(G-G')
  ! has in principle components 2G_max, in the EPS scheme the  
  ! largest value happens for G2=11, therefore we can simply use
  ! the smooth grid...
  ! To perform a true scaling test we need to include the double
  ! grid (and count the operations on the smooth and on the dense
  ! grids).
  !
  ! IMPORTANT: every wfs has the same G-vects, I am not changing
  ! the cutoff for every k-point: |G|^2 is used instead of |k+G|^2
  ! so mind when using k beyond the first BZ
  !
  !----------------------------------------------------------------
  !
  use gspace
  use parameters
  use constants
  use kspace
#ifdef __PARA
  USE para
  USE mp, ONLY: mp_bcast, mp_barrier, mp_end
  USE io_global, ONLY: ionode_id
  USE mp_global,  ONLY : nproc, mpime, nproc_pool, my_pool_id, me_pool
#endif
  implicit none
  !
  ! variables
  !
  integer :: root = 0 ! root node for broadcast
  integer :: ig, ik, ig1, ig2, iq, nk, ik0, i, j, k, ios
  integer :: iw, iwp, iw0, iw0pw, iw0mw, count, ipol, ir
  integer :: recl, unf_recl, irp, igp, rec0, igstart, igstop
  integer :: ngpool, ngr, igs
  integer :: nwcoul, nwgreen, nwalloc, nwsigma
  integer, allocatable :: ind_w0mw (:,:), ind_w0pw (:,:)
  real(dbl) :: ui, uj, uk, wgreenmin, wgreenmax, w0mw, w0pw
  real(dbl) :: gcutm, gcutms, arg, wp
  real(dbl) :: k0mq(3), kplusg(3), xk0(3,nk0), xxq(3), xq0(3)
  real(dbl), allocatable :: ss(:), vs(:)
  real(DP), parameter :: eps8 = 1.0D-8
  real(DP) :: et(nbnd_occ, nq)
  real(DP), allocatable :: g2kin (:)
  real(DP), allocatable :: wtmp(:), wcoul(:), wgreen(:), wsigma(:), w_ryd(:)
  complex(DP) :: cexpp, cexpm
  complex(DP), allocatable :: psi(:,:)
  complex(dbl), allocatable :: vr(:), aux(:)
  complex(dbl), allocatable :: scrcoul (:,:,:), greenf (:,:,:), sigma(:,:,:)
  complex(dbl), allocatable :: scrcoul_g (:,:,:), greenf_g (:,:,:), sigma_g(:,:,:)
  logical :: allowed, foundp, foundm
  CHARACTER (LEN=9)   :: code = 'GWHS'
  CHARACTER(len=3)  :: nd_nmbr = '000' ! node number (used only in parallel case)
  character (len=3) :: nd_nmbr0
  CHARACTER (LEN=6) :: version_number = '0.4.8'
  character (len=256) :: wfcfile
  !
  !
  ! start serial code OR initialize parallel environment
  !
  CALL startup( nd_nmbr, code, version_number )
  !
#ifdef __PARA
  if (me.ne.1.or.mypool.ne.1) open (unit=stdout,file='/dev/null',status='unknown')
#endif

  call start_clock ('GWHS')
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
  !
  !----------------------------------------------------------------
  ! FIND THE G-VECTORS FOR THE SMALL SIGMA CUTOFF
  !----------------------------------------------------------------
  !
  ! the G^2 cutoff in units of 2pi/a_0
  ! Note that in Ry units the kinetic energy is G^2, not G^2/2
  !
  gcutms = four * ecuts / tpiba2      
! gcutms = ecuts / tpiba2      
  !
  ! set the fft grid
  !
  ! estimate nr1 and check if it is an allowed value for FFT
  !
  nr1s = 1 + int (2 * sqrt (gcutms) * sqrt( at(1,1)**2 + at(2,1)**2 + at(3,1)**2 ) ) 
  nr2s = 1 + int (2 * sqrt (gcutms) * sqrt( at(1,2)**2 + at(2,2)**2 + at(3,2)**2 ) ) 
  nr3s = 1 + int (2 * sqrt (gcutms) * sqrt( at(1,3)**2 + at(2,3)**2 + at(3,3)**2 ) ) 
  !
  do while (.not.allowed(nr1s)) 
    nr1s = nr1s + 1
  enddo
  do while (.not.allowed(nr2s)) 
    nr2s = nr2s + 1
  enddo
  do while (.not.allowed(nr3s)) 
    nr3s = nr3s + 1
  enddo
  !
  call ggens ( gcutms )
  !
  !-----------------------------------------------------------------
  ! in the parallel case split the ngms G-vectors across pools
  !-----------------------------------------------------------------
  !
#ifdef __PARA
  !
  npool = nproc / nproc_pool
  if (npool.gt.1) then
    !
    ! number of g-vec per pool and reminder
    ngpool = ngms / npool
    ngr = ngms - ngpool * npool
    ! the reminder goes to the first ngr pools
    if ( my_pool_id < ngr ) ngpool = ngpool + 1
    !
    igs = ngpool * my_pool_id + 1
    if ( my_pool_id >= ngr ) igs = igs + ngr
    !
    !  the index of the first and the last g vec in this pool
    !
    igstart = igs
    igstop = igs - 1 + ngpool
    !
    write (stdout,'(/4x,"Max n. of PW perturbations per pool = ",i5)') igstop-igstart+1
    !
  else
#endif
    !
    igstart = 1
    igstop = ngms
    !
#ifdef __PARA
  endif
#endif
  !
  !----------------------------------------------------------------
  !
#ifdef __PARA
  if (me.eq.1.and.mypool.eq.1) then
#endif
    !
    ! read from file the k-points for the self-energy \Sigma(k)
    !
    open ( 44, file = "./klist.dat", form = 'formatted', status = 'unknown')
    do ik = 1, nk0
      read (44,*) xk0(1,ik), xk0(2,ik), xk0(3,ik) 
    enddo
    !
#ifdef __PARA
  endif
  !
  !  bcast everything to all nodes
  !
  call mp_bcast ( xk0, root)
#endif
  !
  ! construct the empirical pseudopotential
  !
  allocate ( ss(ngm), vs(ngm) )
  do ig = 1, ngm
    arg = twopi * ( g(1,ig) * tau( 1, 1) + g(2,ig) * tau( 2, 1) + g(3,ig) * tau( 3, 1) )
    ss (ig) = cos ( arg )
  enddo
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
  ! set to zero top of valence band by shifting the
  ! local potential
  vr = vr - eshift
  !

  !
  allocate ( g2kin (ngm) )
  allocate ( aux (nrs) )
  allocate ( xq(3,nq), wq(nq), eval_occ(nbnd_occ,nq) )
  allocate ( gmap(ngm,27) )
  !
  ! generate the uniform {q} grid for the Coulomb interaction
  ! no symmetry-reduction for now - uniform and Gamma-centered
  ! (I was going insane with the folding of the MP mesh, I am not sure
  ! it's self-contained)
  !
  count = 0
  do i = 1, nq1
    ui = (i - 1.d0) / float (nq1)
!@    ui = (q1 + 2.d0 * i - nq1 - 1.d0) / (2.d0 * nq1)
    do j = 1, nq2
      uj = (j - 1.d0) / float (nq2)
!@      uj = (q2 + 2.d0 * j - nq2 - 1.d0) / (2.d0 * nq2)
      do k = 1, nq3
        uk = (k - 1.d0) / float (nq3)
!@        uk = (q3 + 2.d0 * k - nq3 - 1.d0) / (2.d0 * nq3)
        count = count + 1
        xq (:, count) = ui * bg(:,1) + uj * bg(:,2) + uk * bg(:,3)
      enddo
    enddo
  enddo
  wq = one / float ( count )
  if (count.ne.nq) call error ('gwhs','q-point count',count)
  ! the {k} grid is taken to coincide with the {q} grid
  ! nks = 2 * nq
  allocate ( xk (3,nks), wk(nks) )
  !
  write(stdout,'(/4x,a)') repeat('-',67)
  write(stdout,'(4x,"Uniform q-point grid for the screened Coulomb interaction"/)') 
  do iq = 1, nq
     write ( stdout, '(4x,"q(",i3," ) = (",3f12.7," ), wq =",f12.7)') &
         iq, (xq (ipol, iq) , ipol = 1, 3) , wq (iq)
  enddo
  write(stdout,'(4x,a/)') repeat('-',67)
  !
  ! generate the occupied eigenstates on the uniform grid
  ! this will be needed for the screened Coulomb below
  !
  recl = 2 * nbnd_occ * ngm  ! 2 stands for complex
  unf_recl = DIRECT_IO_FACTOR * recl
  wfcfile = './silicon'//'.wfc'
#ifdef __PARA
  call set_ndnmbr ( mypool, me_pool, nprocp, npool, nd_nmbr0)
  wfcfile =  trim(wfcfile)//'.'//nd_nmbr0
#endif
  !
  open ( iunwfc, file = wfcfile, iostat = ios, form = 'unformatted', &
       status = 'unknown', access = 'direct', recl = unf_recl)
  !
  write(stdout,'(/4x,a)') repeat('-',67)
  write(stdout,'(4x,"Occupied eigenvalues (eV)")') 
  write(stdout,'(4x,a/)') repeat('-',67)
  !
  allocate ( psi (ngm, nbnd_occ) )
  do iq = 1, nq
    !
    ! the k-dependent kinetic energy in Ry
    ! [Eq. (14) of Ihm,Zunger,Cohen J Phys C 12, 4409 (1979)]
    !
    do ig = 1, ngm
      kplusg = xq(:, iq) + g(:,ig)
      g2kin ( ig ) = tpiba2 * dot_product ( kplusg, kplusg )
    enddo
    !
    call eigenstates2 ( xq(:, iq), vr, g2kin, psi, eval_occ(:,iq) ) 
    !
    !  direct write to file - take into account the k/k+q alternation
    !
    write ( iunwfc, rec = 2 * iq - 1, iostat = ios) psi
    !
    write ( stdout, '(4x,"k(",i3," )",10(3x,f7.3))') iq, eval_occ(:,iq)*ryd2ev
    !
  enddo
  deallocate( psi )
  write(stdout,'(4x,a/)') repeat('-',67)
  !
  ! here we generate the G-map for the folding into the first BZ
  !
  call refold ( )
  !
  !  ------------------------------------------------
  !  MAIN: CALCULATE G AND W AND PERFORM CONVOLUTIONS
  !  ------------------------------------------------
  !
  ! ----------------------------------------------------------------
  ! generate frequency bins 
  ! ----------------------------------------------------------------
  !
  ! Here I assume Sigma is needed for w0 between wsigmamin and wsigmamax
  ! The convolution requires W for positive frequencies w up to wcoulmax
  ! (even function - cf Shishkin and Kress) and the GF spanning w0+-w.
  ! Therefore the freq. range of GF is 
  ! from (wsigmamin-wcoulmax) to (wsigmamax+wcoulmax)
  ! the freq. dependence of the GF is inexpensive, so we use the same spacing
  ! 
  ! NB: I assume wcoulmax>0, wsigmamin=<0, wsigmamax>0 and zero of energy at the Fermi level
  !
  wgreenmin = wsigmamin-wcoulmax
  wgreenmax = wsigmamax+wcoulmax
  !
  nwalloc = 1 + ceiling( (wgreenmax-wgreenmin)/deltaw )
  allocate(wtmp(nwalloc), wcoul(nwalloc), wgreen(nwalloc), wsigma(nwalloc), w_ryd(nwalloc) )
  wcoul = zero
  wgreen = zero
  wsigma = zero
  !
  do iw = 1, nwalloc
    wtmp(iw) = wgreenmin + (wgreenmax-wgreenmin)/float(nwalloc-1)*float(iw-1)
  enddo
  ! align the bins with the zero of energy
  wtmp = wtmp - minval ( abs ( wgreen) )
  !
  nwgreen = 0
  nwcoul = 0
  nwsigma = 0
  !
  do iw = 1, nwalloc
    if ( ( wtmp(iw) .ge. wgreenmin ) .and. ( wtmp(iw) .le. wgreenmax) ) then
       nwgreen = nwgreen + 1
       wgreen(nwgreen) = wtmp(iw)
    endif
    if ( ( wtmp(iw) .ge. zero ) .and. ( wtmp(iw) .le. wcoulmax) ) then
       nwcoul = nwcoul + 1
       wcoul(nwcoul) = wtmp(iw)
    endif
    if ( ( wtmp(iw) .ge. wsigmamin ) .and. ( wtmp(iw) .le. wsigmamax) ) then
       nwsigma = nwsigma + 1
       wsigma(nwsigma) = wtmp(iw)
    endif
  enddo
  ! 
  ! now find the correspondence between the arrays
  ! This is needed for the convolution G(w0-w)W(w) at the end
  !
  allocate ( ind_w0mw (nwsigma,nwcoul), ind_w0pw (nwsigma,nwcoul) )
  !
  do iw0 = 1, nwsigma
    do iw = 1, nwcoul
      !
      w0mw = wsigma(iw0)-wcoul(iw)
      w0pw = wsigma(iw0)+wcoul(iw)
      !
      foundp = .false.
      foundm = .false.
      !
      do iwp = 1, nwgreen
        if ( abs(w0mw-wgreen(iwp)) .lt. 1.d-10 ) then
          foundm = .true.
          iw0mw = iwp
        endif
        if ( abs(w0pw-wgreen(iwp)) .lt. 1.d-10 ) then
          foundp = .true.
          iw0pw = iwp
        endif
      enddo
      !
      if ( ( .not. foundm ) .or. ( .not. foundp ) ) then
         call errore ('gwhs','frequency correspondence not found',1)
      else
         ind_w0mw(iw0,iw) = iw0mw 
         ind_w0pw(iw0,iw) = iw0pw 
      endif
      !
    enddo
  enddo
  !
  allocate ( scrcoul (nrs, nrs, nwcoul) )
  allocate ( greenf (nrs, nrs, nwgreen) )
  allocate ( sigma (nrs, nrs, nwsigma) )
  allocate ( scrcoul_g (ngms, ngms, nwcoul) )
  allocate ( greenf_g (ngms, ngms, nwgreen) )
  allocate ( sigma_g (ngms, ngms, nwsigma) )
  !
  ! prepare the unit to write the Coulomb potential
  ! each q-point is associated with one record
  !
! recl = 2 * nrs * nrs * nwcoul
  recl = 2 * ngms * ngms * nwcoul
  unf_recl = DIRECT_IO_FACTOR * recl
  open ( iuncoul, file = "./silicon.coul", iostat = ios, form = 'unformatted', &
       status = 'unknown', access = 'direct', recl = unf_recl)
  !
  ! prepare the unit to write the Green's function 
  ! each (k0-q)-point is associated with one record
  !
! recl = 2 * nrs * nrs * nwgreen
  recl = 2 * ngms * ngms * nwgreen
  unf_recl = DIRECT_IO_FACTOR * recl
  open ( iungreen, file = "./silicon.green", iostat = ios, form = 'unformatted', &
       status = 'unknown', access = 'direct', recl = unf_recl)
  !
  ! prepare the unit to write the self-energy 
  ! each k0-point is associated with one record
  !
  recl = 2 * ngms * ngms * nwsigma
  unf_recl = DIRECT_IO_FACTOR * recl
  open ( iunsigma, file = "./silicon.sigma", iostat = ios, form = 'unformatted', &
       status = 'unknown', access = 'direct', recl = unf_recl)
  !
  write(stdout,'(4x,"Screened Coulomb interaction:")')
  !
  ! loop over {q} for the screened Coulomb interaction
  !
!@  do iq = 1, nq
  do iq = 1, 2
    !
    write(stdout,'(4x,3x,"iq = ",i3)') iq
    scrcoul = czero
    !
    if (igstart.eq.1) then
      !
      ! In the case (q=0, G=0) we perform a separate
      ! calculation for scrcoul(ig=1,:,:)
      ! (in the parallel case: only the processor having the G=0 vec)
      !
      xq0 = (/ 0.01 , 0.00, 0.00 /) ! this should be set from input
      if ( ( xq(1,iq)*xq(1,iq) + xq(2,iq)*xq(2,iq) + xq(3,iq)*xq(3,iq) ) .lt. 1.d-10 ) &
      call coulomb_q0G0 ( vr, xq0, nwcoul, wcoul, scrcoul )
    endif
    !
    ! the grids {k} and {k+q} for the dVscf will be obtained
    ! by shuffling the {q} grid
    !
    call coulomb ( vr, xq(:,iq), nwcoul, wcoul, scrcoul, igstart, igstop )
    !
#ifdef __PARA
    !
    ! use poolreduce to bring together the results from each pool
    !
    call poolreduce ( 2 * nrs * nrs * nwcoul, scrcoul)
    !
    if (me.eq.1.and.mypool.eq.1) then
#endif
      !
      scrcoul_g = scrcoul(1:ngms,1:ngms,:)
      write ( iuncoul, rec = iq, iostat = ios) scrcoul_g
      !
#ifdef __PARA
    endif
#endif
    write (stdout,'(4x,"Written scrcoul for iq = ",i3)') iq
    !
  enddo 
  !
  write(stdout,'(4x,"Green''s function:")')
  ! loop over the {k0} set for the Self-Energy
  !
  do ik0 = 1, 1 !@ nk0
    !
    write(stdout,'(4x,"ik0 = ",i3)') ik0
    !
    ! loop over the {k0-q} grid for the Green's function
    !
!@    do iq = 1, nq
    do iq = 1, 2
      !
      write(stdout,'(4x,3x,"iq = ",i3)') iq
      greenf = czero
      !
      !  k0mq = k0 - q
      !
      k0mq = xk0(:,ik0) - xq(:,iq)
      !
      ! the k-dependent kinetic energy in Ry
      ! [Eq. (14) of Ihm,Zunger,Cohen J Phys C 12, 4409 (1979)]
      !
      do ig = 1, ngm
        kplusg = k0mq + g(:,ig)
        g2kin ( ig ) = tpiba2 * dot_product ( kplusg, kplusg )
      enddo
      !
      ! need to use multishift in green
      call green_linsys ( vr, g2kin, k0mq, nwgreen, wgreen, greenf, igstart, igstop )
      !
      ! now greenf contains the Green's function
      ! for this k0mq point, all frequencies, and in G-space
      !
#ifdef __PARA
      !
      ! use poolreduce to bring together the results from each pool
      !
      call poolreduce ( 2 * nrs * nrs * nwgreen, greenf)
      !
      if (me.eq.1.and.mypool.eq.1) then
#endif
        !
        rec0 = (ik0-1) * nq + (iq-1) + 1 
        greenf_g = greenf(1:ngms,1:ngms,:)
        write ( iungreen, rec = rec0, iostat = ios) greenf_g 
        !
#ifdef __PARA
      endif
#endif
      !
      ! end loop on {k0-q} and {q}
    enddo 
    !
    ! end loop on {k0}
  enddo 
  !
  ! G TIMES W PRODUCT
  !
  call start_clock ('GW product')
  !
  w_ryd = wcoul / ryd2ev
  do ik0 = 1, 1 !@ nk0 
    !
    write(stdout,'(4x,"Direct product GW for k0(",i3," ) = (",3f12.7," )")') &
      ik0, (xk0 (ipol, ik0) , ipol = 1, 3)
    !
    ! now sum over {q} the products G(k0-q)W(q) 
    !
    sigma = czero
    !
#ifdef __PARA
    ! only proc 0 reads from file and does the product
    ! (need some sort of parallelization here)
    if (me.eq.1.and.mypool.eq.1) then
#endif
!@    do iq = 1, nq
    do iq = 1, 2 
      !
      write(stdout,'(4x,"Summing iq = ",i4)') iq
      !
      ! go to R-space to perform the direct product with W
      ! both indeces are in SIZE order of G-vectors
      !
      rec0 = (ik0-1) * nq + (iq-1) + 1 
      read ( iungreen, rec = rec0, iostat = ios) greenf_g 
      greenf = czero
      greenf(1:ngms,1:ngms,:) = greenf_g
      !
      do iw = 1, nwgreen
        do ig = 1, ngms
          aux = czero
          do igp = 1, ngms
            aux(nls(igp)) = greenf(igp,ig,iw)
          enddo
          call cfft3s ( aux, nr1s, nr2s, nr3s,  1)
          greenf(:,ig,iw) = aux
        enddo
        do ir = 1, nrs
          aux = czero
          do igp = 1, ngms
            aux(nls(igp)) = greenf(ir,igp,iw)
          enddo
          call cfft3s ( aux, nr1s, nr2s, nr3s,  1)
          greenf(ir,:,iw) = aux
        enddo
      enddo
!     write(stdout,'(4x,"FFTW of Green''s function passed"/)') 
      !
      read ( iuncoul, rec = iq, iostat = ios) scrcoul_g
      scrcoul = czero
      scrcoul(1:ngms,1:ngms,:) = scrcoul_g
      !
      ! now scrcoul(nrs,nrs,nwcoul) contains the screened Coulomb
      ! interaction for this q point, all frequencies, and in G-space:
      ! go to R-space to perform the direct product with G
      ! both indeces are in SIZE order of G-vectors
      !
      do iw = 1, nwcoul
        do ig = 1, ngms
          aux = czero
          do igp = 1, ngms
            aux(nls(igp)) = scrcoul(igp,ig,iw)
          enddo
          call cfft3s ( aux, nr1s, nr2s, nr3s,  1)
          scrcoul(:,ig,iw) = aux
        enddo 
        do ir = 1, nrs
          aux = czero
          do igp = 1, ngms
            aux(nls(igp)) = scrcoul(ir,igp,iw)
          enddo
          call cfft3s ( aux, nr1s, nr2s, nr3s,  1)
          scrcoul(ir,:,iw) = aux
        enddo 
      enddo
!     write(stdout,'(4x,"FFTW of Coulomb passed"/)') 
      !
      ! combine Green's function and screened Coulomb ( sum_q wq = 1 )
      !
      do iw = 1, nwcoul
        !
        cexpm = exp ( -ci * eta * w_ryd(iw) )
        cexpp = exp (  ci * eta * w_ryd(iw) )
        !
        do iw0 = 1, nwsigma
          !
          iw0mw = ind_w0mw (iw0,iw)
          iw0pw = ind_w0pw (iw0,iw)
          !
          do ir = 1, nrs
            do irp = 1, nrs
              !
              sigma(ir,irp,iw0) = sigma(ir,irp,iw0) + wq (iq) * ci / twopi        &
                * ( cexpm * greenf(ir,irp,iw0mw) + cexpp * greenf(ir,irp,iw0pw) ) &
                * scrcoul (ir,irp,iw)
              !
            enddo
          enddo
          !
        enddo
        !
      enddo
      ! deltaw is in eV here - so if everythin is ok also sigma is in eV
      sigma = sigma * deltaw 
      !
      ! end loop on {k0-q} and {q}
    enddo 
    !
    ! Now we have summed over q in G(k0-q)W(q) and we can go back
    ! to G-space before calculating the sandwitches with the wavefunctions
    ! note: we go to SIZE order of G-vectors
    !
    do iw = 1, nwsigma
      do ir = 1, nrs
        aux = czero
        do irp = 1, nrs
          aux(irp) = sigma(ir,irp,iw)
        enddo
        call cfft3s ( aux, nr1s, nr2s, nr3s, -1)
        do ig = 1, ngms
          sigma (ir,ig,iw) = aux(nls(ig))
        enddo
      enddo
      do ig = 1, ngms
        aux = czero
        do irp = 1, nrs
          aux(irp) = sigma(irp,ig,iw)
        enddo
        call cfft3s ( aux, nr1s, nr2s, nr3s, -1)
        do igp = 1, ngms
          sigma (igp,ig,iw) = aux(nls(igp))
        enddo
      enddo
    enddo
    !
    ! everything beyond ngms is garbage
    !
    do ig = ngms+1, nrs
     do igp = ngms+1, nrs
      do iw = 1, nwsigma
         sigma (ig,igp,iw) = czero
      enddo
     enddo
    enddo
    !
    ! At this point, since we did FFT back and forth,
    ! the self-energy should be in (G,G') space within
    ! the same convention adopted for G and W and explained
    ! in the paper.
    !
#ifdef __PARA
    endif
    !
    ! use poolreduce to bring together the results from each pool
    !
    call poolreduce ( 2 * nrs * nrs * nwsigma, sigma)
    !
    if (me.eq.1.and.mypool.eq.1) then
#endif
      sigma_g = sigma(1:ngms,1:ngms,:)
      write ( iunsigma, rec = ik0, iostat = ios) sigma_g
#ifdef __PARA
    endif
#endif
    !
    ! end loop on {k0}
  enddo 
  ! 
  call stop_clock ('GW product')
  !
  ! CALCULATION OF THE MATRIX ELEMENTS
  !
  do ik0 = 1, 1 !@ nk0 
    call sigma_matel ( ik0, vr, xk0, nwsigma, wsigma)
  enddo
  !
!
! it looks like I have a problem in closing these
! files in parallel - tried several things (only headnode
! or everybody; only keep; keep some and delete some; all delete)
! should not matter that much as long as it finishes smoothly
!
! close (iuncoul, status = 'delete')
! close (iungreen, status = 'delete')
! close (iunsigma, status = 'keep')

  close (iunwfc, status = 'delete')
  !
  call stop_clock ('GWHS')
  !
  call print_clock('GWHS')
  call print_clock('coulomb')
  call print_clock('green_linsys')
  call print_clock('GW product')
  call print_clock('sigma_matel')
  !
  write(stdout,'(/4x,"End of program GWHS")')
  write(stdout,'(4x,a/)') repeat('-',67)
#ifdef __PARA
  call mp_barrier()
#endif
  stop
  !
  stop
  end program gwhs
  !----------------------------------------------------------------
  !
