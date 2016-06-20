!
! Copyright (C) 2001-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE sigmadens (filplot,plot_num)
  !-----------------------------------------------------------------------
  !      Writes the screened Coulomb interaction or Self-Energy Operator
  !      into a file format suitable for plotting
  !-----------------------------------------------------------------------
  !
  !      DESCRIPTION of the INPUT: see file INPUT_PP in Doc/
  !
  USE kinds,       ONLY : dp
  USE io_global,   ONLY : stdout, ionode, ionode_id
  USE io_files,    ONLY : nd_nmbr, prefix, diropn, tmp_dir
  USE mp_pools,    ONLY : nproc_pool
  USE mp_world,    ONLY : world_comm
  USE mp,          ONLY : mp_bcast
  USE parameters,  ONLY : ntypx
  USE constants,   ONLY : pi, fpi
  USE cell_base,   ONLY : at, bg, celldm, ibrav, alat, omega, tpiba, tpiba2
  USE ions_base,   ONLY : nat, ityp, atm, ntyp => nsp, tau, zv
  USE lsda_mod,    ONLY : nspin
  USE fft_base,    ONLY : scatter_grid, dfftp, dffts
  USE fft_interfaces,   ONLY : fwfft, invfft
  USE grid_subroutines, ONLY : realspace_grid_init
  USE gvect,         ONLY: ngm, nl, g, gcutm
  USE gvecs,         ONLY: gcutms, doublegrid, dual, ecuts 
  USE recvec_subs,   ONLY: ggen 
  USE wvfct,         ONLY: ecutwfc
  USE run_info,      ONLY: title
  USE control_flags, ONLY: gamma_only
  USE wavefunctions_module,  ONLY: psic
  USE symm_base,  ONLY : nsym, s, time_reversal, t_rev, ftau, invs, nrot
  USE klist,         ONLY : wk, xk, nkstot, nks, lgauss


  IMPLICIT NONE
  character (len=256), INTENT(in) :: filplot
  !
  ! If plot_num=-1 the dimensions and structural data are read from the charge
  ! or potential file, otherwise it uses the data already read from
  ! the files in outdir.
  !
  integer, INTENT(in) :: plot_num
  !
  integer, PARAMETER :: nfilemax = 7
  ! maximum number of files with charge

  integer :: ounit, iflag, ios, ipol, nfile, ifile, nx, ny, nz, &
       na, i, output_format, idum, direction

  real(DP) :: e1(3), e2(3), e3(3), x0 (3), radius, m1, m2, m3, &
       weight (nfilemax), isovalue,heightmin,heightmax

  real(DP), allocatable :: aux(:)

  character (len=256) :: fileout
  character (len=13), DIMENSION(0:7) :: formatname = &
       (/ 'gnuplot      ', &
          'contour.x    ', &
          'plotrho.x    ', &
          'XCrySDen     ', &
          'gOpenMol     ', &
          'XCrySDen     ', &
          'Gaussian cube', & 
          'gnuplot x,y,f' /)
  character (len=20), DIMENSION(0:4) :: plotname = &
       (/ '1D spherical average', &
          '1D along a line     ', &
          '2D contour          ', &
          '3D                  ', &
          '2D polar on a sphere'/)

  real(DP) :: celldms (6), gcutmsa, duals, zvs(ntypx), ats(3,3)
  real(DP), allocatable :: taus (:,:), rhor(:), rhos(:)

  integer :: ibravs, nr1sxa, nr2sxa, nr3sxa, nr1sa, nr2sa, nr3sa, &
       ntyps, nats
  integer, allocatable :: ityps (:)
  character (len=3) :: atms(ntypx)
  character (len=256) :: filepp(nfilemax)
  character (len=20) :: interpolation
  real(DP) :: rhotot
  complex(DP), allocatable:: rhog (:)
!rho or polarization in G space
  logical :: fast3d, isostm_flag
!\Sigma(\r,\r';\omega) variables.
  complex(DP), allocatable :: scrcoul_g(:,:,:), sigma(:,:), sigma_g(:,:,:), &
                              sigma_q(:,:,:), sigma_ex(:,:)
  complex(DP), allocatable :: sigmar(:)

!SIGMADENS
  integer   :: ngmpol, nwsigma, iunsigma, iqs, ir
  integer   :: isym, nqs, ig
  integer   :: iq, lrcoul
  integer*8 :: unf_recl
  integer, allocatable     :: gmapsym(:,:)
  complex(DP), allocatable ::  eigv(:,:)
  complex(DP), allocatable :: eigx (:), eigy (:), eigz (:)
  complex(DP), allocatable:: sigmag (:)
  complex(DP)   :: eigx0
  real(DP)      :: aq(3)
  character(len=256) :: tempfile, filename
!  USE disp,  ONLY : nq1, nq2, nq3, x_q, nqs, wq
  real(DP)      :: xq(3,64), wq(64)

#define DIRECT_IO_FACTOR 8 

  allocate ( gmapsym  (ngm, nrot)   )
  allocate ( eigv     (ngm, nrot)   )




  NAMELIST /plot/  &
       nfile, filepp, weight, iflag, e1, e2, e3, nx, ny, nz, x0, &
       radius, output_format, fileout, interpolation, &
       isostm_flag, isovalue, heightmin, heightmax, direction

  !
  !   set the DEFAULT values
  !
  nfile         = 1
  filepp(1)     = filplot
  weight(1)     = 1.0d0
  iflag         = 0
  radius        = 1.0d0
  output_format = -1
  fileout       = ' '

  e1(:)         = 0.d0
  e2(:)         = 0.d0
  e3(:)         = 0.d0
  x0(:)         = 0.d0
  nx            = 0
  ny            = 0
  nz            = 0
  interpolation = 'fourier'
  isostm_flag   = .false.
  isovalue      = 0.d0
  heightmin     = 0.0d0
  heightmax     = 1.0d0
  direction     = 1
  !
  !    read and check input data
  !
  ! reading the namelist 'plot'
  !
  if (ionode) READ (5, plot, iostat = ios)
  !
  CALL mp_bcast( ios, ionode_id, world_comm )
  CALL mp_bcast( nfile, ionode_id, world_comm )

  if (ios /= 0) THEN
     if (nfile > nfilemax) THEN
        ! if this happens the reading of the namelist will fail
        ! tell to user why
        CALL infomsg('chdens ', 'nfile is too large, exiting')
     ELSE
        CALL infomsg ('chdens', 'namelist plot not found or invalid, exiting')
     ENDif
     return
  ENDif

  CALL mp_bcast( filepp, ionode_id, world_comm )
  CALL mp_bcast( weight, ionode_id, world_comm )
  CALL mp_bcast( iflag, ionode_id, world_comm )
  CALL mp_bcast( radius, ionode_id, world_comm )
  CALL mp_bcast( output_format, ionode_id, world_comm )
  CALL mp_bcast( fileout, ionode_id, world_comm )
  CALL mp_bcast( e1, ionode_id, world_comm )
  CALL mp_bcast( e2, ionode_id, world_comm )
  CALL mp_bcast( e3, ionode_id, world_comm )
  CALL mp_bcast( x0, ionode_id, world_comm )
  CALL mp_bcast( nx, ionode_id, world_comm )
  CALL mp_bcast( ny, ionode_id, world_comm )
  CALL mp_bcast( nz, ionode_id, world_comm )
  CALL mp_bcast( interpolation, ionode_id, world_comm )
  CALL mp_bcast( isostm_flag, ionode_id, world_comm )
  CALL mp_bcast( isovalue, ionode_id, world_comm )
  CALL mp_bcast( heightmin, ionode_id, world_comm )
  CALL mp_bcast( heightmax, ionode_id, world_comm )
  CALL mp_bcast( direction, ionode_id, world_comm )  

  if (output_format == -1 .or. iflag == -1) THEN
     CALL infomsg ('chdens', 'output format not set, exiting' )
     return
  ENDif
  !
  ! check for number of files
  !
  if (nfile < 1 .or. nfile > nfilemax) &
       CALL errore ('chdens ', 'nfile is wrong ', 1)

  ! check for iflag make sure input parameters are logical

  if (iflag <= 1) THEN

     ! 1D plot : check variables

     if (e1(1)**2 + e1(2)**2 + e1(3)**2 < 1d-6) &
         CALL errore ('chdens', 'missing e1 vector', 1)
     if (nx <= 0 )   CALL errore ('chdens', 'wrong nx', 1)

  ELSEif (iflag == 2) THEN

     ! 2D plot : check variables

     if (e1(1)**2 + e1(2)**2 + e1(3)**2 <  1d-6 .or. &
         e2(1)**2 + e2(2)**2 + e2(3)**2 <  1d-6)     &
         CALL errore ('chdens', 'missing e1/e2 vectors', 1)
     if (abs(e1(1)*e2(1) + e1(2)*e2(2) + e1(3)*e2(3)) > 1d-6) &
         CALL errore ('chdens', 'e1 and e2 are not orthogonal', 1)
     if (nx <= 0 .or. ny <= 0 )   CALL errore ('chdens', 'wrong nx/ny', 2)

  ELSEif (iflag == 3) THEN

     ! 3D plot : check variables

     if ( abs(e1(1)*e2(1) + e1(2)*e2(2) + e1(3)*e2(3)) > 1d-6 .or. &
          abs(e1(1)*e3(1) + e1(2)*e3(2) + e1(3)*e3(3)) > 1d-6 .or. &
          abs(e2(1)*e3(1) + e2(2)*e3(2) + e2(3)*e3(3)) > 1d-6 )    &
         CALL errore ('chdens', 'e1, e2, e3 are not orthogonal', 1)

     if ((iflag==3) .and.(output_format < 3 .or. output_format > 6)) &
        CALL errore ('chdens', 'incompatible iflag/output_format', 1)
     if ((iflag/=3) .and. ((output_format == 5) .or. (output_format == 6))) &
        CALL errore ('chdens', 'output_format=5/6, iflag<>3', 1)

  ELSEif (iflag  == 4) THEN

     if (nx <= 0 .or. ny <= 0 )   CALL errore ('chdens', 'wrong nx/ny', 4)

  ELSE

     CALL errore ('chdens', 'iflag not implemented', 1)

  ENDif

  !END INPUT CHECK
  ! check interpolation
  if (trim(interpolation) /= 'fourier' .and. trim(interpolation) /= 'bspline') &
     call errore('chdens', 'wrong interpolation: ' // trim(interpolation), 1)

  ! if isostm_flag checks whether the input variables are set
  if (isostm_flag) THEN
     if (heightmax > 1.0 .or. heightmin > 1.0 .or. heightmin < 0.0 &
               .or. heightmax < 0.0 ) THEN
         CALL errore('isostm','problem with heightmax/min',1)
     ENDif
      
     if (direction /= 1 .and. direction /= -1) THEN
         CALL errore('isostm','direction not equal to +- 1',1)
     ENDif
  END if
  !
  ! Read the header and allocate objects
  !
  if (plot_num==-1) THEN
     if (ionode) &
        CALL read_io_header(filepp (1), title, dfftp%nr1x, dfftp%nr2x, &
                dfftp%nr3x, dfftp%nr1, dfftp%nr2, dfftp%nr3, nat, ntyp,&
                ibrav, celldm, at, gcutm, dual, ecutwfc, idum )
     CALL mp_bcast( title, ionode_id, world_comm )
     CALL mp_bcast( dfftp%nr1x, ionode_id, world_comm )
     CALL mp_bcast( dfftp%nr2x, ionode_id, world_comm )
     CALL mp_bcast( dfftp%nr3x, ionode_id, world_comm )
     CALL mp_bcast( dfftp%nr1, ionode_id, world_comm )
     CALL mp_bcast( dfftp%nr2, ionode_id, world_comm )
     CALL mp_bcast( dfftp%nr3, ionode_id, world_comm )
     CALL mp_bcast( nat, ionode_id, world_comm )
     CALL mp_bcast( ntyp, ionode_id, world_comm )
     CALL mp_bcast( ibrav, ionode_id, world_comm )
     CALL mp_bcast( celldm, ionode_id, world_comm )
     CALL mp_bcast( at, ionode_id, world_comm )
     CALL mp_bcast( gcutm, ionode_id, world_comm )
     CALL mp_bcast( dual, ionode_id, world_comm )
     CALL mp_bcast( ecutwfc, ionode_id, world_comm )
     !
     ! ... see comment above
     !
     allocate(tau (3, nat))
     allocate(ityp(nat))
     !
     CALL latgen (ibrav, celldm, at(1,1), at(1,2), at(1,3), omega )
     alat = celldm (1) ! define alat
     at = at / alat    ! bring at in units of alat
     tpiba = 2.d0 * pi / alat
     tpiba2 = tpiba**2
     doublegrid = dual>4.0d0
     if (doublegrid) THEN
        gcutms = 4.d0 * ecutwfc / tpiba2
     ELSE
        gcutms = gcutm
     ENDif

     nspin = 1

     CALL recips (at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
     CALL volume (alat, at(1,1), at(1,2), at(1,3), omega)
     CALL realspace_grid_init ( dfftp, at, bg, gcutm )
     CALL realspace_grid_init ( dffts, at, bg, gcutms)
  ENDif

  allocate  (rhor(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))

  allocate  (rhos(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
  allocate  (taus( 3 , nat))
  allocate  (ityps( nat))
  !
  rhor (:) = 0.0_DP
  !
  ! Read files, verify consistency
  ! Note that only rho is read; all other quantities are discarded
  !
  do ifile = 1, nfile
     !
     CALL plot_io (filepp (ifile), title, nr1sxa, nr2sxa, nr3sxa, &
          nr1sa, nr2sa, nr3sa, nats, ntyps, ibravs, celldms, ats, gcutmsa, &
          duals, ecuts, idum, atms, ityps, zvs, taus, rhos, - 1)

     if (ifile==1.and.plot_num==-1) THEN
        atm=atms
        ityp=ityps
        zv=zvs
        tau=taus
     ENDif
     !
     if (nats>nat) CALL errore ('chdens', 'wrong file order? ', 1)
     if (dfftp%nr1x/=nr1sxa.or.dfftp%nr2x/=nr2sxa) CALL &
          errore ('chdens', 'incompatible nr1x or nr2x', 1)
     if (dfftp%nr1/=nr1sa.or.dfftp%nr2/=nr2sa.or.dfftp%nr3/=nr3sa) CALL &
          errore ('chdens', 'incompatible nr1 or nr2 or nr3', 1)
     if (ibravs/=ibrav) CALL errore ('chdens', 'incompatible ibrav', 1)
     if (abs(gcutmsa-gcutm)>1.d-8.or.abs(duals-dual)>1.d-8.or.&
         abs(ecuts-ecutwfc)>1.d-8) &
          CALL errore ('chdens', 'incompatible gcutm or dual or ecut', 1)
     if (ibravs /= 0 ) THEN
        do i = 1, 6
           if (abs( celldm (i)-celldms (i) ) > 1.0d-7 ) &
              CALL errore ('chdens', 'incompatible celldm', 1)
        enddo
     ENDif
     !
     rhor (:) = rhor (:) + weight (ifile) * rhos (:)
  enddo
  deallocate (ityps)
  deallocate (taus)
  deallocate (rhos)
  !
  ! open output file, i.e., "fileout"
  !
  if (ionode) THEN
     if (fileout /= ' ') THEN
        ounit = 1
        OPEN (unit=ounit, file=fileout, form='formatted', status='unknown')
        write( stdout, '(/5x,"Writing data to be plotted to file ",a)') &
             trim(fileout)
     ELSE
        ounit = 6
     ENDif
  ENDif
  ! the isostm subroutine is called only when isostm_flag is true and the
  ! charge density is related to an STM image (5) or is read from a file 
  if ( (isostm_flag) .AND. ( (plot_num == -1) .OR. (plot_num == 5) ) ) THEN
     if ( .NOT. (iflag == 2))&
        CALL errore ('chdens', 'isostm should have iflag = 2', 1)
        CALL isostm_plot(rhor, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, &
             isovalue, heightmin, heightmax, direction)     
  END if
  !
  !    At this point we start the calculations, first we normalize the
  !    vectors defining the plotting region.
  !    If these vectors have 0 length, replace them with crystal axis
  !
  m1 = sqrt (e1 (1)**2 + e1 (2)**2 + e1 (3)**2)
  if (abs(m1) < 1.d-6) THEN
     e1 (:) = at(:,1)
     m1 = sqrt (e1 (1)**2 + e1 (2)**2 + e1 (3)**2)
  ENDif
  e1 (:) = e1 (:) / m1
  !
  m2 = sqrt (e2 (1)**2 + e2 (2)**2 + e2 (3)**2)
  if (abs(m2) < 1.d-6) THEN
     e2 (:) = at(:,2)
     m2 = sqrt (e2 (1)**2 + e2 (2)**2 + e2 (3)**2)
  ENDif
  e2 (:) = e2 (:) / m2
  !
  m3 = sqrt (e3 (1)**2 + e3 (2)**2 + e3 (3)**2)
  if (abs(m3) < 1.d-6) THEN
     e3 (:) = at(:,3)
     m3 = sqrt (e3 (1)**2 + e3 (2)**2 + e3 (3)**2)
  ENDif
  e3 (:) = e3 (:) / m3
  !
  ! are vectors defining the plotting region aligned along xyz ?
  !
  fast3d = ( e1(2) == 0.d0  .and.  e1(3) == 0.d0) .and. &
           ( e2(1) == 0.d0  .and.  e2(3) == 0.d0) .and. &
           ( e3(1) == 0.d0  .and.  e3(2) == 0.d0)
  !
  ! are crystal axis aligned along xyz ?
  !
  fast3d = fast3d .and. &
       ( at(2,1) == 0.d0  .and.  at(3,1) == 0.d0) .and. &
       ( at(1,2) == 0.d0  .and.  at(3,2) == 0.d0) .and. &
       ( at(1,3) == 0.d0  .and.  at(2,3) == 0.d0)

     fast3d = fast3d .and. (trim(interpolation) == 'fourier')

     if (output_format == 5.and.ionode) THEN

        CALL gmap_sym(nsym, s, ftau, gmapsym, eigv, invs)

!RGO Correlation grid (should be set on input
!        ngmpol   = 2451
!        nwsigma  = 2
!        iunsigma = 32
!Si exchange grid (should be set on input)
!        ngmpol   = 283
!        nwsigma  = 1
!        iunsigma = 32
        ngmpol   = 59
        nwsigma  = 7
        iunsigma = 32
!        CALL q_points()

     call kpoint_grid(nsym, time_reversal, .false., s, t_rev,& 
                      bg, 64, 0,0,0, 4,4,4, nqs, xq, wq )

      xq(:,1) = 0.0d0

      do iq = 1, nqs
         write(stdout, '(5x,i3, 3f14.9)') iq, xq(1,iq), xq(2,iq), xq(3,iq)
      end do


        !allocate (sigmar(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))

        allocate (sigma_g(ngmpol, ngmpol, nwsigma))

!FOR INPUT FILE
        filename = trim(prefix)//"."//"sigma1"
        !filename = trim(prefix)//"."//"coul1"
        !filename = trim(prefix)//"."//"sigma_ex1"

        tempfile = trim(tmp_dir) //"_gw0/"// trim(filename)
        lrcoul = 2*ngmpol*ngmpol*nwsigma
        unf_recl = DIRECT_IO_FACTOR * int(lrcoul, kind=kind(unf_recl))
        open ( iunsigma, file = trim(adjustl(tempfile)), iostat = ios, form ='unformatted', & 
        status = 'unknown', access = 'direct', recl = unf_recl)
        write(stdout, '("iostat ", i4)'), ios
       !!! x0 and e1 are in alat units !!!
        allocate (rhog   ( ngm ))
        allocate (sigmag ( ngm ))
        sigmag = dcmplx(0.0d0, 0.0d0)
        rhog   = dcmplx(0.0d0, 0.0d0)
        sigmar = dcmplx(0.0d0, 0.0d0)

!        nx = dfftp%nr1x
!        ny = dfftp%nr2x
!        nz = dfftp%nr3x

        nx = 4*dfftp%nr1x
        ny = 4*dfftp%nr2x
        nz = 4*dfftp%nr3x

        allocate (sigmar(nx*ny*nz))

        write(stdout, '(/5x, "ngmpol ", i4, " nwsigma", i4 )') ngmpol, nwsigma
        write(stdout, '(/5x, "nsym ", i4, " nks", i4, "omega, " f12.4 )')   nsym,   nks , omega
        write(stdout, '(/5x, "x0(1) ", f7.4, " x0(2) ", f7.4, " x0(3) ", f7.4)') x0(1), x0(2), x0(3) 
        write(stdout, '(/5x, "e1(1) ", f7.4, " e1(2) ", f7.4, " e1(3) ", f7.4)') e1(1), e1(2), e1(3) 
        write(stdout, '(/5x, "e2(1) ", f7.4, " e2(2) ", f7.4, " e2(3) ", f7.4)') e2(1), e2(2), e2(3) 
        write(stdout, '(/5x, "e3(1) ", f7.4, " e3(2) ", f7.4, " e3(3) ", f7.4)') e3(1), e3(2), e3(3) 

        do iq = 1, nks 
           sigma_q = dcmplx(0.0,0.0)
           !read( iunsigma, rec = iq, iostat = ios) sigma_g
           call davcio( sigma_g, lrcoul, iunsigma, iq, -1) 
           write(stdout, '(/5x, "ios ", i4)') ios
           if(iq.eq.1) sigma_g(1,1,:) = dcmplx(0.0d0,0.0d0)
           do ig = 1, ngmpol
              sigma_g(ig,ig,:) = 1.0 + sigma_g(ig,ig,:)
           enddo
           write(stdout, '(/5x, "Mag Sig ", f12.7)') sum(sigma_g(:,:,:))
           do isym = 1, nsym
           !"single_point FFT to get \Sigma_{x0}(\r'; omega)
           !CALL rotate(xk(1,iq), aq, s, nsym, isym)
              CALL rotate(xq(1,iq), aq, s, nsym, isym)
              do ig = 1, ngmpol
                 eigx0 = exp((0.d0, -1.d0)*2.d0*pi*(&
                              !x0(1)*4.0*(g(1,gmapsym(ig, isym)) + aq(1)) + &
                              !x0(2)*4.0*(g(2,gmapsym(ig,isym))  + aq(2)) + &
                              !x0(3)*4.0*(g(3,gmapsym(ig,isym))  + aq(3))))
                              x0(1)*(g(1,gmapsym(ig, isym)) + aq(1)) + &
                              x0(2)*(g(2,gmapsym(ig,isym))  + aq(2)) + &
                              x0(3)*(g(3,gmapsym(ig,isym))  + aq(3))))
           !reduce r co-ordinate
                 sigmag(1:ngmpol) = sigmag(1:ngmpol) + (1.0/omega)*wk(iq)*eigx0*sigma_g(gmapsym(ig,isym), 1:ngmpol, 1)
              enddo
              CALL plot_3d_sig (celldm (1), at, nat, tau, atm, ityp, ngm, g, sigmag,&
                   nx, ny, nz, m1, m2, m3, aq, e1, e2, e3, output_format, &
                   ounit, rhotot, gmapsym(1,1), isym, sigmar)
           enddo
        enddo
        sigmar = (1.0/float(nsym))*sigmar
        CALL xsf_struct (alat, at, nat, tau, atm, ityp, ounit)
        CALL xsf_fast_datagrid_3d &
             (sigmar, dfftp%nr1, dfftp%nr2, dfftp%nr3, nx,&
              ny, nz, at, alat, ounit)
             !(sigmar, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x,&
             ! dfftp%nr2x, dfftp%nr3x, at, alat, ounit)
     ENDif
  !
  write(stdout, '(5x,"Plot Type: ",a,"   Output format: ",a)') &
       plotname(iflag), formatname(output_format)
  !
  if (allocated(rhog)) deallocate(rhog)
  deallocate(rhor)
  deallocate(tau)
  deallocate(ityp)
END SUBROUTINE sigmadens

SUBROUTINE plot_3d_sig (alat, at, nat, tau, atm, ityp, ngm, g, rhog, &
     nx, ny, nz, m1, m2, m3, xk, e1, e2, e3, output_format, ounit, &
     rhotot, gmapsym, isym, carica)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE constants, ONLY:  pi
  USE io_global, ONLY : stdout, ionode
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_sum
  USE symm_base,  ONLY : nsym, s, time_reversal, t_rev, ftau, invs, nrot
  IMPLICIT NONE
  integer :: nat, ityp (nat), ngm, nx, ny, nz, output_format, ounit
  ! number of atoms
  ! type of atoms
  ! number of G vectors
  ! number of points along x, y, z
  ! output format
  ! output unit
  character(len=3) :: atm(*)
  real(DP) :: alat,  tau(3,nat), at(3,3), g(3,ngm), xk(3), &
              e1(3), e2(3), e3(3), m1, m2, m3
  ! lattice parameter
  ! atomic positions
  ! lattice vectors
  ! G-vectors
  ! origin
  ! vectors e1,e2,e3 defining the parallelepiped
  ! moduli of e1,e2,e3
  complex(DP) :: rhog (ngm)
  ! rho or polarization in G space
  integer :: i, j, k, ig
  ! min, max value of the charge, total charge, total absolute charge
  ! steps along e1, e2, e3
  complex(DP), allocatable :: eigx (:), eigy (:), eigz (:)
  integer  :: gmapsym(ngm, nrot), isym
  real(DP) :: rhomin, rhomax, rhotot, rhoabs, deltax, deltay, deltaz
  real(DP) :: omega
!HL
! real(DP), allocatable :: carica (:,:,:)
  real(DP)  :: carica (nx,ny,nz)

! HL:
! allocate (carica( nx , ny , nz))
  allocate (eigx(  nx))
  allocate (eigy(  ny))
  allocate (eigz(  nz))

!  deltax = m1 / nx
!  deltay = m2 / ny
!  deltaz = m3 / nz

  deltax = 4.0*m1 / nx
  deltay = 4.0*m2 / ny
  deltaz = 4.0*m3 / nz

 !HL want to accumulate sum
 !carica = 0.d0
  do ig = 1, ngm
     !
     ! eigx=exp(iG*e1+iGx0), eigy=exp(iG*e2), eigz=exp(iG*e3)
     ! These factors are calculated and stored in order to save CPU time
     !
     do i = 1, nx
        eigx (i) = exp( (0.d0,1.d0) * 2.d0 * pi * ( (i-1) * deltax * &
             (e1(1)*(g(1,gmapsym(ig,isym)) +  xk(1)) +  &
              e1(2)*(g(2,gmapsym(ig,isym)) + xk(2))  +  &
              e1(3)*(g(3,gmapsym(ig,isym)) + xk(3)))))
     enddo
     do j = 1, ny
        eigy (j) = exp( (0.d0,1.d0) * 2.d0 * pi * (j-1) * deltay * &
             (e2(1)*(g(1,gmapsym(ig,isym))+xk(1)) + e2(2)*(g(2,gmapsym(ig,isym))+xk(2)) + &
              e2(3)*(g(3,gmapsym(ig,isym))+xk(3))))
     enddo
     do k = 1, nz
        eigz (k) = exp((0.d0,1.d0)*2.d0*pi*(k-1)*deltaz* &
             (e3(1)*(g(1,gmapsym(ig,isym))+xk(1)) + e3(2)*(g(2,gmapsym(ig,isym))+xk(2)) + &
              e3(3)*(g(3,gmapsym(ig,isym))+xk(3))))
     enddo
     do k = 1, nz
        do j = 1, ny
           do i = 1, nx
              carica (i, j, k) = carica (i, j, k) + &
                    dble (rhog (gmapsym(ig,isym)) * eigz (k) * eigy (j) * eigx (i) )
           enddo
        enddo
     enddo
  enddo
  !
  !
  CALL mp_sum( carica, intra_bgrp_comm )
  !
  ! HL:
  ! Here we check the value of the resulting exchange correlation 
  ! potential.
  !
  CALL volume(alat,e1(1),e2(1),e3(1),omega)

  rhomin = max    ( minval (carica), 1.d-10 )
  rhomax = maxval (carica)
  rhotot = sum (carica(:,:,:))      * omega * deltax * deltay * deltaz
  rhoabs = sum (abs(carica(:,:,:))) * omega * deltax * deltay * deltaz

  write(stdout, '(/5x, "RUNNING THIS ROUTINE!")')
  write(stdout, '(/5x,"Min, Max, Total, Abs charge: ",2f10.6,2x, 2f10.4)')&
     rhomin, rhomax, rhotot, rhoabs
!
!  if (ionode) THEN
!     if (output_format == 4) THEN
!       "gOpenMol" file
!        CALL write_openmol_file (alat, at, nat, tau, atm, ityp, x0, &
!             m1, m2, m3, nx, ny, nz, rhomax, carica, ounit)
!     ELSE
        ! user has calculated for very long, be nice and write some output even
        ! if the output_format is wrong; use XSF format as default
        !
        ! XCRYSDEN's XSF format
        !
!        CALL xsf_struct      (alat, at, nat, tau, atm, ityp, ounit)
!        CALL xsf_datagrid_3d &
!             (carica, nx, ny, nz, m1, m2, m3, x0, e1, e2, e3, alat, ounit)
!     ENDif
!  ENDif
!  deallocate (carica)
  deallocate (eigz)
  deallocate (eigy)
  deallocate (eigx)
  return
END SUBROUTINE plot_3d_sig
