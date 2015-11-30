PROGRAM diel_vis
  USE environment,ONLY : environment_start
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE mp_global,  ONLY : nproc_pool
  USE mp,         ONLY : mp_bcast
  USE parameters, ONLY : ntypx
  USE constants,  ONLY :  pi, fpi, tpi
  USE cell_base
  USE ions_base,         ONLY : nat, ityp, atm, ntyp => nsp, tau, zv
  USE lsda_mod,          ONLY: nspin
  USE gvect
  USE gsmooth
  USE wavefunctions_module,  ONLY: psic
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc, find_free_unit, outdir
  USE fft_base,              ONLY : grid_scatter
  USE printout_base,         ONLY : title
  USE control_flags,         ONLY : gamma_only
  USE scf,                   ONLY : rho, vltot, v
  USE klist,                 ONLY : nks, xk, wk
  USE wvfct,                ONLY : et, nbnd, npwx, npw, igk, g2kin
  USE wavefunctions_module, ONLY : evc

IMPLICIT NONE

integer, parameter :: DIRECT_IO_FACTOR = 8
real(DP), parameter :: ryd2ev = 13.6058
complex(DP), allocatable :: a(:)
REAL(DP), allocatable    ::  wim(:), rhor(:)
REAL(DP)                 :: ar, ai
INTEGER :: nwim, iq, ios, iw, nqs, unitcell
INTEGER :: iuncoul, iunsex, iunsigma
INTEGER :: ig, igp, iwim, recl
integer :: i, j
INTEGER*8 :: unf_recl
REAL (DP) :: eta, w_ryd, RYTOEV
integer :: glist(3), gplist(3), ngcount, ngpcount
integer :: ounit, style, output_format
complex(DP), allocatable :: invepsz(:)
complex(DP) :: arg, phase
integer :: nx, ny, ngmpol, ipol, idum, na
integer :: ngmsex, nwsigma
complex(DP), allocatable :: rhog (:), fiu(:)
integer  :: plot_num
CHARACTER(len=256) :: filplot
COMPLEX(DP), ALLOCATABLE :: scrcoul_g(:,:,:), sigma(:,:,:), sigma_ex(:,:)
COMPLEX(DP), ALLOCATABLE :: fxyz(:,:,:)
INTEGER                  :: x, y, z
LOGICAL                  :: sigorcoul, corr
LOGICAL                  :: plot_re
LOGICAL                  :: plot_rho
real(DP) :: rhomin, rhomax, rhoim, deltax, deltay
REAL (DP) :: xi,yi,zi
  ! minimum value of the charge
  ! maximum value of the charge
  ! integrated imaginary charge
  ! steps along e1
  ! steps along e2
complex(DP), allocatable :: eigx (:), eigy (:), carica(:,:), carical(:)
complex(DP), allocatable :: evctr(:)

!e1,e2, e3 are vectors that
real(DP) :: e1(3), e2(3), e3(3), x0 (3), R(3), radius, m1, m2, m3, e(3)

!ounit just standard out.
  ounit = 400
  CALL environment_start ( 'POST-PROC' )

  CALL input_from_file ( )

  prefix = 'si'
  outdir = './tmp'

  !    This subroutine reads the data for the output file produced by pw.x
  !    extracts and calculates the desired quantity (rho, V, ...)
  !    writes it to a file for further processing or plotting
  call extract (filplot, plot_num) 
!!!!!!!!! SHOULD FORMULATE INPUT FILE !!!!!!!!!!!!!!!!!
!ngmpol: number of G-vectors used to describe dielectric matrix.
!nwim: number of imaginary frequencies
!for Si:
  ngmpol    = 137
  ngmsex    = 411
  nwsigma   = 3
!sigorcoul
  sigorcoul = .true.
  plot_rho  = .false.
!output format
  output_format = 0
!vectors for 
 !!! x0, e1, e2 are in alat units !!!
     e1(:) = (/1.0, 1.0, 0.0/) !!!
     e2(:) = (/0.0, 0.0, 1.0/) !!!
     e3(:) = (/0.0, 0.0, 0.0/)
     x0(:) = (/0.0, 0.0, 0.0/)
!number of points in the plane
!     rho(i,j) = rho( x0 + e1 * (i-1)/(nx-1)
!              + e2 * (j-1)/(ny-1) ), i=1, nx ; j=1, ny
     nx    = 80
     ny    = 40

     if (e1(1)**2 + e1(2)**2 + e1(3)**2 <  1d-6 .or. &
         e2(1)**2 + e2(2)**2 + e2(3)**2 <  1d-6)     &
         call errore ('chdens', 'missing e1/e2 vectors', 1)

     if (abs(e1(1)*e2(1) + e1(2)*e2(2) + e1(3)*e2(3)) > 1d-6) &
         call errore ('chdens', 'e1 and e2 are not orthogonal', 1)

     if (nx <= 0 .or. ny <= 0 )   call errore ('chdens', 'wrong nx/ny', 2)

  !
  !    At this point we start the calculations, first we normalize the 
  !    vectors defining the plotting region. 
  !    If these vectors have 0 length, replace them with crystal axis
  !

  !    At this point we start the calculations, first we normalize the 
  !    vectors defining the plotting region. 
  !    If these vectors have 0 length, replace them with crystal axis

  m1 = sqrt (e1 (1)**2 + e1 (2)**2 + e1 (3)**2)
  if (abs(m1) < 1.d-6) then
     e1 (:) = at(:,1)
     m1 = sqrt (e1 (1)**2 + e1 (2)**2 + e1 (3)**2)
  end if
  e1 (:) = e1 (:) / m1
  !
  m2 = sqrt (e2 (1)**2 + e2 (2)**2 + e2 (3)**2)
  if (abs(m2) < 1.d-6) then
     e2 (:) = at(:,2)
     m2 = sqrt (e2 (1)**2 + e2 (2)**2 + e2 (3)**2)
  end if
  e2 (:) = e2 (:) / m2
  !
  m3 = sqrt (e3 (1)**2 + e3 (2)**2 + e3 (3)**2)
  if (abs(m3) < 1.d-6) then
     e3 (:) = at(:,3)
     m3 = sqrt (e3 (1)**2 + e3 (2)**2 + e3 (3)**2)
  end if
  e3 (:) = e3 (:) / m3

!READ in dielectric matrix...
     open (45, file = "./imfreq.dat", form = 'formatted', status = 'unknown')
           read (45, *)
           read (45, *) nwim

     allocate(fiu(nwim))

     if (sigorcoul) then 
         ALLOCATE (sigma_ex (ngmsex, ngmsex))
         ALLOCATE (sigma(ngmpol,ngmpol,nwsigma))
     else
         ALLOCATE ( scrcoul_g (ngmpol, ngmpol, nwim) )
     endif

     do iw = 1, nwim
        read (45, *) ar, ai
        fiu(iw) = dcmplx(ar, ai)
        write(6,*) fiu(iw)
     enddo
     close(45)

     if(sigorcoul) then
        iunsex = 33
        recl = 2 * ngmsex * ngmsex
        unf_recl = 8 * int(recl, kind=kind(unf_recl))
        open ( iunsex, file = "./_gw0si.sigma_ex1", iostat = ios, form = 'unformatted', &
               status = 'unknown', access = 'direct', recl = unf_recl)

        iunsigma = 32
        recl = 2 * ngmpol * ngmpol * nwsigma
        unf_recl = 8 * int(recl, kind=kind(unf_recl))
        open ( iunsigma, file = "./_gw0si.sigma1", iostat = ios, form = 'unformatted', &
               status = 'unknown', access = 'direct', recl = unf_recl)

     else
       iuncoul = 28
       recl = 2 * ngmpol * ngmpol * nwim
       unf_recl = 8 * int(recl, kind=kind(unf_recl))
       open ( iuncoul, file = "./_gw0si.coul1", iostat = ios, form = 'unformatted', &
              status = 'unknown', access = 'direct', recl = unf_recl)
     endif

     write(6,*) ios

     iq = 1
     iw = 1
!Number of K-points:

     if(sigorcoul) then
     !   READ ( iunsigma, rec = iq, iostat = ios) sigma
     !   READ ( iunsex,   rec = iq, iostat = ios) sigma_ex
     else
        READ ( iuncoul,  rec = iq, iostat = ios) scrcoul_g
     endif

     if(.not.sigorcoul) then
    !add one back since we store: \epsilon^{-1} - 1.
         do ig = 1, ngmpol
            scrcoul_g(ig,ig,iw) = scrcoul_g(ig,ig,iw) + 1.0d0
         enddo
     endif

ALLOCATE  (rhog( ngm))
ALLOCATE  (rhor(nrx1*nrx2*nrx3))

rhog(:) = (0.0d0,0.0d0)

style=7
if(style.eq.1) then
! response in cubefile format to a point test charge at a particular point in real space.
       do igp   = 1, ngmpol
          do ig = 1, ngmpol
            !arg = tpi*dcmplx(0.0d0, 1.0d0)*(e3(1)*(g(1,ig) + e(2)*g(2,ig) + e(3)*g(3,ig))
             arg = -tpi*(e3(1)*g(1,ig) + e3(2)*g(2,ig) + e3(3)*g(3,ig))
            !rhog(igp) = rhog(igp) + scrcoul_g(ig,igp,iw)*exp(dcmplx(0.0d0, 1.0d0)*arg)
            !to plot sigma_{\r}(\r').
             rhog(igp) = rhog(igp) + sigma(ig,igp,6)*exp(dcmplx(0.0d0, 1.0d0)*arg)
          enddo
       enddo

       rhog(:) = (1.0d0/sqrt(omega))*rhog(:)

       ALLOCATE  (invepsz(nrx1*nrx2*nrx3))

       invepsz(nl(:)) = rhog(:)

       call cft3 (invepsz, nr1, nr2, nr3, nrx1, nrx2, nrx3, +1)

       rhor(:) = (1.0d0/sqrt(omega))*abs(real(invepsz(:)))

       call write_cubefile (alat, at, bg, nat, tau, atm, ityp, rhor, &
             nr1, nr2, nr3, nrx1, nrx2, nrx3, ounit)

else if (style.eq.2) then
!GNUPLOT FORMAT OUTPUT FOR function_{\r}(\r').
     do igp   = 1, ngmpol
        do ig = 1, ngmpol
           arg = -tpi*(e3(1)*g(1,ig) + e3(2)*g(2,ig) + e3(3)*g(3,ig))
          !rhog(igp) = rhog(igp) + scrcoul_g(ig,igp,iw)*exp(dcmplx(0.0d0, 1.0d0)*arg)
          !to plot sigma_{\r}(\r').
           rhog(igp) = rhog(igp) + sigma(ig,igp, 11)*exp(dcmplx(0.0d0, 1.0d0)*arg)
        enddo
     enddo
     rhog(:) = (1.0d0/sqrt(omega))*rhog(:)
     ALLOCATE  (invepsz(nrx1*nrx2*nrx3))
     invepsz(nl(:)) = rhog(:)
     call cft3 (invepsz, nr1, nr2, nr3, nrx1, nrx2, nrx3, +1)
     rhor(:) = (1.0d0/sqrt(omega))*(real(invepsz(:)))

     allocate(fxyz(nr1,nr2,nr3))
  
     DO x = 1, nr1
          DO y = 1, nr2
             DO z = 1, nr3
                fxyz(x,y,z) = fxyz(x,y,z) + rhor(x + (y-1)*nr2 + (z-1)*nr3*nr2)
             ENDDO
         ENDDO
     ENDDO
!writing potential in xyz format should be middle of cell:
     write(6,*) int(nr1/2)
     !z = int(nr3/4)
     z = 10
     do x = int(nr1/2)+1, nr1
        do y = int(nr2/2)+1, nr2
              write (450, '(2i4, 1f12.7)') x-nr1, y-nr2, real(fxyz(x,y,z))
        enddo
        do y = 1, int(nr2/2)
              write (450, '(2i4, 1f12.7)') x-nr1, y, real(fxyz(x,y,z))
        enddo
           write (450,*)
     enddo

     do x = 1, int(nr1/2)
        do y = int(nr2/2)+1, nr2
              write (450, '(2i4, 1f12.7)') x, y-nr2, real(fxyz(x,y,z))
        enddo
        do y = 1, int(nr2/2)
              write (450, '(2i4, 1f12.7)') x, y, real(fxyz(x,y,z))
        enddo
            write (450,*)
     enddo

else if (style.eq.3) then
     !\inveps_{\G}(G').
     !ig = 9
     !  do igp = 1, ngmpol
     !      rhog(igp) = rhog(igp) + scrcoul_g(ig,igp,iw)
     !  enddo
       ALLOCATE  (invepsz(nrx1*nrx2*nrx3))
       invepsz(nl(:)) = rhog(:)
       call cft3 (invepsz, nr1, nr2, nr3, nrx1, nrx2, nrx3, +1)
       rhor = (1.0d0/sqrt(omega))*real(invepsz(:))
       call write_cubefile (alat, at, bg, nat, tau, atm, ityp, rhor, &
             nr1, nr2, nr3, nrx1, nrx2, nrx3, ounit)

else if (style.eq.4) then
!planar fxns.
! If one wanted to plot charge density they would do this.
     !   READ ( iunsigma, rec = iq, iostat = ios) sigma
     !   READ ( iunsex,   rec = iq, iostat = ios) sigma_ex
  if(plot_rho) then
       ALLOCATE  (invepsz(nrx1*nrx2*nrx3))
       invepsz = rho%of_r(:,1)
       call cft3 (invepsz, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1)
       rhog (:) = invepsz(nl(:))
  else
! If one wanted to plot Sigma they would do...
       do igp   = 1, ngmpol
          do ig = 1, ngmpol
             arg = -tpi*(x0(1)*g(1,ig) + x0(2)*g(2,ig) + x0(3)*g(3,ig))
             rhog(igp) = rhog(igp) + sigma(ig,igp,1)*exp(dcmplx(0.0d0, 1.0d0)*arg)
          enddo
       enddo
       rhog(:) = (1.0d0/sqrt(omega))*rhog(:)
!or sigma_x
!       do igp   = 1, ngmsex
!          do ig = 1, ngmsex
!             arg = -tpi*(x0(1)*g(1,ig) + x0(2)*g(2,ig) + x0(3)*g(3,ig))
!             rhog(igp) = rhog(igp) + sigma_ex(ig,igp)*exp(dcmplx(0.0d0, 1.0d0)*arg)
!          enddo
!       enddo
!       rhog(:) = (1.0d0/sqrt(omega))*rhog(:)
  endif!plotrho
       call plot_2d (nx, ny, m1, m2, x0, e1, e2, ngm, g, rhog, alat, &
                     at, nat, tau, atm, ityp, output_format, ounit)
else if (style.eq.5) then
     allocate (eigx(nx))    
     allocate (eigy(ny))    
     allocate (carica(nx, ny))    
     write(6, *) nks
!x0 determines the origin in crystal co-ordinates.
     corr  = .true.
!    x0(:) = (/0.0, 0.0, 0.0/)
!    x0(:) = (/0.25, 0.25, 0.25/)
!    x0(:) = (/0.125, 0.125, 0.125/)
!    e3(:) = (/0.125, 0.125, 0.125/)
     x0(:) = (/0.0, 0.0, 0.0/)
     e3(:) = (/0.0, 0.0, 0.0/)
!m3 = sqrt (e3 (1)**2 + e3 (2)**2 + e3 (3)**2)
!e3 (:) = e3 (:) / m3
     carica(:,:) = dcmplx(0.0d0, 0.0d0)
     rhog(:)     = dcmplx(0.0d0, 0.0d0)
     do iq = 1, 10
        rhog(:)  = dcmplx(0.0d0, 0.0d0)
        sigma    = dcmplx(0.0d0, 0.0d0)
        sigma_ex = dcmplx(0.0d0, 0.0d0)
        READ ( iunsigma, rec = iq, iostat = ios) sigma
        READ ( iunsex,   rec = iq, iostat = ios) sigma_ex
        if(ios.ne.0) write(6, '("i/o error")')
        if(ios.ne.0) stop
        write(6,*) xk(:,iq)
!Real space lattice displacement in cartesian 
!       R(:) = (/1.0*unitcell, 1.0*unitcell, 1.0*unitcell/)
        R(:) = (/0.0, 0.0, 0.0/)
        write(6,*) phase
           do igp = 1, ngmpol
              do ig = 1, ngmpol
                 arg       = -tpi*(e3(1)*(g(1,ig) + xk(1,iq)) + e3(2)*(g(2,ig) + xk(2,iq)) + e3(3)*(g(3,ig) + xk(3,iq)))
                 phase     = exp(dcmplx(0.0d0,1.0d0)*(-tpi*(x0(1)*(g(1,igp)-g(1,ig))+x0(2)*(g(2,igp)-g(2,ig))+x0(3)*(g(3,igp)-g(3,ig)))))
                 rhog(igp) = rhog(igp) + wk(iq)*sigma(ig,igp,1)*exp(dcmplx(0.0d0,1.0d0)*arg)*phase
              enddo
           enddo
           do igp = 1, ngmsex
              do ig = 1, ngmsex
                 arg       = -tpi*(e3(1)*(g(1,ig) + xk(1,iq)) + e3(2)*(g(2,ig) + xk(2,iq)) + e3(3)*(g(3,ig) + xk(3,iq)))
                 phase     = exp(dcmplx(0.0d0,1.0d0)*(-tpi*(x0(1)*(g(1,igp)-g(1,ig))+x0(2)*(g(2,igp)-g(2,ig))+x0(3)*(g(3,igp)-g(3,ig)))))
                 rhog(igp) = rhog(igp) + wk(iq)*sigma_ex(ig,igp)*exp(dcmplx(0.0d0,1.0d0)*arg)*phase
              enddo
           enddo
       !if(iq.eq.1) rhog(1) =  1.0d0 
        rhog(:) = (1.0d0/omega)*rhog(:) 
        deltax  = m1 / (nx - 1)
        deltay  = m2 / (ny - 1)
        do ig = 1, ngm
           ! eigx=exp(iG*e1+iGx0), eigy=(iG*e2)
           ! These factors are calculated and stored in order to save CPU time
           do i = 1, nx
           !phase in two point fxn is already included.
              eigx (i) = exp ( (0.d0, 1.d0) * 2.d0 * pi * ( (i - 21) * deltax * &
              (e1(1)*(g(1,ig)+xk(1,iq)) + e1(2)*(g(2,ig)+xk(2,iq)) + e1(3)*(g(3,ig)+xk(3,iq))  ) ) )
           !(e1(1) * (g(1,ig)+xk(1,iq)) + e1(2) * (g(2,ig) + xk(2,iq)) + e1(3)  * (g(3,ig)+xk(3,iq)) ) + &
           !(x0 (1) *(g(1,ig))+ x0 (2) * (g(2,ig)) + x0 (3) * (g(3,ig)) ) ) )
           enddo
           do j = 1, ny
              eigy (j) = exp ( (0.d0, 1.d0) * 2.d0 * pi * (j - 21) * deltay * &
                   (e2(1) * (g(1,ig)+xk(1,iq)) + e2(2) * (g(2,ig)+xk(2,iq))+ e2(3) * (g(3,ig)+xk(3,iq))) )
           enddo
           do j = 1, ny
              do i = 1, nx
                 carica (i, j) = carica (i, j) + rhog (ig) * eigx (i) * eigy (j)
              enddo
           enddo
        enddo
     enddo
     do i = 1, nx
           !write (500, '(e25.14)') (DBLE(carica(i,j)), j = 1, ny)
           write (500+unitcell, '(e25.14)') (real(carica(i,j)), j = 1, ny)
           write (500+unitcell, *)
     enddo
     carica(:,:) = (0.0d0, 0.0d0)
!   enddo
!if(corr) then
!else
!endif
else if (style.eq.6) then
!Plot 1d
     nx = 50
     allocate (eigx(  nx))    
     allocate (eigy(  ny))    
     allocate (carical(nx))    
     write(6,*) nks
     corr  = .true.
     x0(:) = (/0.0, 0.0, 0.0/)
     e3(:) = (/0.0, 0.0, 0.0/)
     e(:) = (/1.0, 1.0, 1.0/)
     carica(:,:) = dcmplx(0.0d0, 0.0d0)
     rhog(:)     = dcmplx(0.0d0, 0.0d0)
     carical(:) = (0.d0,0.d0)
     m1    = sqrt (e (1)**2 + e (2)**2 + e (3)**2)
     e (:) = e (:) / m1

     do iq = 1, 10
        rhog(:)  = dcmplx(0.0d0, 0.0d0)
        sigma    = dcmplx(0.0d0, 0.0d0)
        sigma_ex = dcmplx(0.0d0, 0.0d0)
        READ ( iunsigma, rec = iq, iostat = ios) sigma
        READ ( iunsex,   rec = iq, iostat = ios) sigma_ex
        if(ios.ne.0) write(6, '("i/o error")')
        if(ios.ne.0) stop
        write(6,*) xk(:,iq)
!Real space lattice displacement in cartesian 
!       R(:) = (/1.0*unitcell, 1.0*unitcell, 1.0*unitcell/)
        R(:) = (/0.0, 0.0, 0.0/)
        write(6,*) phase
           do igp = 1, ngmpol
              do ig = 1, ngmpol
                 arg       = -tpi*(e3(1)*(g(1,ig) + xk(1,iq)) + e3(2)*(g(2,ig) + xk(2,iq)) + e3(3)*(g(3,ig) + xk(3,iq)))
                 phase     = exp(dcmplx(0.0d0,1.0d0)*(-tpi*(x0(1)*(g(1,igp)-g(1,ig))+x0(2)*(g(2,igp)-g(2,ig))+x0(3)*(g(3,igp)-g(3,ig)))))
                 rhog(igp) = rhog(igp) + wk(iq)*sigma(ig,igp,1)*exp(dcmplx(0.0d0,1.0d0)*arg)*phase
              enddo
           enddo
           do igp = 1, ngmsex
              do ig = 1, ngmsex
                 arg       = -tpi*(e3(1)*(g(1,ig) + xk(1,iq)) + e3(2)*(g(2,ig) + xk(2,iq)) + e3(3)*(g(3,ig) + xk(3,iq)))
                 phase     = exp(dcmplx(0.0d0,1.0d0)*(-tpi*(x0(1)*(g(1,igp)-g(1,ig))+x0(2)*(g(2,igp)-g(2,ig))+x0(3)*(g(3,igp)-g(3,ig)))))
                 rhog(igp) = rhog(igp) + wk(iq)*sigma_ex(ig,igp)*exp(dcmplx(0.0d0,1.0d0)*arg)*phase
              enddo
           enddo
       !if(iq.eq.1) rhog(1) =  1.0d0 
        rhog(:) = (1.0d0/omega)*rhog(:) 
        deltax = m1 / (nx - 1)
        do i = 1, nx
           xi = x0 (1) + (i - 1) * deltax * e (1)
           yi = x0 (2) + (i - 1) * deltax * e (2)
           zi = x0 (3) + (i - 1) * deltax * e (3)
        !
        !     for each point we compute the charge from the Fourier components
        !
          do ig = 1, ngm
        !
        !     NB: G are in 2pi/alat units, r are in alat units
        !
             arg = 2.d0 * pi * ( xi*(g(1,ig)+xk(1,iq)) + yi*(g(2,ig)+xk(2,iq)) + zi*(g(3,ig) + xk(3,iq)))
             carical(i) = carical(i) + rhog (ig) * exp(arg*(0.0d0,1.0d0))
          enddo
        enddo
     enddo !iq
     write (500, '(e25.14)') (real(carical(:)))
     write (500, *)
else if (style.eq.7) then
     x0(:) = (/0.0, 0.0, 0.0/)
     e3(:) = (/0.0, 0.0, 0.0/)
     e(:)  = (/1.0, 1.0, 1.0/)
     m1    = sqrt (e (1)**2 + e (2)**2 + e (3)**2)
     e (:) = e (:) / m1

!CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin) 
     write(6,*) nwordwfc, iunwfc
     CALL davcio (evc, nwordwfc, iunwfc, 1, - 1)
     allocate (evctr(npw))
!Plot real part of wave function along a line.
     evctr(:)= (0.0d0,0.0d0)
     do i = 2, 4
        evctr(:) = evctr(:) + (1.0/3.0)*evc(:,i) 
     enddo
  
     nx = 50
     allocate (carical(nx))    
     deltax = m1 / (nx - 1)
     do i = 1, nx
          xi = x0 (1) + (i - 1) * deltax * e (1)
          yi = x0 (2) + (i - 1) * deltax * e (2)
          zi = x0 (3) + (i - 1) * deltax * e (3)
         do ig = 1, npw
       !    NB: G are in 2pi/alat units, r are in alat units
            arg = 2.d0 * pi * ( xi*(g(1,ig)+xk(1,iq)) + yi*(g(2,ig)+xk(2,iq)) + zi*(g(3,ig) + xk(3,iq)))
            !carical(i) = carical(i) + evc(ig,6) * exp(arg*(0.0d0,1.0d0))
            carical(i) = carical(i) + evctr(ig) * exp(arg*(0.0d0,1.0d0))
         enddo
     enddo
     write(600,*)evc(:,3)
     write (500+unitcell, '(e25.14)') (real(carical(:)))
     write (500+unitcell, *)
endif!style

        write (ounit, '(i4)') nat
        write (ounit, '(3f8.4,i3)') ( (tau(ipol,na), ipol=1,3), 1, na=1,nat)
        write (ounit, '(f10.6)') celldm (1)
        write (ounit, '(3(3f12.6/))') at
END PROGRAM diel_vis


subroutine plot_2d (nx, ny, m1, m2, x0, e1, e2, ngm, g, rhog, alat, &
     at, nat, tau, atm, ityp, output_format, ounit)
  !-----------------------------------------------------------------------
  !
  USE kinds, only : DP
  use constants, only : pi
  use io_global, only : stdout, ionode
  USE mp_global,  ONLY : intra_pool_comm
  USE mp,         ONLY : mp_sum
  implicit none
  integer :: nx, ny, ngm, nat, ityp (nat), output_format, ounit
  ! number of points along x
  ! number of points along y
  ! number of G vectors
  ! number of atoms
  ! types of atoms
  ! output unit
  ! output format
  character(len=3) :: atm(*) ! atomic symbols
  real(DP) :: e1(3), e2(3), x0(3), m1, m2, g(3,ngm), alat, &
       tau(3,nat), at(3,3)
  ! vectors e1, e2 defining the plane
  ! origin
  ! modulus of e1
  ! modulus of e2
  ! G-vectors

  complex(DP) :: rhog (ngm)
  ! rho or polarization in G space
  integer :: i, j, ig

  real(DP) :: rhomin, rhomax, rhoim, deltax, deltay
  ! minimum value of the charge
  ! maximum value of the charge
  ! integrated imaginary charge
  ! steps along e1
  ! steps along e2
  complex(DP), allocatable :: eigx (:), eigy (:), carica(:,:)

  allocate (eigx(  nx))    
  allocate (eigy(  ny))    
  allocate (carica( nx , ny))    

  deltax = m1 / (nx - 1)
  deltay = m2 / (ny - 1)

  carica(:,:) = (0.d0,0.d0)
  do ig = 1, ngm
     !
     ! eigx=exp(iG*e1+iGx0), eigy=(iG*e2)
     ! These factors are calculated and stored in order to save CPU time
     !
     do i = 1, nx
        eigx (i) = exp ( (0.d0, 1.d0) * 2.d0 * pi * ( (i - 1) * deltax * &
             (e1(1) * g(1,ig) + e1(2) * g(2,ig) + e1(3) * g(3,ig) ) + &
             (x0 (1) * g(1,ig) + x0 (2) * g(2,ig) + x0 (3) * g(3,ig) ) ) )
     enddo
     do j = 1, ny
        eigy (j) = exp ( (0.d0, 1.d0) * 2.d0 * pi * (j - 1) * deltay * &
             (e2(1) * g(1,ig) + e2(2) * g(2,ig) + e2(3) * g(3,ig) ) )
     enddo
     do j = 1, ny
        do i = 1, nx
           carica (i, j) = carica (i, j) + rhog (ig) * eigx (i) * eigy (j)
        enddo
     enddo
  enddo
  !HL
  !call mp_sum( carica, intra_pool_comm ) 
  !
  !    Here we check the value of the resulting charge
  !
  rhomin =  1.0d10
  rhomax = -1.0d10

  rhoim = 0.d0
  do i = 1, nx
     do j = 1, ny
        rhomin = min (rhomin,  DBLE (carica (i, j) ) )
        rhomax = max (rhomax,  DBLE (carica (i, j) ) )
        rhoim = rhoim + abs (AIMAG (carica (i, j) ) )
     enddo

  enddo

  rhoim = rhoim / nx / ny
  write(stdout, '(5x,"Min, Max, imaginary charge: ",3f12.6)') &
                 rhomin, rhomax, rhoim

  !
  !     and we print the charge on output
  !
  if (ionode) then
     if (output_format == 0) then
        !
        !     gnuplot format
        !
        !         write(ounit,'(2i6)') nx,ny
        !HL
        !original
        !do i = 1, nx
        !   write (ounit, '(e25.14)') (  DBLE(carica(i,j)), j = 1, ny )
        !   write (ounit, * )
        !enddo
        !with folding along y
        do i = 1, nx
           write (ounit, '(e25.14)') (  DBLE(carica(i,j)), j = int(ny/2)+1, ny )
           write (ounit, '(e25.14)') (  DBLE(carica(i,j)), j = 1, (ny/2) )
           write (ounit, * )
        enddo
     elseif (output_format == 1) then
        !
        !     contour.x format
        !
        write (ounit, '(3i5,2e25.14)') nx, ny, 1, deltax, deltay
        write (ounit, '(4e25.14)') ( (  DBLE(carica(i,j)), j = 1, ny ), i = 1, nx )
     elseif (output_format == 2) then
        !
        !     plotrho format
        !
        write (ounit, '(2i4)') nx - 1, ny - 1
        write (ounit, '(8f8.4)') (deltax * (i - 1) , i = 1, nx)
        write (ounit, '(8f8.4)') (deltay * (j - 1) , j = 1, ny)
        write (ounit, '(6e12.4)') ( (  DBLE(carica(i,j)), i = 1, nx ), j = 1, ny )
        write (ounit, '(3f8.4)') x0
        write (ounit, '(3f8.4)') (m1 * e1 (i) , i = 1, 3)
        write (ounit, '(3f8.4)') (m2 * e2 (i) , i = 1, 3)
    ! elseif (output_format == 3) then
    !
    ! XCRYSDEN's XSF format
    !
    !    call xsf_struct (alat, at, nat, tau, atm, ityp, ounit)
    !    call xsf_datagrid_2d (carica, nx, ny, m1, m2, x0, e1, e2, alat, ounit)
     else
        call errore('plot_2d', 'wrong output_format', 1)
     endif
  endif

  deallocate (carica)
  deallocate (eigy)
  deallocate (eigx)
  return
end subroutine plot_2d
