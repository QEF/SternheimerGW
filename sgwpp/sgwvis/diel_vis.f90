PROGRAM diel_vis
! diel_vis reads in the dielectric matrix, stored in
! G-space and transforms it in to real space. Also with frequency
! dependence we can look at plasmon modes.
  USE environment,ONLY : environment_start
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE mp_global,  ONLY : nproc_pool
  USE mp,         ONLY : mp_bcast
  USE parameters, ONLY : ntypx
  USE constants,  ONLY :  pi, fpi, tpi
  USE cell_base
  USE ions_base,  ONLY : nat, ityp, atm, ntyp => nsp, tau, zv
  USE lsda_mod,   ONLY: nspin
  USE gvect
  USE gsmooth
  USE wavefunctions_module,  ONLY: psic
  USE io_files, ONLY: nd_nmbr, prefix, outdir
  USE fft_base,   ONLY: grid_scatter
  USE printout_base, ONLY: title
  USE control_flags, ONLY: gamma_only

IMPLICIT NONE

integer, parameter :: DIRECT_IO_FACTOR = 8
real(DP), parameter :: ryd2ev = 13.6058
complex(DP), allocatable :: z(:), a(:)
REAL(DP), allocatable    ::  wim(:), rhor(:)
REAL(DP)                 :: ar, ai
INTEGER :: nwim, iq, ios, iuncoul, iw
INTEGER :: ig, igp, iwim, recl
INTEGER*8 :: unf_recl
COMPLEX(DP), ALLOCATABLE :: scrcoul_g(:,:,:)
REAL (DP) :: eta, w_ryd, RYTOEV
integer :: glist(3), gplist(3), ngcount, ngpcount
integer :: ounit, style, output_format
complex(DP), allocatable :: invepsz(:)
complex(DP) :: arg
integer :: nx, ny, ngmpol, ipol, idum, na
complex(DP), allocatable :: rhog (:), fiu(:)
integer  :: plot_num
CHARACTER(len=256) :: filplot



LOGICAL :: plot_re

!e1,e2, e3 ar vector that
  real(DP) :: e1(3), e2(3), e3(3), x0 (3), radius, m1, m2, m3

!ounit just standard out.
  ounit = 6
  CALL environment_start ( 'POST-PROC' )

  CALL input_from_file ( )

  prefix = 'MoS2'
  outdir = './tmp'

  !    This subroutine reads the data for the output file produced by pw.x
  !    extracts and calculates the desired quantity (rho, V, ...)
  !    writes it to a file for further processing or plotting
  call extract (filplot, plot_num) 

!ngmpol: number of G-vectors used to describe dielectric matrix.
!nwim: number of imaginary frequencies
     iuncoul = 28
     ngmpol    = 125
!vectors for 
   !!! x0, e1, e2 are in alat units !!!
     e1(:) = (/1.0, 0.0, 0.0/)
     e2(:) = (/0.0, 1.0, 0.0/)
     e3(:) = (/0.5, 0.5, 0.5/)
   !!!origin of plane!!!
     x0(:) = (/0.0, 0.0, 0.0/)

!number of points in the plane
!     rho(i,j) = rho( x0 + e1 * (i-1)/(nx-1)
!              + e2 * (j-1)/(ny-1) ), i=1,nx ; j=1,ny

     nx    = 20
     ny    = 20

     if (e1(1)**2 + e1(2)**2 + e1(3)**2 <  1d-6 .or. &
         e2(1)**2 + e2(2)**2 + e2(3)**2 <  1d-6)     &
         call errore ('chdens', 'missing e1/e2 vectors', 1)

     if (abs(e1(1)*e2(1) + e1(2)*e2(2) + e1(3)*e2(3)) > 1d-6) &
         call errore ('chdens', 'e1 and e2 are not orthogonal', 1)

     if (nx <= 0 .or. ny <= 0 )   call errore ('chdens', 'wrong nx/ny', 2)

!READ in dielectric matrix...
!all this could apply to sigma as well...
     open (45, file = "./imfreq.dat", form = 'formatted', status = 'unknown')
           read (45, *)
           read (45, *) nwim

     allocate (fiu(nwim))
     allocate ( scrcoul_g (ngmpol, ngmpol, nwim) )

     do iw = 1, nwim
        read (45, *) ar, ai
        fiu(iw) = dcmplx(ar, ai)
        write(6,*) fiu(iw)
     enddo
     close(45)

     recl = 2 * ngmpol * ngmpol * nwim
     unf_recl = 8 * int(recl, kind=kind(unf_recl))
     open ( iuncoul, file = "./_gw0MoS2.coul1", iostat = ios, form = 'unformatted', &
            status = 'unknown', access = 'direct', recl = unf_recl)

      write(6,*) ios
      iq = 1
      READ ( iuncoul, rec = iq, iostat = ios) scrcoul_g

!Generate FFT grid unless this is already performed in extract?
!      call allocate_fft()
!      call ggen()

      iw = 1

      write(6,*) ngm
      write(6,*) g(:,1:3)
      write(6,*) scrcoul_g(1:3,1:3,1)

      do ig = 1, ngmpol
         scrcoul_g(ig,ig,iw) = scrcoul_g(ig,ig,iw) + 1.0d0
      enddo

allocate (rhog( ngm))
rhog(:) = (0.0d0,0.0d0)

style=2
if(style.eq.1) then
! response in cubefile format to a point test charge at a particular point in real space.
       do igp   = 1, ngmpol
          do ig = 1, ngmpol
            !arg = tpi*dcmplx(0.0d0, 1.0d0)*(e3(1)*(g(1,ig) + e(2)*g(2,ig) + e(3)*g(3,ig))
             arg = -tpi*(e3(1)*g(1,ig) + e3(2)*g(2,ig) + e3(3)*g(3,ig))
             rhog(igp) = rhog(igp) + scrcoul_g(ig,igp,iw)*exp(dcmplx(0.0d0, 1.0d0)*arg)
          enddo
       enddo

       rhog(:) = (1.0d0/sqrt(omega))*rhog(:)

       ALLOCATE  (invepsz(nrx1*nrx2*nrx3))

       invepsz(nl(:)) = rhog(:)

       call cft3 (invepsz, nr1, nr2, nr3, nrx1, nrx2, nrx3, +1)
       rhor(:) = (1.0d0/sqrt(omega))*real(invepsz(:))


       call write_cubefile (alat, at, bg, nat, tau, atm, ityp, rhor, &
             nr1, nr2, nr3, nrx1, nrx2, nrx3, ounit)

else if (style.eq.2) then
          
     !\inveps_{\G}(G').
       ig = 9
       do igp = 1, ngmpol
           rhog(igp) = rhog(igp) + scrcoul_g(ig,igp,iw)
       enddo

!add one back to the diagonal. 
       ALLOCATE  (invepsz(nrx1*nrx2*nrx3))
       ALLOCATE  (rhor(nrx1*nrx2*nrx3))

       invepsz(nl(:)) = rhog(:)

       call cft3 (invepsz, nr1, nr2, nr3, nrx1, nrx2, nrx3, +1)

       rhor = (1.0d0/sqrt(omega))*real(invepsz(:))

       call write_cubefile (alat, at, bg, nat, tau, atm, ityp, rhor, &
             nr1, nr2, nr3, nrx1, nrx2, nrx3, ounit)

else if (style.eq.3) then
! response in a particular plane to a plane wave perturbation. 
        call plot_2d (nx, ny, m1, m2, x0, e1, e2, ngm, g, rhog, alat, &
                      at, nat, tau, atm, ityp, output_format, ounit)
endif

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
  call mp_sum( carica, intra_pool_comm ) 
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
        do i = 1, nx
           write (ounit, '(e25.14)') (  DBLE(carica(i,j)), j = 1, ny )
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
