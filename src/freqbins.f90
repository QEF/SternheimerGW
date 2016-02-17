!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! 
! Copyright (C) 2010 - 2016 
! Henry Lambert, Martin Schlipf, and Feliciano Giustino
!
! Sternheimer-GW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Sternheimer-GW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Sternheimer-GW. If not, see
! http://www.gnu.org/licenses/gpl.html .
!
!------------------------------------------------------------------------------ 
SUBROUTINE freqbins()
  !
  ! generate frequency bins
  ! ----------------------------------------------------------------
  ! Here I assume Sigma is needed for w0 between wsigmamin and wsigmamax
  ! The convolution requires W for positive frequencies w up to wcoulmax
  ! (even function - cf Shishkin and Kress) and the GF spanning w0+-w.
  ! Therefore the freq. range of GF is
  ! from (wsigmamin-wcoulmax) to (wsigmamax+wcoulmax)
  ! the freq. dependence of the GF is inexpensive, so we use the same spacing
  ! NB: I assume wcoulmax>0, wsigmamin=<0, wsigmamax>0 and zero of energy at the Fermi level
  ! TODO HL: should set two frequency windows one fine grid for the range around the fermi level
  ! say  ef +/- 60 eV  down to the lowest pseudo state included! and a second 
  ! course window for everything outside this range. 
 
  USE freq_gw,    ONLY : nwcoul, nwgreen, nwalloc, nwsigma, wtmp, wcoul,& 
                         wgreen, wsigma, wsigmamin, wsigmamax,&
                         deltaw, wcoulmax, ind_w0mw, ind_w0pw, wgreenmin,&
                         wgreenmax, fiu, nfs, greenzero, w0pmw, wgtcoul
  USE io_global,  ONLY :  stdout, ionode, ionode_id
  USE kinds,      ONLY : DP
  USE constants,  ONLY : RYTOEV, pi
  USE control_gw, ONLY : eta, godbyneeds, padecont, freq_gl, do_imag
  USE disp,       ONLY : num_k_pts, w_of_k_start, w_of_k_stop, xk_kpoints,& 
                         nq1, nq2, nq3
  USE cell_base,  ONLY : omega, tpiba2, at, bg, tpiba, alat
           

  IMPLICIT NONE 

  REAL(DP)  :: zero, w0mw, w0pw, rcut
  REAL(DP)  :: x1, x2       !ranges
  REAL(DP), ALLOCATABLE  :: x(:), w(:)   !abcissa and weights
  INTEGER  :: n            !number of points
  INTEGER  :: iw, iw0, iwp, iw0mw, iw0pw, i
  LOGICAL  :: foundp, foundm

  zero = 0.0d0
  IF(.not.freq_gl) THEN
   !wgreenmin = wsigmamin-wcoulmax
   !wgreenmax = wsigmamax+wcoulmax
!for change of variables in Green's function
   wgreenmin = -wcoulmax
   wgreenmax = wcoulmax

   nwalloc = 1 + ceiling( (wgreenmax-wgreenmin) / deltaw )
   allocate(wtmp(nwalloc), wcoul(nwalloc), wgreen(nwalloc), wsigma(nwalloc) )
   wcoul = zero
   wgreen = zero
   wsigma = zero
   do iw = 1, nwalloc
      wtmp(iw) = wgreenmin + (wgreenmax-wgreenmin)/float(nwalloc-1)*float(iw-1)
   enddo
   nwgreen = 0
   nwcoul = 0
   nwsigma = 0
   do iw = 1, nwalloc
    if ( ( wtmp(iw) .ge. zero ) .and. ( wtmp(iw) .le. wcoulmax) ) then
      nwcoul = nwcoul + 1
      wcoul(nwcoul) = wtmp(iw)
    endif
    if ( ( wtmp(iw) .ge. wsigmamin ) .and. ( wtmp(iw) .le. wsigmamax) ) then
      nwsigma = nwsigma + 1
      wsigma(nwsigma) = wtmp(iw)
    endif
   enddo

   nwgreen = 2*nwcoul
   allocate ( w0pmw (nwsigma, 2*nwcoul) )
   DO iw0 = 1, 1
     DO iw = 1, nwcoul
        w0pmw(iw0, iw)        = wsigma(iw0) + wcoul(iw)
        w0pmw(iw0, iw+nwcoul) = wsigma(iw0) - wcoul(iw)
     ENDDO
   ENDDO
  ELSE
! We generate Sigma on a uniform grid:
   nwgreen = 2*nwcoul
   nwsigma = 1 + ceiling( (wsigmamax-wsigmamin) / deltaw )
   allocate(wsigma(nwsigma))
   allocate(wcoul(nwcoul), wgtcoul(nwcoul))
   wcoul   = zero
   wgtcoul = zero
   wsigma  = zero
   do iw = 1, nwsigma
      wsigma(iw) = wsigmamin + (wsigmamax-wsigmamin)/float(nwsigma-1)*float(iw-1)
   enddo
! We generate W(iw) on a gauss legendre grid:
   wcoul   = zero
   wgtcoul = zero
   CALL gauleg_grid(0.0d0, wcoulmax, wcoul, wgtcoul, nwcoul)
! Correpspondence arrays for convolution arguments for green's function.
   allocate ( w0pmw (nwsigma, 2*nwcoul) )
   DO iw0 = 1, nwsigma
     DO iw = 1, nwcoul
        w0pmw(iw0, iw)        = wsigma(iw0) + wcoul(iw)
        w0pmw(iw0, iw+nwcoul)        = wsigma(iw0) - wcoul(iw)
     ENDDO
   ENDDO
   WRITE(stdout, '(7x, "Gauss-Legendre grid: ")')
   DO i = 1, nwcoul
      WRITE(stdout,'(8x, i4, 4x, 2f12.3)') i, wgtcoul(i), wcoul(i)
   ENDDO
  ENDIF
! Print out Frequencies on Imaginary Axis for reference.
  WRITE(stdout, '(//5x,"Frequency Grids (eV):")')
  WRITE(stdout, '(/5x, "wsigmamin, wsigmamax, deltaw")')
  WRITE(stdout, '(5x, 3f10.4 )') wsigmamin, wsigmamax, deltaw 
  WRITE(stdout, '(/5x, "wcoulmax:", 1f10.4, " eV")'), wcoulmax
  WRITE(stdout, '(5x, "nwgreen:", i5)'), nwgreen
  WRITE(stdout, '(5x, "nwcoul:", i5)'), nwcoul
  WRITE(stdout,'(//5x, "Dynamic Screening Model:")')
  IF(godbyneeds) then
      WRITE(stdout, '(/6x, "Godby Needs Plasmon-Pole")')
  else if (padecont) then
      WRITE(stdout, '(/6x, "Analytic Continuation")')
  else if (.not.padecont.and..not.godbyneeds) then
      WRITE(stdout, '(/6x, "No screening model chosen!")')
  ENDIF
  WRITE(stdout, '(/7x, "Imag. Frequencies: ")')
  DO i = 1, nfs
       WRITE(stdout,'(8x, i4, 4x, 2f9.4)')i, fiu(i)*RYTOEV
  ENDDO
  WRITE(stdout, '(/5x, "Broadening: ", 1f10.4)'), eta

  rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))
  WRITE(stdout, '(/5x, "Spherical Cutoff: ", 1f10.4)'), rcut

  WRITE(stdout, '(/7x, "K-points: ", i4)') num_k_pts
  DO i = 1, num_k_pts
       WRITE(stdout,'(8x, i4, 4x, 3f9.4)')i, xk_kpoints(:, i)
  ENDDO

END SUBROUTINE freqbins
