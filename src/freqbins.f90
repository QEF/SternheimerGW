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
  USE constants,  ONLY : RYTOEV
  USE control_gw, ONLY : eta, godbyneeds, padecont, freq_gl
           

  IMPLICIT NONE 

  LOGICAL  :: foundp, foundm
  REAL(DP) :: zero, w0mw, w0pw
  INTEGER  :: iw, iw0, iwp, iw0mw, iw0pw, i
  REAL(DP)               :: x1, x2       !ranges
  INTEGER                :: n            !number of points
  REAL(DP), ALLOCATABLE  :: x(:), w(:)   !abcissa and weights



  IF(.not.freq_gl) THEN
   zero = 0.0d0
   wgreenmin = wsigmamin-wcoulmax
   wgreenmax = wsigmamax+wcoulmax
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
  ! now find the correspondence between the arrays
  ! This is needed for the convolution G(w0-w)W(w) at the end
   allocate ( ind_w0mw (nwsigma,nwcoul), ind_w0pw (nwsigma,nwcoul) )
   do iw0 = 1, nwsigma
     do iw = 1, nwcoul
       w0mw = wsigma(iw0)-wcoul(iw)
       w0pw = wsigma(iw0)+wcoul(iw)
       foundp = .false.
       foundm = .false.
       do iwp = 1, nwgreen
         if (abs(w0mw-wgreen(iwp)) .lt. 1.d-10 ) then
           foundm = .true.
           iw0mw = iwp
         endif
         if ( abs(w0pw-wgreen(iwp)) .lt. 1.d-10 ) then
           foundp = .true.
           iw0pw = iwp
          endif
        enddo
       if ( ( .not. foundm ) .or. ( .not. foundp ) ) then
          call errore ('gwhs','frequency correspondence not found',1)
       else
          ind_w0mw(iw0,iw) = iw0mw
          ind_w0pw(iw0,iw) = iw0pw
       endif
     enddo
    enddo
  ELSE
! We generate Sigma on a uniform grid:

    nwgreen = 2*nwcoul
!   nwgreen = nwcoul !if we  use G(iw) = G^{*}(-iw)
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
      WRITE(stdout, '(6x, "Godby Needs Plasmon-Pole")')
  else if (padecont) then
      WRITE(stdout, '(6x, "Analytic Continuation")')
  else if (.not.padecont.and..not.godbyneeds) then
      WRITE(stdout, '(6x, "No screening model chosen!")')
  ENDIF
  WRITE(stdout, '(7x, "Imag. Frequencies: ")')
  DO i = 1, nfs
       WRITE(stdout,'(8x, i4, 4x, 2f9.4)')i, fiu(i)*RYTOEV
  ENDDO
  WRITE(stdout, '(5x, "Broadening: ", 1f10.4)'), eta
END SUBROUTINE freqbins
