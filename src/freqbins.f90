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

  USE freq_gw,    ONLY : nwcoul, nwgreen, nwalloc, nwsigma, wtmp, wcoul,& 
                         wgreen, wsigma, wsigmamin, wsigmamax,&
                         deltaw, wcoulmax, ind_w0mw, ind_w0pw, wgreenmin,&
                         wgreenmax, fiu, nfs, greenzero
  USE kinds,      ONLY : DP
  USE constants,  ONLY : RYTOEV
  USE control_gw, ONLY : eta
           

  IMPLICIT NONE 

  LOGICAL  :: foundp, foundm
  REAL(DP) :: zero, w0mw, w0pw
  INTEGER  :: iw, iw0, iwp, iw0mw, iw0pw, i

!  wsigmamin = -14.d0 
!  wsigmamax =  28.d0 
!  deltaw    = 0.25d0
!  wcoulmax  = 80.d0


   greenzero      = 0.0d0

   wgreenmin = wsigmamin-wcoulmax
   wgreenmax = wsigmamax+wcoulmax

   nwalloc = 1 + ceiling( (wgreenmax-wgreenmin) / deltaw )

   allocate(wtmp(nwalloc), wcoul(nwalloc), wgreen(nwalloc), wsigma(nwalloc) )

   wcoul = greenzero
   wgreen = greenzero
   wsigma = greenzero

  do iw = 1, nwalloc
    wtmp(iw) = wgreenmin + (wgreenmax-wgreenmin)/float(nwalloc-1)*float(iw-1)
  enddo

 !align the bins with the zero of energy

  wtmp = wtmp - minval ( abs ( wgreen) )

 !
  nwgreen = 0
  nwcoul = 0
  nwsigma = 0

 
  do iw = 1, nwalloc
   if ( ( wtmp(iw) .ge. wgreenmin ) .and. ( wtmp(iw) .le. wgreenmax) ) then
     nwgreen = nwgreen + 1
     wgreen(nwgreen) = wtmp(iw)
   endif

   if ( ( wtmp(iw) .ge. greenzero ) .and. ( wtmp(iw) .le. wcoulmax) ) then
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
        if ( abs(w0mw-wgreen(iwp)) .lt. 1.d-10 ) then
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
   
 !Print out Frequencies on Imaginary Axis for reference...
  WRITE(6,'("wcoulmax:")')
  write(6,*) wcoulmax
  WRITE(6,'("eta")')
  WRITE(6,'(1f10.4 )') eta
  WRITE(6,'("wsigmamin, wsigmamax, deltaw:")')
  WRITE(6,'(3f10.4 )') wsigmamin, wsigmamax, deltaw 
  WRITE(6,'(//5x, "Imag. Frequencies: ")')
  DO i = 1, nfs
       WRITE(6,'(i4, 4x, 2f9.4)')i, fiu(i)*RYTOEV
  ENDDO

END SUBROUTINE freqbins
