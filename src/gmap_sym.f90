!------------------------------------------------------------------------------
!
! This file is part of the SternheimerGW code.
! 
! Copyright (C) 2010 - 2018 Jesse Noffsinger, Brad Malone,
! Henry Lambert, Martin Schlipf, and Feliciano Giustino
!
! SternheimerGW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! SternheimerGW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with SternheimerGW. If not, see
! http://www.gnu.org/licenses/gpl.html .
!
!------------------------------------------------------------------------------ 
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE gmap_sym (num_g_corr, nsym, s, ftau, gmapsym, eigv, invs)
  !-----------------------------------------------------------------------
  !
  !  For every G vector, find S(G) for all the symmetry operations
  !  of the crystal. Construct the matrix
  !  eigv(ig,isym) = e^{i G v(S)} where v(S) is the (possible) 
  !  fractional translation associated with the symmetry operation
  !
  !  No parallelization on G-vecs at the moment  
  !  (actually this is done on the global array, but in elphel2.f90
  !  every processor has just a chunk of the array, I may need some
  !  communication)
  !
  !  No ultrasoft now
  !  HL: Fractional rotations not applicable if FFT grid is different from full wave
  !  function grid!
  !
  !
  !----------------------------------------------------------------------

  USE fft_base,      ONLY : dfftp
  USE gvect,         ONLY : mill, ngm
  USE kinds,         ONLY : DP
  USE mp_world,      ONLY : mpime

  implicit none
  !
  ! input variables
  !
  !> number of G vectors in correlation grid
  INTEGER, INTENT(IN) :: num_g_corr
  !
  integer :: nsym, s(3,3,48), ftau(3,48), invs(48)
  ! the number of symmetries of the crystal
  ! the symmetry matrices
  ! the fractional traslations
  !
  ! output variables
  !
  !integer :: gmapsym (ngm, 48)
  integer :: gmapsym (num_g_corr, nsym)
  ! the map S(G) = gmapsym (G,S) 1...nsym
  !complex(kind=DP) :: eigv (ngm, 48)
  complex(kind=DP) :: eigv (num_g_corr, nsym)
  ! e^{ iGv} for 1...nsym
  !
  ! local variables
  !
  integer :: ig, jg, i, j, k, notfound, isym, ism1
  logical :: tfound
  real(kind=DP), PARAMETER :: twopi = 6.28318530717959
  real(kind=DP) :: rdotk
  complex(kind=DP), PARAMETER :: ci = (0.d0,1.d0), &
     cone = (1.d0, 0.d0)

  !
  !  loop on the symmetries of the crystal
  !
  DO isym = 1, nsym
    !
    ism1 = invs(isym)
    !
    ! loop on the G vectors 
    !
    notfound = 0
    !
    DO ig = 1, num_g_corr
      !
      !  the rotated G-vector
      !
      i = s (1, 1, isym) * mill(1,ig) + s (1, 2, isym) * mill(2,ig) + s (1, 3, isym) * mill(3,ig)
      j = s (2, 1, isym) * mill(1,ig) + s (2, 2, isym) * mill(2,ig) + s (2, 3, isym) * mill(3,ig)
      k = s (3, 1, isym) * mill(1,ig) + s (3, 2, isym) * mill(2,ig) + s (3, 3, isym) * mill(3,ig)
      !
      jg = 0
      tfound = .false.
      DO while ((.not.tfound).and.(jg.lt.ngm))
        jg = jg + 1
        tfound = (i.eq.mill(1,jg)) .and. (j.eq.mill(2,jg)) .and. (k.eq.mill(3,jg))
      ENDDO
      !
      IF (tfound) THEN
        gmapsym ( ig, isym ) = jg
      ELSE
        gmapsym ( ig, isym ) = 0
        notfound = notfound + 1
      ENDIF
      !
      ! now the phase factors e^{iGv}
      !
      IF ( ftau (1, isym).ne.0 .or. ftau (2, isym).ne.0 .or. ftau (3, isym).ne.0 ) THEN
!      if(ig.eq.1)  write(6,'("FTAU")')
!      if(ig.eq.1)  write(6,*) isym, ftau(:,isym)
!
!      fractional traslation in crystal coord is ftau/nr*
!      for cart/crys transform of the G-vecctors have a look at the bottom
!
      rdotk = float( mill(1,ig) * ftau (1, isym) ) / float (dfftp%nr1) &
            + float( mill(2,ig) * ftau (2, isym) ) / float (dfftp%nr2) &
            + float( mill(3,ig) * ftau (3, isym) ) / float (dfftp%nr3)

        !     ft(:)    = at(:,1)*ftau(1,ns)/nr1 + &
        !                at(:,2)*ftau(2,ns)/nr2 + &
        !                at(:,3)*ftau(3,ns)/nr3

        !     rdotk    = g(1,ig) * ft(1) + &
        !                g(2,ig) * ft(2) + &
        !                g(3,ig) * ft(3)        !

        ! the actual translation is -v (have a look at ruota_ijk.f90)
        ! 
        eigv (ig, isym) = exp( - ci*twopi*rdotk ) 
      ELSE
        eigv (ig, isym) = cone
      ENDIF
        !
    ENDDO
    !
    IF (notfound.gt.0) then
       write(6,'("WRITE GVEC not rotated properly gmap_sym.")')
       write(1000+mpime,'("WRITE GVEC not rotated properly gmap_sym.")')
       stop
    ENDIF
    ! CALL errore ('gmap_sym','incomplete mapping of G vectors: notfound = ',notfound)
    !
  ENDDO
  !
  END SUBROUTINE gmap_sym
