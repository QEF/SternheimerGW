!------------------------------------------------------------------------------
!
! This file is part of the SternheimerGW code.
! Parts of this file are taken from the Quantum ESPRESSO software
! P. Giannozzi, et al, J. Phys.: Condens. Matter, 21, 395502 (2009)
!
! Copyright (C) 2010 - 2018 Quantum ESPRESSO group,
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
subroutine star_q (xq, at, bg, nsym, s, invs, nq, sxq, isq, nsq, imq, verbosity )
  !FROM DRESSELHAUS:
  !The  irreducible representation D(S) of a space group S is
  !characterized by the star of wave vectors nad the irreducible
  !representation Gamma(H) of the point group H(K_{1}). Its dimension is
  !m = qd where q is the index of the star of wave vectors (HL nqs), and
  !d is the dimension of the representation \Gamma.  
  !for instance the 0 0 0 kpoint has nqs = 1 and its representation is the
  !same as the whole crystal.
  !-----------------------------------------------------------------------
  ! generate the star of q vectors that are equivalent to the input one
  ! NB: input s(:,:,1:nsym) must contain all crystal symmetries,
  ! i.e. not those of the small-qroup of q only
  !
  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  implicit none
  !
  real(DP), parameter :: accep=1.e-5_dp

  integer, intent(in) :: nsym, s (3, 3, 48), invs(48)
  ! nsym matrices of symmetry operations
  ! invs: list of inverse operation indices
  real(DP), intent(in) :: xq (3), at (3, 3), bg (3, 3)
  ! xq: q vector
  ! at: direct lattice vectors
  ! bg: reciprocal lattice vectors
  !
  integer, intent(out) :: nq, isq (48), imq
  ! nq  : degeneracy of the star of q
  ! isq : index of q in the star for a given sym
  ! imq : index of -q in the star (0 if not present)

  real(DP), intent(out) :: sxq (3, 48)
  ! list of vectors in the star of q
  logical, intent(in) :: verbosity
  ! if true prints several messages.
  !
  integer :: nsq (48), isym, ism1, iq, i
  ! number of symmetry ops. of bravais lattice
  ! counters on symmetry ops.
  ! index of inverse of isym
  ! counters
  real(DP) :: saq (3, 48), aq (3), raq (3), zero (3)
  ! auxiliary list of q (crystal coordinates)
  ! input q in crystal coordinates
  ! rotated q in crystal coordinates
  ! coordinates of fractionary translations
  ! a zero vector: used in eqvect
  logical, external :: eqvect
  ! function used to compare two vectors
  !
  zero(:) = 0.d0
  !
  ! go to  crystal coordinates
  !
  do i = 1, 3
     aq(i) = xq(1) * at(1,i) + xq(2) * at(2,i) + xq(3) * at(3,i)
  enddo
  !
  ! create the list of rotated q
  !
  do i = 1, 48
     nsq (i) = 0
     isq (i) = 0
  enddo
  nq = 0
  do isym = 1, nsym
     ism1 = invs (isym)
     do i = 1, 3
        raq (i) = s (i, 1, ism1) * aq (1) &
                + s (i, 2, ism1) * aq (2) &
                + s (i, 3, ism1) * aq (3)
     enddo
     do i = 1, 3
        sxq (i, 48) = bg (i, 1) * raq (1) &
                    + bg (i, 2) * raq (2) &
                    + bg (i, 3) * raq (3)
     enddo
     do iq = 1, nq
        if (eqvect (raq, saq (1, iq), zero, accep) ) then
           isq (isym) = iq
           nsq (iq) = nsq (iq) + 1
        endif
     enddo
     if (isq (isym) == 0) then
        nq = nq + 1
        nsq (nq) = 1
        isq (isym) = nq
        saq(:,nq) = raq(:)
        do i = 1, 3
           sxq (i, nq) = bg (i, 1) * saq (1, nq) &
                       + bg (i, 2) * saq (2, nq) &
                       + bg (i, 3) * saq (3, nq)
        enddo
     endif
  enddo
  !
  ! set imq index if needed and check star degeneracy
  !
  raq (:) = - aq(:)
  imq = 0
  do iq = 1, nq
     if (eqvect (raq, saq (1, iq), zero, accep) ) imq = iq
     if (nsq(iq)*nq /= nsym) call errore ('star_q', 'wrong degeneracy', iq)
  enddo
  !
  ! writes star of q
  !
  IF (verbosity) THEN
  WRITE( stdout, * )
  WRITE( stdout, '(5x,a,i4)') 'Number of q in the star = ', nq
  WRITE( stdout, '(5x,a)') 'List of q in the star:'
  WRITE( stdout, '(7x,i4,3f14.9)') (iq, (sxq(i,iq), i=1,3), iq=1,nq)
  if (imq == 0) then
     WRITE( stdout, '(5x,a)') 'In addition there is the -q list: '
     WRITE( stdout, '(7x,i4,3f14.9)') (iq, (-sxq(i,iq), i=1,3), iq=1,nq)
  endif
  ENDIF
  return
end subroutine star_q

