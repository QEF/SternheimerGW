  !-----------------------------------------------------------------------
  ! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !-----------------------------------------------------------------------
  !
  !----------------------------------------------------------------
subroutine refold ( )
  !----------------------------------------------------------------
  ! the map is defined as follows
  ! g(:,gmap(ig,ig0)) = g(:,ig) - g0vec(:,ig0)
  ! at the exit, the folding vectors are in cartesian coordinates
  !----------------------------------------------------------------
  !
  USE kinds,           ONLY : DP
  USE gvect,           ONLY : ngm, g
  USE gsmooth,         ONLY : ngms
  USE disp,            ONLY : g0vec, gmap, eval_occ
  USE cell_base,       ONLY : at
  USE klist,           ONLY : nkstot
  USE wvfct,           ONLY : nbnd

  IMPLICIT NONE
  !
  integer ::  notfound, g2(3), count
  integer :: ig0, ig, igp, ig1, ig2, ig3, ng0vec
  real(DP) :: g1(3), bg (3, 3) 

  !
  ! the 3^3 possible translations (in very anisotropic materials
  ! this could go to 5^3 or even more - see EPW code)
  ! crystal coordinates

   ALLOCATE (gmap (ngm, 27))
   ALLOCATE (eval_occ(nbnd, nkstot))

  ng0vec = 0
  do ig1 = -1,1
    do ig2 = -1,1
      do ig3 = -1,1
        ng0vec = ng0vec + 1
        g0vec(1,ng0vec) = ig1
        g0vec(2,ng0vec) = ig2
        g0vec(3,ng0vec) = ig3
      enddo
    enddo
  enddo
  !
  ! bring all the G-vectors in crystal coordinates
  ! (in real codes we use the array of the miller indices mill_ )
  !  
  call cryst_to_cart (ngm, g, at, -1)

  count = 0
  notfound = 0
  do ig = 1, ngm
    !
    do ig0 = 1, ng0vec
      !
      gmap (ig,ig0) = 0 
      count = count + 1
      !
      g2 = nint ( g(:,ig) - g0vec(:,ig0) )
      !
      do igp = 1, ngm
        if ( ( g2(1) .eq. nint(g(1,igp)) ) .and. &
             ( g2(2) .eq. nint(g(2,igp)) ) .and. &
             ( g2(3) .eq. nint(g(3,igp)) ) )     &
              gmap (ig,ig0) = igp
      enddo
      !write(6,*) gmap (ig,ig0)
      if (gmap (ig,ig0).eq.0) notfound = notfound + 1
    enddo
  enddo
  write(6,'(4x,"refold: notfound = ",i6," out of ",i6)') notfound,count

!HLDEBUG
!  do ig=1,ngm
!     do ig0=1, ng0vec
!         write(6,*), gmap (ig,ig0) 
!     enddo
!  enddo
!ENDEBUG

  !
  ! back to cartesian
  !  
  call cryst_to_cart (ngm, g, bg, 1)
  call cryst_to_cart (ng0vec, g0vec, bg, 1)
  !
end subroutine refold
