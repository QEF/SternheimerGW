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
subroutine deallocate_gwq
!----------========-----------------------
!
!  deallocates the variables allocated by allocate_gwq
!
  USE becmod,                ONLY: bec_type, becp, deallocate_bec_type
  USE control_lr,            ONLY: lgamma, nbnd_occ
  USE eqv,                   ONLY: dmuxc, vlocq, dpsi, dvpsi, evq
  USE eqv_gw,                ONLY: dpsim, dpsip, dvbare
  USE lrus,                  ONLY: becp1
  USE qpoint,                ONLY: eigqts, igkq, ikks, ikqs

  IMPLICIT NONE
  INTEGER :: ik


  !IMPORTANT igk arrays need to be nullified/deallocated???
  if (lgamma) then
     if(associated(evq)) nullify(evq)
     if(associated(igkq)) nullify(igkq)
  else
     if(associated(evq)) deallocate(evq)
     if(associated(igkq)) deallocate(igkq)
  end if

  if(allocated(dvpsi)) deallocate (dvpsi)    
  if(allocated(dpsi)) deallocate ( dpsi)    
  if(allocated(vlocq)) deallocate (vlocq)
  if(allocated(dmuxc)) deallocate (dmuxc)
  if(allocated(ikks)) deallocate (ikks)
  if(allocated(ikqs)) deallocate (ikqs)
  if(allocated(eigqts)) deallocate (eigqts)
  !HL \psi\pm W,G, dvbare 
  if(allocated(dpsim)) deallocate(dpsim)
  if(allocated(dpsip)) deallocate(dpsip)
  if(allocated(dvbare)) deallocate(dvbare)

  if(allocated(becp1))  then
     do ik=1,size(becp1)
        call deallocate_bec_type ( becp1(ik) )
     end do
     deallocate(becp1)
  end if
  call deallocate_bec_type ( becp )

  if(allocated(nbnd_occ)) deallocate(nbnd_occ)

  return
end subroutine deallocate_gwq
