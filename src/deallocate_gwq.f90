!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! Parts of this file are taken from the Quantum ESPRESSO software
! P. Giannozzi, et al, J. Phys.: Condens. Matter, 21, 395502 (2009)
!
! Copyright (C) 2010 - 2016 Quantum ESPRESSO group,
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
subroutine deallocate_gwq
!----------========-----------------------
!
!  deallocates the variables allocated by allocate_gwq
!
  USE noncollin_module,      ONLY: m_loc
  USE becmod,                ONLY: bec_type, becp, deallocate_bec_type
  USE wavefunctions_module,  ONLY: evc
  USE qpoint,                ONLY: eigqts, igkq, ikks, ikqs, nksq
  USE gc_lr,                 ONLY: grho, gmag, dvxc_rr,  dvxc_sr,  dvxc_ss, dvxc_s, &
                                   vsgga, segni
  USE gamma_gamma,           ONLY : with_symmetry, has_equivalent, equiv_atoms, &
                                    n_equiv_atoms
  USE eqv_gw,                ONLY : dmuxc, vlocq, dpsi, dvpsi, evq, eprec, eprectot, &
                                    dpsim, dpsip, dvbare
  USE nlcc_gw,               ONLY : drc
  USE units_gw,              ONLY : this_dvkb3_is_on_file, this_pcxpsi_is_on_file
  USE control_gw,            ONLY : lgamma
  USE lrus,                  ONLY : becp1

  IMPLICIT NONE
  INTEGER :: ik, ipol

  if(allocated(dvpsi)) deallocate (dvpsi)    
  if(allocated(dpsi)) deallocate ( dpsi)    
  if(allocated(vlocq)) deallocate (vlocq)
  if(allocated(dmuxc)) deallocate (dmuxc)
  if(allocated(eprec)) deallocate (eprec)
  if(allocated(eprectot)) deallocate (eprectot)
  if(allocated(ikks)) deallocate (ikks)
  if(allocated(ikqs)) deallocate (ikqs)
  if(allocated(eigqts)) deallocate (eigqts)
  !HL \psi\pm W,G, dvbare 
  if(allocated(dpsim)) deallocate(dpsim)
  if(allocated(dpsip)) deallocate(dpsip)
  if(allocated(dvbare)) deallocate(dvbare)
  if(allocated(drc)) deallocate(drc)
  if(allocated(this_dvkb3_is_on_file)) deallocate (this_dvkb3_is_on_file)    
  if(allocated(this_pcxpsi_is_on_file)) deallocate (this_pcxpsi_is_on_file)

  if(allocated(becp1))  then
     do ik=1,size(becp1)
        call deallocate_bec_type ( becp1(ik) )
     end do
     deallocate(becp1)
  end if
  call deallocate_bec_type ( becp )


  return
end subroutine deallocate_gwq
