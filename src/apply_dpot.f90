!
! Copyright (C) 2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE apply_dpot(aux1, dvscfins, current_spin)
!
!  This routine applies the change of the self consistent potential to
!  one wavefunction
!
USE kinds, ONLY : DP
USE noncollin_module, ONLY : noncolin, npol, nspin_mag
USE spin_orb, ONLY : domag
!USE gsmooth, ONLY : dffts%nnr 
USE fft_base, ONLY : dffts
! USE gvect, ONLY : nrxx

IMPLICIT NONE
COMPLEX(DP), INTENT(IN) :: dvscfins(dffts%nnr,nspin_mag)
!COMPLEX(DP), INTENT(IN) :: dvbare(dffts%nnr)
COMPLEX(DP), INTENT(INOUT) :: aux1(dffts%nnr,npol)
INTEGER, INTENT(IN) :: current_spin

COMPLEX(DP) :: sup, sdwn
INTEGER :: ir

DO ir = 1, dffts%nnr
   aux1(ir,1)=aux1(ir,1)*dvscfins(ir,current_spin)
ENDDO

RETURN
END SUBROUTINE apply_dpot
