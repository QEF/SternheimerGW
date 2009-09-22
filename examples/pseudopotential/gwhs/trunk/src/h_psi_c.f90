  !
  !-----------------------------------------------------------------------
  subroutine h_psi_c ( psi, hpsi, g2kin, vr)
  !-----------------------------------------------------------------------
  !
  ! hpsi = H_k | psi >
  ! complex version
  !
  !-----------------------------------------------------------------------
  !
  use gspace
  use constants
  use parameters
  implicit none
  !
  real(DP) :: g2kin(ngm)
  complex(kind=DP) :: vr (nr), psic (nr), psi (ngm), hpsi (ngm)
  integer :: ig, ir
  !
  ! Here we apply the kinetic energy (k+G)^2 psi
  !
  do ig = 1, ngm
    hpsi (ig) = g2kin (ig) * psi (ig)
  enddo
  !
  ! the empirical pseudopotential potential V * psi
  ! We compute the product in real space 
  !
  psic = czero
  do ig = 1, ngm
    psic ( nl ( ig ) ) = psi (ig)
  enddo
  call cfft3 ( psic, nr1, nr2, nr3,  1)
  !
  do ir = 1, nr
    psic (ir) = psic(ir) * vr (ir)
  enddo
  !
  ! back to G-space (fft order of G-vectors)
  call cfft3 ( psic, nr1, nr2, nr3, -1)
  !
  ! addition to the total product
  ! (automatically switch to magnitude-order of G-vectors)
  !
  do ig = 1, ngm
    hpsi(ig) = hpsi(ig) + psic( nl(ig) ) 
  enddo
  !
  return
  end subroutine h_psi_c
  !
