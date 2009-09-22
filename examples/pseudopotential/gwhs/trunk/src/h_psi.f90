  !
  !-----------------------------------------------------------------------
  subroutine h_psi ( psi, hpsi, g2kin, vr)
  !-----------------------------------------------------------------------
  !
  ! hpsi = H_k | psi >
  ! FG - tested vs. \sum_G' H(G,G') * psi(G') - works fine
  !
  use gspace
  use constants
  use parameters
  implicit none
  !
  real(kind=DP) :: psi (ngm), hpsi (ngm), g2kin(ngm)
  complex(kind=DP) :: vr (nr), psic (nr), rhoca(nr1,nr2,nr3)
  integer :: ig, ir, i, j, k
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
    psic ( nl ( ig ) ) = dcmplx( psi (ig), zero )
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
    hpsi(ig) = hpsi(ig) + real ( psic( nl(ig) ) ) 
  enddo
  !
  return
  end subroutine h_psi
  !
