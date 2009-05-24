  !
  !-----------------------------------------------------------------------
  subroutine cgdiag (nbnd, psi, e, g2kin, vr)
  !-----------------------------------------------------------------------
  !
  ! FG - wrapper to the complex ccgdiagg of espresso
  !
  use gspace, only : ngm, nr
  use parameters, only : DP, eps, maxter
  use constants, only : pi, tpiba2
  implicit none
  !
  logical :: reorder
  integer :: notconv, nbnd
  real(kind=DP) :: psi (ngm, nbnd)
  real(kind=DP) :: e (nbnd), precondition (ngm), g2kin(ngm)
  complex(kind=DP) :: vr(nr), psic (ngm, nbnd)
  real(kind=DP) :: avg_iter
  !
  ! preconditioning matrix
  !
  precondition = max( 1.d0, g2kin ) 
  !
  psic = dcmplx (psi, 0.d0)
  !
  call ccgdiagg (ngm, ngm, nbnd, psic, e, precondition, eps, &
     maxter, .true., notconv, avg_iter, g2kin, vr) 
  !
  psi = dreal (psic)
  !
!  write(6,*) avg_iter
  !
  return
  end subroutine cgdiag
  !
