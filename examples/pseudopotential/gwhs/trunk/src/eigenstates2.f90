  !
  !----------------------------------------------------------------
  subroutine eigenstates2 ( xxk, vr, g2kin, psi, et )
  !----------------------------------------------------------------
  ! 
  ! find occupied eigenstates through conjugate gradients
  !
  !----------------------------------------------------------------
  !
  use parameters
  use constants
  use gspace
  implicit none
  !
  complex(dbl) :: psi(ngm, nbnd_occ), vr (nr)
  real(dbl) :: g2kin (ngm), xxk(3), et(nbnd_occ)
  real(dbl) :: precondition (ngm)
  !
  integer :: ideltag, ig1, ig2, ig3, ik, ig, ierr, ibnd
  real(dbl) :: arg, deltag(3)
  real(kind=DP), allocatable :: eval(:), u(:,:), fv1(:), fv2(:), hk(:,:), ss(:), vs(:)
  logical :: equiv
  integer :: notconv
  real(kind=DP) :: avg_iter
  !
  allocate ( eval(ngm0), u(ngm0,ngm0), fv1(ngm0), fv2(ngm0), & 
    ss(ngm), vs(ngm), hk (ngm0, ngm0) )
  !
  !
  ! STARTING WAVEFUNCTIONS FOR CONJUGATE GRADIENTS FROM DIRECT DIAGONALIZATION
  !
  ! the structure factor
  !
  do ig = 1, ngm
    arg = twopi * ( g(1,ig) * tau( 1, 1) + g(2,ig) * tau( 2, 1) + g(3,ig) * tau( 3, 1) )
    ss (ig) = cos ( arg )
  enddo
  !
  ! the empirical pseudopotential
  !
  vs = zero
  ! integer comparison - careful with other structures
  do ig = 1, ngm
    if     ( int ( gl(igtongl(ig)) ) .eq.  3 ) then
      vs (ig) =  v3
    elseif ( int ( gl(igtongl(ig)) ) .eq.  8 ) then
      vs (ig) =  v8
    elseif ( int ( gl(igtongl(ig)) ) .eq. 11 ) then
      vs (ig) =  v11
    endif
  enddo
  !
  do ig = 1, ngm0
    hk ( ig, ig) = g2kin ( ig )
  enddo
  !
  do ig1 = 1, ngm0
   do ig2 = ig1+1, ngm0
     !
     ! define ideltag
     ! 
     deltag = g(:,ig1) - g(:,ig2)
     ideltag = 1 
     do while ( .not. equiv ( deltag, g(:,ideltag) ) .and. ideltag.lt.ngm )
       ideltag = ideltag + 1
     enddo   
     !
     if (ideltag.ne.ngm) then
       hk ( ig1, ig2) = ss (ideltag) * vs (ideltag) 
       hk ( ig2, ig1) = hk ( ig1, ig2) 
     endif
     !
   enddo
  enddo
  !
  call rs ( ngm0, ngm0, hk, eval, 1, u, fv1, fv2, ierr)
  psi = czero
  do ibnd = 1, nbnd_occ
   do ig = 1, ngm0
      psi(ig,ibnd) = dcmplx ( u(ig,ibnd), zero)
   enddo
  enddo
  !
  ! CONJUGATE GRADIENTS DIAGONALIZATION - OCCUPIED STATES
  !
  precondition = max( 1.d0, g2kin )
  !
  call ccgdiagg (ngm, ngm, nbnd_occ, psi, et, precondition, eps, &
     maxter, .true., notconv, avg_iter, g2kin, vr)
  !
  deallocate ( eval, hk, u, fv1, fv2, ss, vs )
  !
  return
  end subroutine eigenstates2
  !----------------------------------------------------------------
  !





