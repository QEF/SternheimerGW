  !
  !-----------------------------------------------------------------------
  subroutine green_coeff ( ik, g2kin, vr, nw, w)
  !-----------------------------------------------------------------------
  !
  ! Green's function from Haydock's recursion: G_k (G,G',w) for a given k
  ! the k-dependence is inside the kinetic energy g2kin 
  !
  ! for some reasons the dot_product(conjg(),()) does not work
  ! need to check this out
  !
  ! ik is the index of this k-point in the {k0-q} list
  !
  !-----------------------------------------------------------------------
  !
  use gspace
  use constants
  use parameters
  implicit none
  !
  complex(kind=DP) :: vr (nr)
  real(kind=DP) :: g2kin (ngm)
  integer :: ig1, ig2, i, j, k, ir, ih, ik, nw
  real(DP) :: w(nw)
  !
  ! notation as in Haydock - review (pag 227)
  ! e.g.: bnpunp = b_{n+1} * u_{n+1}, hun = H * u_n
  !
  complex(DP) :: d1 (ngm), d2 (ngm), phi (ngm)
  !
  real(DP) :: a(nstep+2,4), b(nstep+2,4), norm(4), a_term(4), b_term(4)
  ! a(:,1) are the coefficients for <u|G|u>
  ! a(:,2) are the coefficients for <v|G|v>
  ! ...
  ! norm(1) is the norm <u|u>^1/2 
  ! ...
  ! a_term(1) is the constant coefficient for the terminator for <u|G|u>
  ! ...
  !
  complex(DP) :: prefac(4)
  !
  complex(DP) :: ZDOTC
  real(DP) :: norm1
  integer :: ig, iw, ibnd
  complex(DP) :: gr, gr_tmp
  !
  ! variables for debugging
  !
  real(DP) :: eval(nbnd)
  complex(DP) :: ctmp1, ctmp2, gr_exp(nw), psi(ngm,nbnd)
  !
  prefac = (/ cone, -cone, ci, -ci/)
  !
  a = zero
  b = zero
  a_term = zero
  b_term = zero
  norm = zero
  !
  do ig1 = 1, 2 !@ ngms
   do ig2 = 1, 2 !@ ngms
     !
     ! define planewaves for projecting the Green's function
     !
     ! d1 = exp(i*G1*r), d2 =  exp(i*G2*r)
     !
     d1 = czero
     d1 ( ig1 ) = cone 
     d2 = czero
     d2 ( ig2 ) = cone 
     !
     ! check normalization, should be 1 (OK)
     !
     ! norm1 = DREAL ( ZDOTC (ngm, d1, 1, d1, 1) )
     ! norm1 = DREAL ( ZDOTC (ngm, d2, 1, d2, 1) )
     !
     ! DEBUG - calculate the Green's function using the sum over states
     ! <d1| G(r,r',k,w) |d2>
     !
     call  eigenstates_all ( vr, g2kin, psi, eval )
     gr_exp = czero
     do ibnd = 1, nbnd
       ctmp1 = ZDOTC (ngm, psi(:,ibnd), 1, d1, 1)
       ctmp2 = ZDOTC (ngm, psi(:,ibnd), 1, d2, 1)
       do iw = 1, nw
         gr_exp(iw) = gr_exp(iw) + conjg(ctmp1)*ctmp2 / ( w(iw) - eval(ibnd)*ryd2ev + ci * eta*ryd2ev)
       enddo
     enddo
     do iw = 1, nw
       write(200,'(3f15.10)') w(iw), gr_exp(iw) 
     enddo
     !
     ! Haydock's formula for off-diagonal elements [insert ref. to book]
     !
     ! G_12 = G_uu - G_vv - i * (G_ww - G_zz) 
     ! G_21 = G_uu - G_vv + i * (G_ww - G_zz) 
     !                   ^^^
     ! u = 0.5 * (d1+ d2)
     ! v = 0.5 * (d1- d2)
     ! w = 0.5 * (d1+id2)
     ! z = 0.5 * (d1-id2)
     !
     ! NOTE: when ig1=ig2 the procedure at (d1-d2)/2
     ! finds a nihil eigenvalue and the
     ! Green's function is peaked at 0 - spurious
     !
     ! NORMALIZATION: If I use a non-normalized function
     ! Haydock does not work. Keep this in mind for the
     ! local-orbital implementation
     !
     ! the 4 Haydock sequences
     !
     do ih = 1, 4
       if (ih.ne.2 .or. ig1.ne.ig2) then
         !
         ! in the case ih=2 we need to remove the spurious 
         ! solution corresponding to ig1=ig2 (phi = d1-d2 = 0)
         !
         phi = 0.5d0 * ( d1 + prefac(ih) * d2 ) 
         call normalize ( phi, norm(ih) )
         call haydock ( g2kin, vr, phi, a(:,ih), b(:,ih), a_term(ih), b_term(ih))
       endif
     enddo
     !
     call rw_haydock ( ik, ig1, ig2, a, b, a_term, b_term, norm, +1)
     !
     !
   enddo
  enddo
  !
  return
  end subroutine green_coeff
  !
