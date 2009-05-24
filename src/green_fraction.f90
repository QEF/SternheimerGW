  !
  !-----------------------------------------------------------------------
  subroutine green_fraction ( ik, nw, w, greenf )
  !-----------------------------------------------------------------------
  !
  ! Green's function from Haydock's recursion: G_k (G,G',w) for a given k
  ! the k-dependence is inside the kinetic energy g2kin 
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
  integer :: ig1, ig2, ih, ik
  !
  ! notation as in Haydock - review (pag 227)
  ! e.g.: bnpunp = b_{n+1} * u_{n+1}, hun = H * u_n
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
  complex(dbl) :: greenf (nrs,nrs,nw)
  !
  integer :: ig, iw, nw
  real(DP) :: w(nw), w_ryd(nw)
  complex(DP) :: gr, gr_tmp
  !
  prefac = (/ cone, -cone, ci, -ci/)
  w_ryd = w/ryd2ev
  !
  greenf = czero 
  !
  do ig1 = 1, 2 !@ ngms
   do ig2 = 1, 2 !@ ngms
     !
     !
     call rw_haydock ( ik, ig1, ig2, a, b, a_term, b_term, norm, -1)
     !
     do iw = 1, nw
       !
       gr = czero
       !
       do ih = 1, 4  
         !
         call recfrac ( w_ryd(iw), a(:,ih), b(:,ih), a_term(ih), b_term(ih), gr_tmp)
         !
         ! in the case ih=2 we need to remove the spurious 
         ! solution corresponding to ig1=ig2 (phi = d1-d2 = 0)
         !
         if ( ih.ne.2 .or. ig1.ne.ig2 ) &
           gr = gr + prefac(ih) * gr_tmp * norm(ih) * norm(ih)
         !
       enddo
       !
       ! output in (eV ; 1/eV)
       write(201,'(3f15.10)')  w(iw), gr / ryd2ev
       !
       greenf(ig1,ig2,iw) = gr
       !
     enddo
     !
   enddo
  enddo
  !
  return
  end subroutine green_fraction
  !
