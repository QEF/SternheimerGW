  !
  !-----------------------------------------------------------------------
  subroutine ch_psi_all ( h, ah, e, ik, nb, g2kin, vr, evq)
  !-----------------------------------------------------------------------
  !
  ! This routine applies the operator ( H - \epsilon S + alpha_pv P_v)
  ! to a vector h. The result is given in Ah.
  !
  !-----------------------------------------------------------------------
  !
  use parameters
  use gspace
  use constants
  implicit none
  !
  integer :: ik, nb
  ! input: the dimension of h
  ! input: the number of bands
  ! input: the k point
  real(kind=DP) :: e (nb)
  ! input: the eigenvalue
  complex(kind=DP) :: h (ngm, nb), ah (ngm, nb)
  ! input: the vector
  ! output: the operator applied to the vector
  !
  !   local variables
  !
  integer :: ibnd, ikq, ig
  ! counter on bands
  ! the point k+q
  ! counter on G vetors
  complex(kind=DP) :: ps (nb,nb)
  ! scalar products
  complex(kind=DP), allocatable :: hpsi (:,:)
  ! the product of the Hamiltonian and h
  complex(kind=DP) :: vr(nr)
  real(kind=DP) :: g2kin(ngm)
  complex(kind=DP) :: evq (ngm, nbnd_occ)
  !
  allocate ( hpsi( ngm , nb) )    
  hpsi = czero
  !
  !   compute the product of the hamiltonian with the h vector
  !
  do ibnd = 1, nb
    call h_psi_c ( h(:,ibnd), hpsi(:,ibnd), g2kin, vr)
  enddo
  !
  !   then we compute the operator H-epsilon S
  !
  do ibnd = 1, nb
     do ig = 1, ngm
        ah (ig, ibnd) = hpsi (ig, ibnd) - e (ibnd) * h (ig, ibnd)
     enddo
  enddo
  !
  !   Here we compute the projector in the valence band
  !  
  ! important: don't forget dcmplx and czero in the second call 
  !
  call ZGEMM ('C', 'N', nb , nb, ngm, (1.d0, 0.d0) , evq, &
       ngm, h, ngm, czero , ps, nb)
  call ZGEMM ('N', 'N', ngm, nb, nb, dcmplx(alpha_pv,0.d0), evq, &
       ngm, ps, nb, czero, hpsi, ngm)   
  !
  do ibnd = 1, nb 
     do ig = 1, ngm
        ah (ig, ibnd) = ah (ig, ibnd) + hpsi (ig, ibnd)
     enddo
  enddo
  !
  deallocate (hpsi)
  !
  return
  end subroutine ch_psi_all
  !-----------------------------------------------------------------------
  !
