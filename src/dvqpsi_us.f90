!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine dvqpsi_us (dvbarein, ik, uact, addnlcc)
  !----------------------------------------------------------------------
  !
  ! This routine calculates dV_bare/dtau * psi for one perturbation
  ! with a given q. The displacements are described by a vector u.
  ! The result is stored in dvpsi. The routine is called for each k point
  ! and for each pattern u. It computes simultaneously all the bands.
  ! It implements Eq. B29 of PRB 64, 235118 (2001). The contribution
  ! of the local pseudopotential is calculated here, that of the nonlocal
  ! pseudopotential in dvqpsi_us_only.
  !
  ! HL- FORMERLY this routine does what it says up there. I've commented out the stuff that pertains
  ! to initializing a phonon perturbation. And simply multiply the perturbing potential generated in coulomb
  ! a single plane wave by the eigenvector psi_nk
  !

  USE kinds, only : DP
  USE ions_base, ONLY : nat, ityp
  USE cell_base, ONLY : tpiba
  USE gvect,     ONLY : nrxx, eigts1, eigts2, eigts3, ig1,ig2,ig3, g, nl, &
                        ngm, nr1,nr2,nr3,nrx1,nrx2,nrx3
  USE gsmooth,   ONLY : nrxxs, ngms, doublegrid, nls, &
                        nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s
  USE lsda_mod,  ONLY : lsda, isk
  USE noncollin_module, ONLY : npol
  use uspp_param,ONLY : upf
  USE wvfct,     ONLY : nbnd, npw, npwx, igk
  USE wavefunctions_module,  ONLY: evc
  USE nlcc_gw,    ONLY : nlcc_any, drc
  USE eqv,        ONLY : dvpsi, dmuxc, vlocq
  USE qpoint,     ONLY : npwq, igkq, xq, eigqts, ikks

  implicit none
  !
  !   The dummy variables
  !

  integer :: ik, igpert
  ! input: the k point
  ! counter for G-vector for plane wave perturbation  

  complex(DP) :: uact (3 * nat)
  ! input: the pattern of displacements

  logical :: addnlcc
  !
  !   And the local variables
  !

  integer :: na, mu, ikk, ig, nt, ibnd, ir, is, ip
  ! counter on atoms
  ! counter on modes
  ! the point k
  ! counter on G vectors
  ! the type of atom
  ! counter on bands
  ! counter on real mesh

  complex(DP) :: gtau, gu, fact, u1, u2, u3, gu0
  complex(DP) , allocatable, target :: aux (:)
  complex(DP) , allocatable :: aux1 (:), aux2 (:)
  complex(DP) , pointer :: auxs (:)

  complex(DP) dvbarein(nrxxs)

  ! work space

  call start_clock ('dvqpsi_us')
  if (nlcc_any.and.addnlcc) then
     allocate (aux( nrxx))    
     if (doublegrid) then
        allocate (auxs( nrxxs))    
     else
        auxs => aux
     endif
  endif


  allocate (aux1( nrxxs))    
  allocate (aux2( nrxxs))    
  !
  !    We start by computing the contribution of the local potential.
  !    The computation of the derivative of the local potential is done in
  !    reciprocal space while the product with the wavefunction is done in
  !    real space
  !
  ikk = ikks(ik)
  dvpsi(:,:) = (0.d0, 0.d0)
  aux1(:) = (0.d0, 0.d0)

#if 0
  do na = 1, nat
     fact = tpiba * (0.d0, -1.d0) * eigqts (na)
     mu = 3 * (na - 1)
     if (abs (uact (mu + 1) ) + abs (uact (mu + 2) ) + abs (uact (mu + &
          3) ) .gt.1.0d-12) then
        nt = ityp (na)
        u1 = uact (mu + 1)
        u2 = uact (mu + 2)
        u3 = uact (mu + 3)
        gu0 = xq (1) * u1 + xq (2) * u2 + xq (3) * u3
        do ig = 1, ngms
           gtau = eigts1 (ig1 (ig), na) * eigts2 (ig2 (ig), na) * eigts3 ( &
                ig3 (ig), na)
           gu = gu0 + g (1, ig) * u1 + g (2, ig) * u2 + g (3, ig) * u3
           aux1 (nls (ig) ) = aux1 (nls (ig) ) + vlocq (ig, nt) * gu * &
                fact * gtau
        enddo
     endif
  enddo
#endif



   !HL INITIALIZE PERTURBATION dvbarein is in real space.

   aux1 (:) = dvbarein(:)

!   if (nlcc_any.and.addnlcc) then
!      aux(:) = (0.d0, 0.d0)
!      do na = 1,nat
!         fact = tpiba*(0.d0,-1.d0)*eigqts(na)
!         mu = 3*(na-1)
!         if (abs(uact(mu+1))+abs(uact(mu+2))  &
!                         +abs(uact(mu+3)).gt.1.0d-12) then
!            nt=ityp(na)
!            u1 = uact(mu+1)
!            u2 = uact(mu+2)
!            u3 = uact(mu+3)
!            gu0 = xq(1)*u1 +xq(2)*u2+xq(3)*u3
!            if (upf(nt)%nlcc) then
!               do ig = 1,ngm
!                  gtau = eigts1(ig1(ig),na)*   &
!                         eigts2(ig2(ig),na)*   &
!                         eigts3(ig3(ig),na)
!                  gu = gu0+g(1,ig)*u1+g(2,ig)*u2+g(3,ig)*u3
!                  aux(nl(ig))=aux(nl(ig))+drc(ig,nt)*gu*fact*gtau
!               enddo
!            endif
!         endif
!      enddo
!      call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,+1)
!      if (.not.lsda) then
!         do ir=1,nrxx
!            aux(ir) = aux(ir) * dmuxc(ir,1,1)
!         end do
!      else
!         is=isk(ikk)
!         do ir=1,nrxx
!            aux(ir) = aux(ir) * 0.5d0 *  &
!                 (dmuxc(ir,is,1)+dmuxc(ir,is,2))
!         enddo
!      endif
!     call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
!      if (doublegrid) then
!         auxs(:) = (0.d0, 0.d0)
!         do ig=1,ngms
!            auxs(nls(ig)) = aux(nl(ig))
!         enddo
!      endif
!      aux1(:) = aux1(:) + auxs(:)
!   endif
!HL- Compute Delta V_q(r')*Psi_nk(r') 

  do ibnd = 1, nbnd
     do ip=1,npol
        aux2(:) = (0.d0, 0.d0)
        if (ip==1) then
           do ig = 1, npw
              aux2 (nls (igk (ig) ) ) = evc (ig, ibnd)
           enddo
        else
           do ig = 1, npw
              aux2 (nls (igk (ig) ) ) = evc (ig+npwx, ibnd)
           enddo
        end if
        !
        !  This wavefunction is transformed into real space
        !

        call cft3s (aux2, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 2)
        do ir = 1, nrxxs
           aux2 (ir) = aux2 (ir) * aux1 (ir)
        enddo

        !
        ! and finally dV_loc/dtau * psi is transformed in reciprocal space
        !

        call cft3s (aux2, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, - 2)

        if (ip==1) then
           do ig = 1, npwq
              dvpsi (ig, ibnd) = aux2 (nls (igkq (ig) ) )
           enddo
        else
           do ig = 1, npwq
              dvpsi (ig+npwx, ibnd) = aux2 (nls (igkq (ig) ) )
           enddo
        end if
     enddo
  enddo

!HL FFT how it works? 
! This should recover the perturbing potential in G space. 
!         write (6, '("")' )
!         do ig = 1, 3
!         write (6, '("dvbare(nl(igpert)) ", 3f7.4, 3f7.3)' )aux1(nl(ig))
!         end do
!         call cft3 (aux1, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1)
!         call cft3 (aux1, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -1)
!         write(6, '("")')
!         do ig = 1, 3
!         write (6, '("dvbare(nl(igpert)) ", 3f7.4, 3f7.3)' )aux1(nl(ig))
!         write (6, '("dvbare(nl(igpert)) ", 3f7.4, 3f7.3)' )aux1(ig)
!         end do
! HL This does recover the correct array of fourier components. I think drho should transform in exactly the same way.



  deallocate (aux2)
  deallocate (aux1)
  if (nlcc_any.and.addnlcc) then
     deallocate (aux)
     if (doublegrid) deallocate (auxs)
  endif

!  We add the contribution of the nonlocal potential in the US form
!  First a term similar to the KB case.
!  Then a term due to the change of the D coefficients. 
!  HL We need to add anything else to the right hand side.
!  No change in Q, Beta, or S of the self-consistent potential. Just the shift due to 
!  electrons. 
!  HL call dvqpsi_us_only (ik, uact)

  call stop_clock ('dvqpsi_us')
  return
end subroutine dvqpsi_us
