!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine addusddens (drhoscf, dbecsum, mode0, npe, iflag)
  !----------------------------------------------------------------------
  !
  !  This routine adds to the change of the charge and of the
  !  magnetization densities the part due to the US augmentation.
  !  It assumes that the array dbecsum has already accumulated the
  !  change of the becsum term. It calculates Eq. B31 of Ref [1].
  !  If called from drho (iflag=1), dbecsum and drhoscf contain the
  !  orthogonalization contribution to the change of the wavefunctions
  !  and the terms with alphasum and becsum are added. If called 
  !  from solve_linter (iflag=0) drhoscf and dbecsum contains the contribution 
  !  of the solution of the linear system and the terms due to alphasum
  !  and becsum are not added. In this case the change of the charge 
  !  calculated by drho (called \Delta \rho in [1]) is read from file 
  !  and added. The contribution of the change of 
  !  the Fermi energy is not calculated here but added later by ef_shift.
  !  [1] PRB 64, 235118 (2001).
  !  Eq. A17 
  !  d\rho/du = 2 \sum_{kv}\psi^{*}_{kv}(\r)\Delta\psi_{kv}(\r) + 2 \sum_{kv}\sum_{Inm} Q_{nm}(r-R_{I})
  !              <\psi_{k\v}|\Beta^{I}_{n}> <\Beta^{I}_{m}|\Delta\psi_{v\sigma}> 
  !  Eq. B31
  !  d\rho/du =  2 FT_{q+G} \sum_{kv}psi(\r) \delta \psi(r) +  \sum_{s_{1}}\sum_{nm}Q_{nm}(q+G)
  !                e^{-i(q+G)\tau_{s_{1}}*[a_{s1,nm}]
  !  a_{s_{1}nm} = \frac{2}{N} \sum_{kv}


  USE kinds,                 ONLY : DP
  USE gvect,                 ONLY : gg, ngm, nrxx, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                                    nl, g, eigts1, eigts2, eigts3, ig1, ig2, ig3
  USE uspp,                  ONLY : okvan, becsum
  USE cell_base,             ONLY : tpiba
  USE ions_base,             ONLY : nat, ityp, ntyp => nsp
  USE wavefunctions_module,  ONLY: psic
  USE uspp_param,            ONLY: upf, lmaxq, nh, nhm
  USE paw_variables,         ONLY : okpaw
  USE modes,                 ONLY : u
  USE qpoint,                ONLY : xq, eigqts
  USE gwus,                  ONLY : becsumort, alphasum
  USE units_gw,              ONLY : iudrhous, lrdrhous
  USE control_gw,            ONLY : lgamma
  USE noncollin_module,      ONLY : nspin_mag

  implicit none
  !
  !   the dummy variables
  !

  integer :: iflag, npe
  ! input: if zero does not compute drho
  ! input: the number of perturbations

  complex(DP) :: drhoscf (nrxx, nspin_mag), &
                      dbecsum (nhm*(nhm+1)/2, nat, nspin_mag)

  !inp/out: change of the charge density
  !input: sum over kv of bec
  integer ::  mode0
  ! input:the mode of the representation
  ! here the local variables

  integer :: ig, na, nt, ih, jh, mu, mode, ipert, is, ijh
  ! counter on G vectors
  ! counter on atoms
  ! counter on atomic type
  ! counter on beta functions
  ! counter on beta functions
  ! counter on r vectors
  ! pointer on modes
  ! pointer on the mode
  ! counter on perturbations
  ! counter on spin
  ! counter on combined beta functions

  real(DP), allocatable  :: qmod (:), qpg (:,:), ylmk0 (:,:)
  ! the modulus of q+G
  ! the values of q+G
  ! the spherical harmonics

  complex(DP) :: fact, zsum, bb, alpha, alpha_0, u1, u2, u3
  ! auxiliary variables
  !HL PH
  !complex(DP), allocatable ::  sk (:), qgm(:), drhous (:,:), aux (:,:,:)
  complex(DP), allocatable ::  sk (:), qgm(:), drhous (:,:), aux (:,:)
  ! the structure factor
  ! q_lm(G)
  ! contain the charge of drho
  ! auxiliary variable for drho(G)

  if (.not.okvan) return

  call start_clock ('addusddens')
  allocate (aux(  ngm , nspin_mag))    
  allocate (sk (  ngm))    
  allocate (ylmk0(ngm , lmaxq * lmaxq))    
  allocate (qgm(  ngm))    
  allocate (qmod( ngm))    
  if (.not.lgamma) allocate (qpg( 3  , ngm))    
  !      WRITE( stdout,*) aux, ylmk0, qmod
  !
  !  And then we compute the additional charge in reciprocal space
  !
  if (.not.lgamma) then
     call setqmod (ngm, xq, g, qmod, qpg)
     call ylmr2 (lmaxq * lmaxq, ngm, qpg, qmod, ylmk0)
     do ig = 1, ngm
        qmod (ig) = sqrt (qmod (ig) )
     enddo
  else
     call ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
     do ig = 1, ngm
        qmod (ig) = sqrt (gg (ig) )
     enddo
  endif

  fact = cmplx (0.d0, - tpiba, kind=DP)
  aux(:,:) = (0.d0, 0.d0)

!HL freezing modes perts
  npe = 1
  ipert = 1

  do nt = 1, ntyp
     if (upf(nt)%tvanp  ) then
        ijh = 0
        do ih = 1, nh (nt)
           do jh = ih, nh (nt)
          !qvan2 computes the fourier transform of the Q functions places them in qgm 
              call qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
              ijh = ijh + 1
             do na = 1, nat
                 if (ityp (na).eq.nt) then
                    !
                    ! calculate the structure factor
                    !
                    do ig = 1, ngm
                       sk (ig) = eigts1 (ig1 (ig), na) * &
                                 eigts2 (ig2 (ig), na) * &
                                 eigts3 (ig3 (ig), na) * &
                                 eigqts (na) * qgm (ig)
                    enddo
                    !
                    !  And qgmq and becp and dbecq
                    !
                       do is = 1, nspin_mag
                          if (iflag==1) then
                             zsum = dbecsum (ijh, na, is)
                          else
                             zsum = 2.0_DP*dbecsum (ijh, na, is)
                          endif
                             call zaxpy (ngm, zsum, sk, 1, aux(1,is), 1)
                       enddo
                 endif
              enddo
           enddo
        enddo
     endif
  enddo
 !
 !     convert aux to real space and add to density
 !
     mu = 1
     do is = 1, nspin_mag
        psic(:) = (0.d0, 0.d0)
        do ig = 1, ngm
           psic (nl (ig) ) = aux (ig, is)
        enddo
        call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
        call daxpy (2*nrxx, 1.0_DP, psic, 1, drhoscf(1,is), 1)
     enddo

  if (.not.lgamma) deallocate (qpg)
  deallocate (qmod)
  deallocate (qgm)
  deallocate (ylmk0)
  deallocate (sk)
  deallocate (aux)
  call stop_clock ('addusddens')
  return
end subroutine addusddens
