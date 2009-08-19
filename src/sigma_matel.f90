  !
  !----------------------------------------------------------------
  subroutine sigma_matel (  ik0, vr, xk0, nw, w )
  !----------------------------------------------------------------
  ! 
  use parameters
  use constants
  use gspace
  use kspace
#ifdef __PARA
  USE para
  USE mp_global,  ONLY : nproc, mpime, nproc_pool, my_pool_id, me_pool
  USE mp, ONLY:  mp_barrier
#endif
  implicit none
  !
  integer :: ig, igp, nw, iw, ibnd, jbnd, ios, ipol, ik0
  real(dbl) :: xk0(3), kplusg(3), g2kin(ngm), et(nbnd), w(nw), w_ryd(nw)
  complex(dbl) :: vr(nr), evc(ngm,nbnd), sigma(ngms,ngms,nw), aux(ngms)
  complex(kind=DP) :: ZDOTC, sigma_band(nbnd_sig,nbnd_sig,nw)
  real(dbl) :: resig_diag(nw), imsig_diag(nw), et_qp(nbnd_sig), a_diag(nw)
  !
  call start_clock ('sigma_matel')
  !
  w_ryd = w/ryd2ev
  !
  write(stdout,'(/4x,"k0(",i3," ) = (",3f7.3," )")') ik0, (xk0 (ipol) , ipol = 1, 3)
  !
  ! the k-dependent kinetic energy in Ry
  ! [Eq. (14) of Ihm,Zunger,Cohen J Phys C 12, 4409 (1979)]
  !
  do ig = 1, ngm
    kplusg = xk0 + g(:,ig)
    g2kin ( ig ) = tpiba2 * dot_product ( kplusg, kplusg )
  enddo
  !
  call  eigenstates_all ( vr, g2kin, evc, et )
  !
  write(stdout,'(4x,"LDA eigenval (eV)",8(1x,f7.3))') et(1:nbnd_sig)*ryd2ev
  !
#ifdef __PARA
  ! only proc 0 reads from file and does the product
  ! (need some sort of parallelization here)
  if (me.eq.1.and.mypool.eq.1) then
#endif
  read ( iunsigma, rec = ik0, iostat = ios) sigma
  !
  ! matrix elements of the self-energy
  !
  do ibnd = 1, nbnd_sig 
   do jbnd = 1, nbnd_sig
    do iw = 1, nw
      !
      sigma_band (ibnd, jbnd, iw) = czero
      !
      do ig = 1, ngms
        aux = sigma (ig, 1:ngms, iw) 
        sigma_band (ibnd, jbnd, iw) = sigma_band (ibnd, jbnd, iw) + &
            evc (ig, ibnd) * ZDOTC (ngms, evc (1:ngms, jbnd), 1, aux, 1)
      enddo 
      !
    enddo
   enddo
  enddo
  !
  ! Now calculate the expectation value of the self-energy
  ! using the diagonal matrix elements
  ! NOTE: we cannot calculate the QP correction until we
  ! have the expt value of the Vxc. I could get Vxc from the
  ! density and the LDA expression.
  !
  do ibnd = 1, nbnd_sig
    !
    write(stdout,*)
    do iw = 1, nw
      resig_diag (iw) = real( sigma_band (ibnd, ibnd, iw) )
      imsig_diag (iw) = aimag ( sigma_band (ibnd, ibnd, iw) )
      a_diag (iw) = one/pi * abs ( imsig_diag (iw) ) / &
         ( abs ( w_ryd(iw) - et(ibnd) - resig_diag (iw) )**2.d0 &
          + abs ( imsig_diag (iw) )**2.d0 ) 
!     write(stdout,'(5f15.8)') et(ibnd)*ryd2ev, w_ryd(iw)*ryd2ev, &
!        resig_diag (iw)*ryd2ev, imsig_diag (iw)*ryd2ev, a_diag (iw)*ryd2ev
    enddo
    !
    call qp_eigval ( nw, w_ryd, resig_diag, et(ibnd), et_qp (ibnd) )  
    !
  enddo
  !
#ifdef __PARA
  endif
#endif
  write(stdout,'(4x,"GW  expt val (eV)",8(1x,f7.3)/)') et_qp*ryd2ev
  !
  call stop_clock ('sigma_matel')
  ! 
  return
  end subroutine sigma_matel
  !----------------------------------------------------------------
  subroutine qp_eigval ( nw, w, sig, et, et_qp )  
  !----------------------------------------------------------------
  !
  use parameters
  use constants
  implicit none
  integer :: nw, iw, iw1, iw2
  real(DP) :: w(nw), sig(nw), et, et_qp, dw, w1, w2, et_qp0, sig1, sig2
  !
  dw = w(2)-w(1)
  !
  if ((et.lt.w(1)+dw).or.(et.gt.w(nw)-dw)) &
    call errore ('qp_eigval','original eigenvalues outside the frequency range of the self-energy',1)
  iw = 1
  do while ((iw.lt.nw).and.(w(iw).lt.et))
    iw = iw + 1
    iw1 = iw-1
    iw2 = iw
  enddo
  w1 = w(iw1)
  w2 = w(iw2)
  sig1 = sig(iw1)
  sig2 = sig(iw2)
  !
  et_qp0 = sig1 + ( sig2 - sig1 ) * (et-w1) / (w2-w1)
  !
  ! temporary - until I do not have Vxc
  !
  et_qp = et_qp0
  !
  end subroutine qp_eigval
  !----------------------------------------------------------------
  !







