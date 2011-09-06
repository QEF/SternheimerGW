SUBROUTINE sigma_matel (ik0)
  !
  !----------------------------------------------------------------
  !
  ! Known bug: below (grep @) the code works fine with the convention
  ! <i|Sigma|j> = sum_G,G' [u_i,k(G)]* <G|Sigma|G'> u_j,k(G')
  ! which is NOT the convention set in the paper. The convention of
  ! the paper should be:
  ! <i|Sigma|j> = sum_G,G' u_i,-k(G) <G|Sigma|G'> [u_j,-k(G')]*
  ! but if I use -xk0 below the code gives wrong results. Indeed the
  ! shape of Sigma^c is slightly wrong, but most importantly the
  ! Sigma^ex is totally wrong. A quick sanity check is to calculate
  ! the sandwiches of Sigma^ex with v set to the delta function -
  ! this should give the normalization of the wfs. With -xk0 this
  ! normalization is screwed up, while everything works fine with
  ! +xk0. Also the QP energies are ok when using +xk0.
  ! In particolar, with -xk0 the Gamma point at 1 1 1 does not give
  ! the same Sigma^ex as 0 0 0, while this is the case for +xk0.
  ! The most obvious conclusion is that I messed up somewhere in the
  ! analytical calculation of the matrix elements - need to check this
  ! out once more.
  !
  ! Note that in any case Sigma(-k,-G,-G') = Sigma(k,G,G') for silicon
  ! if we exploit inversion symmetry.
  !
  !----------------------------------------------------------------

  USE io_global,    ONLY : stdout
  USE io_files,             ONLY : prefix, iunigk
  USE kinds,         ONLY : DP
  USE gvect,         ONLY : ngm, nrxx, g, nr1, nr2, nr3, nrx1, nrx2, nrx3, nl
  USE gsmooth,       ONLY : nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nls, ngms
  USE constants,     ONLY : e2, fpi, RYTOEV, tpi, pi
  USE freq_gw,       ONLY : fpol, fiu, nfs, nwsigma, wsigma
  USE klist,         ONLY : xk, wk, nkstot
  USE wvfct,         ONLY : nbnd, npw, npwx, igk, g2kin, et
  USE scf,           ONLY : rho
  USE qpoint,        ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE units_gw,      ONLY : iunsigma, iuwfc, lrwfc, lrsigma
  USE control_gw,    ONLY : nbnd_occ
  USE gwsigma,      ONLY : ngmsig, nbnd_sig   
  USE wavefunctions_module,          ONLY : evc


IMPLICIT NONE

INTEGER    ::   ig, igp, nw, iw, ibnd, jbnd, ios, ipol, ik0, ir
REAL(DP)   ::   w_ryd(nwsigma)
REAL(DP)   :: resig_diag(nwsigma,nbnd_sig), imsig_diag(nwsigma,nbnd_sig), et_qp(nbnd_sig), a_diag(nwsigma,nbnd_sig)
REAL(DP)     ::   dresig_diag(nwsigma,nbnd_sig), vxc_tr, vxc_diag(nbnd_sig)
REAL(DP)     ::   resig_diag_tr(nwsigma), imsig_diag_tr(nwsigma), a_diag_tr(nwsigma), et_qp_tr, z_tr, z(nbnd_sig)
REAL(DP)     ::  v_xc(nrxx), one
COMPLEX(DP)  ::   czero
COMPLEX(DP)  ::   aux(ngmsig), sigma(ngmsig,ngmsig,nwsigma), psic(nrxx), vpsi(ngm)
COMPLEX(DP)  ::   ZDOTC, sigma_band(nbnd_sig,nbnd_sig,nwsigma), vxc(nbnd_sig,nbnd_sig)
LOGICAL      :: do_band, do_iq, setup_pw, exst
INTEGER      :: iman, nman, ndeg(nbnd_sig), ideg, iq 


one   = 1.0d0 
czero = (0.0d0, 0.0d0)
w_ryd = wsigma/RYTOEV

write(stdout,'(/4x,"k0(",i3," ) = (",3f7.3," )")') ik0, (xk (ipol,ik0) , ipol = 1, 3)

  !@FG
  ! NOTE - I calculate the eigenstates of -xk0
  ! in order to have c_k(-G) = [c_-k(G)]*
  ! Because of my convention on the FFTs on G,G' in the paper,
  ! below we mix G and -G in the sandwitches. The easiest way
  ! to perform the calculation is to use the eigenvectors for -xk0
  ! and take their cc to obtain c(-G) for xk0
  !  note the -xk0 for the reason above!
  !@ kplusg = -xk0 + g(:,ig)

  !nbnd for silicon looking at 4 occupied valence states and four un-occupied.

     nbnd = nbnd_sig 

   !do_band = .true.

      iq = 1 
      CALL prepare_kmq(do_band, do_iq, setup_pw, iq, ik0)
      CALL run_pwscf(do_band)
      CALL initialize_gw()
  
  !READ wave function at \psi_{ik0}

     if (nksq.gt.1) rewind (unit = iunigk)

        if (nksq.gt.1) then
           read (iunigk, err = 100, iostat = ios) npw, igk
 100        call errore ('green_linsys', 'reading igk', abs (ios) )
        endif

      CALL davcio (evc, lrwfc, iuwfc, 1, -1)

! Set zero of energy to top of the valence band.
! SHIFT
! et(:,:) = et(:,:) - (6.2/13.605)


! #ifdef __PARA
! only proc 0 reads from file and does the product
! (need some sort of parallelization here)
! if (me.eq.1.and.mypool.eq.1) then
! #endif

! MATRIX ELEMENTS OF THE XC POTENTIAL
! vxc is a just a long list of the potential at different points on the grid
! I think what I will do is generate it in a one shot calc 
! from pw.x and just leave it at that. 

! Should the potential from bigger more complicated systems be required 
! we will treat them on a case by case basis.
! Need to figure out where the above @ bug comes from?

     open(unit=110,file='vxc.dat')
     rewind(110)
     do ir = 1, nrxx
       read(110,*) v_xc(ir)
     enddo
    close(110)

!Following the same convention for fourier transforming wave functions as appears in incdrhoscf.f90

  WRITE(6,'("NBND")')
  WRITE(6,*) nbnd_sig
!
  do jbnd = 1, nbnd_sig

    psic = czero

    do ig = 1, npw
!SGW psic ( nl (ig) ) = evc(ig, jbnd)
!     psic ( nls (igk(ig)) ) = conjg(evc(ig, jbnd))
      psic ( nls (igk(ig)) ) = evc(ig, jbnd)
    enddo

!   call cfft3 ( psic, nr1, nr2, nr3,  2)
!   Remember for wave functions we use the index \pm 2.

    call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)

    do ir = 1, nrxx
      psic (ir) = psic(ir) * v_xc (ir)
  !   evxc = (omega/1728) * psic(ir) * v_xc(ir) * conjg(psic(ir)
    enddo
!   call cfft3 ( psic, nr1, nr2, nr3, -1)
    call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2)

    do ig = 1, npw
      !vpsi(ig) = psic( nl(ig) )
       vpsi(ig) = psic(nls(igk(ig)))
    enddo
    !
    do ibnd = 1, nbnd_sig
       vxc(ibnd,jbnd) = ZDOTC (npw, evc (1, ibnd), 1, vpsi, 1)
       !write(6,*) vxc(ibnd,jbnd)
    enddo
    !
  enddo
!
!  MATRIX ELEMENTS OF THE SELF-ENERGY
!
  CALL davcio (sigma, lrsigma, iunsigma, ik0, -1)

!Need nbnd_sig = 8 i.e. for Si we want to look at the 4 valence states and 4 conduction states.
  do ibnd = 1, nbnd_sig
   do jbnd = 1, nbnd_sig
    do iw = 1, nwsigma
      !   
        sigma_band (ibnd, jbnd, iw) = czero
      !   
        do ig = 1, ngmsig

        !do ig = 1, npw
        !aux = sigma (ig, 1:ngms, iw) ! paper convention

            aux = sigma (1:ngmsig, ig, iw)

!why not conjg evc(ig,ibnd) ?

!            sigma_band (ibnd, jbnd, iw) = sigma_band (ibnd, jbnd, iw) + &
!                evc (ig, ibnd) * ZDOTC (ngmsig, evc (1:ngmsig, jbnd), 1, aux, 1)

            sigma_band (ibnd, jbnd, iw) = sigma_band (ibnd, jbnd, iw) + &
                 evc (ig, ibnd) * ZDOTC (ngmsig, evc (1:ngmsig, jbnd), 1, aux, 1)

!                evc (ig, ibnd) * ZDOTC (ngmsig, evc (1:npw, jbnd), 1, aux, 1)
!                evc (ig, ibnd) * ZDOTC (ngmsig, evc (1, jbnd), 1, aux, 1)

        enddo
      enddo
     enddo
    enddo

  ! Now calculate the expectation value of the self-energy
  ! using the diagonal matrix elements
  !
   do ibnd = 1, nbnd_sig
     !
     do iw = 1, nwsigma
       resig_diag (iw,ibnd) = real( sigma_band (ibnd, ibnd, iw) )
       dresig_diag (iw,ibnd) = resig_diag (iw,ibnd) - real( vxc(ibnd,ibnd) )
       imsig_diag (iw,ibnd) = aimag ( sigma_band (ibnd, ibnd, iw) )
       a_diag (iw,ibnd) = one/pi * abs ( imsig_diag (iw,ibnd) ) / &
          ( abs ( w_ryd(iw) - et(ibnd, ik0) - ( resig_diag (iw,ibnd) - vxc(ibnd,ibnd) ) )**2.d0 &
           + abs ( imsig_diag (iw,ibnd) )**2.d0 )
     enddo
     !
     call qp_eigval ( nwsigma, w_ryd, dresig_diag(1,ibnd), et(ibnd,ik0), et_qp (ibnd), z(ibnd) )
     !
    enddo
  !
  ! Now take the trace (get rid of phase arbitrariness of the wfs)
  ! (alternative and more approrpiate: calculate the nondiag on the
  ! deg subspaces and diagonalize)
  !
  ! count degenerate manifolds and degeneracy...
  !
   nman = 1
   ndeg = 1
   do ibnd = 2, nbnd_sig
     if ( abs( et (ibnd, ik0) - et (ibnd-1, ik0)  ) .lt. 1.d-5 ) then
       ndeg (nman) = ndeg(nman) + 1
     else
       nman = nman + 1
     endif
   enddo

   write(6,'(" Manifolds")')
   write (stdout, *) nman, (ndeg (iman) ,iman=1,nman)
  !
  ! ...and take the trace over the manifold
  !

  ibnd = 0
  jbnd = 0
  do iman = 1, nman
    !
    resig_diag_tr = 0.d0
    imsig_diag_tr = 0.d0
    a_diag_tr = 0.d0
    et_qp_tr = 0.d0
    z_tr = 0.d0
    vxc_tr = 0.d0
    !
    do ideg = 1, ndeg(iman)
      ibnd = ibnd + 1
      resig_diag_tr = resig_diag_tr + resig_diag (:,ibnd)
      imsig_diag_tr = imsig_diag_tr + imsig_diag (:,ibnd)
      a_diag_tr = a_diag_tr + a_diag (:,ibnd)
      et_qp_tr = et_qp_tr + et_qp (ibnd)
      z_tr = z_tr + z (ibnd)
      vxc_tr = vxc_tr + real(vxc(ibnd,ibnd))
    enddo
   !
    do ideg = 1, ndeg(iman)
      jbnd = jbnd + 1
      resig_diag (:,jbnd) = resig_diag_tr / float( ndeg(iman) )
      imsig_diag (:,jbnd) = imsig_diag_tr / float( ndeg(iman) )
      a_diag (:,jbnd) = a_diag_tr / float( ndeg(iman) )
      et_qp (jbnd) = et_qp_tr / float( ndeg(iman) )
      z (jbnd) = z_tr / float( ndeg(iman) )
      vxc_diag (jbnd) = vxc_tr / float( ndeg(iman) )
    enddo
   !
  enddo
   !
  write(stdout,'(/4x,"LDA eigenval (eV)",8(1x,f7.3))') et(1:nbnd_sig, ik0)*RYTOEV
  write(stdout,'(4x,"Vxc expt val (eV)",8(1x,f7.3))') vxc_diag(1:nbnd_sig)*RYTOEV
  write(stdout,'(4x,"GW qp energy (eV)",8(1x,f7.3))') et_qp(1:nbnd_sig)*RYTOEV
  write(stdout,'(4x,"GW qp renorm     ",8(1x,f7.3)/)') z(1:nbnd_sig)
  !
  do iw = 1, nwsigma
    write(stdout,'(9f15.8)') wsigma(iw), (RYTOEV*resig_diag (iw,ibnd), ibnd=1,nbnd_sig)
  enddo
  write(stdout,*)
  do iw = 1, nwsigma
    write(stdout,'(9f15.8)') wsigma(iw), (RYTOEV*imsig_diag (iw,ibnd), ibnd=1,nbnd_sig)
  enddo
  write(stdout,*)
  do iw = 1, nwsigma
    write(stdout,'(9f15.8)') wsigma(iw), (a_diag (iw,ibnd)/RYTOEV, ibnd=1,nbnd_sig)
  enddo
!
! #ifdef __PARA
!   endif
! #endif

    CALL clean_pw_gw(ik0)
RETURN
END SUBROUTINE sigma_matel


!----------------------------------------------------------------
  SUBROUTINE  qp_eigval ( nw, w, sig, et, et_qp, z )
!----------------------------------------------------------------
!
  USE kinds,         ONLY : DP

  IMPLICIT NONE

  integer :: nw, iw, iw1, iw2
  real(DP) :: w(nw), sig(nw), et, et_qp, dw, w1, w2, sig_et, sig1, sig2, z, sig_der, one
  
  one = 1.0d0
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
  sig_et = sig1 + ( sig2 - sig1 ) * (et-w1) / (w2-w1)
!
  sig_der = ( sig2 - sig1 ) / ( w2 - w1 )
  z = one / ( one - sig_der)
!
! temporary - until I do not have Vxc
!
  et_qp = et + z * sig_et
!
  END SUBROUTINE qp_eigval
!----------------------------------------------------------------
!
