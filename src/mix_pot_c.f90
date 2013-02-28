  !
  !-----------------------------------------------------------------------
  subroutine mix_potential_c (ndim, vout, vin, alphamix, dr2, tr2, &
     iter, n_iter, conv)
  !-----------------------------------------------------------------------
  !
  ! slightly modified (simplified) for gwdfpt code - FG
  ! and adapted to complex potentials
  ! -- to make it complex I just considered that the beta matrix must
  ! be hermitian (looking at the paper by Johnson), the method seems
  ! to work like a rocket !! - need to investigate more...
  !
  ! ----------------------------------------------------------------------
  !
  ! Modified Broyden's method for potential/charge density mixing
  !             D.D.Johnson, PRB 38, 12807 (1988)
  ! On input :
  !    ndim      dimension of arrays vout, vin
  !    vout      output potential/rho at current iteration
  !    vin       potential/rho at previous iteration
  !    alphamix  mixing factor (0 < alphamix <= 1)
  !    tr2       threshold for selfconsistency
  !    iter      current iteration number
  !    n_iter    number of iterations used in the mixing
  !    filename  if present save previous iterations on file 'filename'
  !              otherwise keep everything in memory
  ! On output:
  !    dr2       [(vout-vin)/ndim]^2
  !    vin       mixed potential
  !    vout      vout-vin
  !    conv      true if dr2.le.tr2
  !
  ! ----------------------------------------------------------------------
  !
  !USE parameters, only : DP, maxter => nmax_iter 

  USE kinds,       ONLY : DP 
  USE control_gw,  ONLY : nmix_gw
  USE mp_global,   ONLY : inter_pool_comm, intra_pool_comm, mpime, mp_global_end,&
                          intra_image_comm


  !max number of iterations used in mixing: n_iter must be.le.maxter
  implicit none

  !
  !   First the dummy variables
  !
  integer :: ndim, iter, n_iter
  complex(kind=DP) :: vout (ndim), vin (ndim)
  real(kind=DP) :: alphamix, dr2, tr2
  logical :: conv
  !
  !   Here the local variables
  !

  integer :: maxter
  parameter (maxter = 8)

  integer :: iunit, iunmix, n, i, j, iwork (maxter), info, iter_used, &
       ipos, inext, ndimtot
  ! work space containing info from previous iterations:
  ! must be kept in memory and saved between calls 

  complex(kind=DP), allocatable, save :: df (:,:), dv (:,:)
  complex(kind=DP), allocatable :: vinsave (:)
  real(kind=DP) :: norm
  complex(kind=DP) :: gamma, beta (maxter, maxter), work (maxter)
  real(kind=DP) :: DZNRM2
  complex(kind=DP) :: ZDOTC
  external ZDOTC, DZNRM2
  !adjustable parameters as suggested in the original paper
  real(kind=DP) w (maxter), w0
  !data w0 / 0.01d0 /, w / maxter * 1.d0 /
!HL Hard coding the number of iterations
  data w0 /0.01d0/, w/maxter*1.d0/
  !
 
  !if(mpime.eq.0) write(6,*)w

  do n = 1, ndim
     vout (n) = vout (n) - vin (n)
  enddo
  dr2 = DZNRM2 (ndim, vout, 1) **2
  ndimtot = ndim
  dr2 = (sqrt (dr2) / ndimtot) **2
  !
  conv = dr2.lt.tr2
  !
  if (iter.eq.1) then
     if (.not.allocated(df)) allocate (df( ndim , n_iter))
     if (.not.allocated(dv)) allocate (dv( ndim , n_iter))
  endif
  if (conv) then
     deallocate (dv)
     deallocate (df)
     return
  endif
  allocate (vinsave( ndim))
  !
  ! iter_used = iter-1  if iter <= n_iter
  ! iter_used = n_iter  if iter >  n_iter
  !
  iter_used = min (iter - 1, n_iter)
  !
  ! ipos is the position in which results from the present iteraction
  ! are stored. ipos=iter-1 until ipos=n_iter, then back to 1,2,...
  !
  ipos = iter - 1 - ( (iter - 2) / n_iter) * n_iter
  !
  if (iter.gt.1) then
     do n = 1, ndim
        df (n, ipos) = vout (n) - df (n, ipos)
        dv (n, ipos) = vin (n) - dv (n, ipos)
     enddo
     norm = (DZNRM2 (ndim, df (1, ipos), 1) ) **2
     norm = sqrt (norm)
     call ZSCAL (ndim, dcmplx ( 1.d0 / norm, 0.d0), df (1, ipos), 1)
     call ZSCAL (ndim, dcmplx ( 1.d0 / norm, 0.d0), dv (1, ipos), 1)
  endif
  !
  call ZCOPY (ndim, vin, 1, vinsave, 1)
  !
  do i = 1, iter_used
     do j = i + 1, iter_used
        beta (i, j) = w (i) * w (j)  * ZDOTC (ndim, df (1, j), 1, df (1, i), 1)
     enddo
     beta (i, i) = w0**2 + w (i) **2
  enddo
  !
 !The factorisation seems to go screwy here if we enforce the hermiticity
 !going to try it with a generalized inversion routine in case the matrix is no longer hermitian.
  call ZHETRF ('U', iter_used, beta, maxter, iwork, work, maxter, info)
 !HL generalized for non-hermitian matrix:
 !call ZGETRF (iter_used, iter_used, beta, maxter, iwork, info)
  call errore ('broyden', 'factorization', info)
 !HL
  call ZHETRI ('U', iter_used, beta, maxter, iwork, work, info)
 !Generalized back substituion for non-hermitian matrix:
 !call ZGETRI (iter_used, beta, maxter, iwork, work, maxter, info)
  call errore ('broyden', 'ZSYTRI', info)

!HL Trying with the hermiticity condition led to disaster and tears
 do i = 1, iter_used
    do j = i + 1, iter_used
       beta (j, i) = conjg ( beta (i, j) ) 
    enddo
 enddo
!
  do i = 1, iter_used
     work (i) = ZDOTC (ndim, df (1, i), 1, vout, 1)
  enddo
!
  do n = 1, ndim
     vin (n) = vin (n) + alphamix * vout (n)
  enddo
!
  do i = 1, iter_used
     gamma = 0.d0
     do j = 1, iter_used
        gamma = gamma + beta (j, i) * w (j) * work (j)
     enddo
!
     do n = 1, ndim
        vin (n) = vin (n) - w (i) * gamma * (alphamix * df (n, i) + dv (n, i) )
     enddo
  enddo
!
  inext = iter - ( (iter - 1) / n_iter) * n_iter
  call ZCOPY (ndim, vout, 1, df (1, inext), 1)
  call ZCOPY (ndim, vinsave, 1, dv (1, inext), 1)
  deallocate(vinsave)
  !
  return
  end subroutine mix_potential_c
  !-----------------------------------------------------------------------
