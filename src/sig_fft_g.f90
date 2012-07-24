  !
  !----------------------------------------------------------------
  SUBROUTINE sig_fft_g (nr1tmp, nr2tmp, nr3tmp, nrtmp, ecuttmp, corx)
  !----------------------------------------------------------------
  !----------------------------------------------------------------
  ! FIND THE G-VECTORS FOR THE SIGMA CUTOFF
  !----------------------------------------------------------------
  ! the G^2 cutoff in units of 2pi/a_0
  ! Note that in Ry units the kinetic energy is G^2, not G^2/2
  ! FG: determine G-vectors within the cutoff from the
  ! array already created in ggen
  !
  ! HL: This routine sets two Gamma centered G-grids and the corresponding size of the real 
  ! space grid for Sigma_Correlation, and Sigma_exchange. The exchange cutoff of the self-energy
  ! needs to be high enough that the wave function*barecoulomb product doesn't get brutalized. 

  ! corx determines whether it is the correlation grid being generated or the 
  ! exchange within this routine. 

  ! this switch is necessary for the time being. Alternatively I could define this routine as 
  ! being contained in the gw_sigma module, and then I could directly pass the arrays from main
  ! to the sub routine to be allocated and populated. 


  USE kinds,            ONLY : DP
  USE constants,        ONLY : tpi
  USE gvect,            ONLY : gcutm, ecutwfc, dual, nr1, nr2, nr3, ngm, g, igtongl, gl, nl
  USE cell_base,        ONLY : at, tpiba2
  USE fft_scalar,       ONLY : allowed
  USE basis,            ONLY : starting_wfc, starting_pot, startingconfig
  USE qpoint,           ONLY : xq
  USE control_gw,       ONLY : done_bands, reduce_io, recover, tmp_dir_gw, &
                              ext_restart, bands_computed
  USE save_gw,          ONLY : tmp_dir_save
  USE input_parameters, ONLY : pseudo_dir
  USE io_files,         ONLY : prefix, tmp_dir
  USE control_flags,    ONLY : restart
  USE gwsigma,          ONLY : nlsex, nlsco, ngmsex, ngmsco, ngmsig
  
  IMPLICIT NONE

  integer  :: n1, n2, n3, i, j, k, ipol, ig, igl, ng
  integer  :: corx
  REAL(DP) :: ecuttmp
  REAL(DP) :: gcuttmp
  INTEGER  :: nr1tmp, nr2tmp, nr3tmp, nrtmp
  INTEGER  :: ngmtmp
  INTEGER, ALLOCATABLE :: nltmp(:)
  !INTEGER, POINTER :: nltmp(:)

  gcuttmp = 4.D0 * ecuttmp / tpiba2

!HLS 
! gcuttmp = ecuttmp
! gcutsigtmp = 4.D0 * ecuttmp / tpiba2
! nr1tmp = 1 + int (2 * sqrt (gcutsigtmp) * sqrt( at(1,1)**2 + at(2,1)**2 + at(3,1)**2 ) )
! nr2tmp = 1 + int (2 * sqrt (gcutsigtmp) * sqrt( at(1,2)**2 + at(2,2)**2 + at(3,2)**2 ) )
! nr3tmp = 1 + int (2 * sqrt (gcutsigtmp) * sqrt( at(1,3)**2 + at(2,3)**2 + at(3,3)**2 ) )

  nr1tmp = 1 + int (2 * sqrt (gcuttmp) * sqrt( at(1,1)**2 + at(2,1)**2 + at(3,1)**2 ) )
  nr2tmp = 1 + int (2 * sqrt (gcuttmp) * sqrt( at(1,2)**2 + at(2,2)**2 + at(3,2)**2 ) )
  nr3tmp = 1 + int (2 * sqrt (gcuttmp) * sqrt( at(1,3)**2 + at(2,3)**2 + at(3,3)**2 ) )



  do while (.not.allowed(nr1tmp))
    nr1tmp = nr1tmp + 1
  enddo

  do while (.not.allowed(nr2tmp))
    nr2tmp = nr2tmp + 1
  enddo

  do while (.not.allowed(nr3tmp))
    nr3tmp = nr3tmp + 1
  enddo

  !Loop over all g vectors in (order of size) and test that |G|.le.ngmtmp.
  !igtongl should only be defined when init_run is called in pwscf.

  do ng = 1, ngm
    if ( gl( igtongl (ng) ) .le. gcuttmp ) ngmtmp = ng
    !this determines cut off of G, W:
    !if ( gl( igtongl (ng) ) .le. gcuttmp ) ngmwcutoff = ng
  enddo

  !do ng = 1, ngm
  !   this determines cut off of \Sigma and hence the FFT grid:
  !   if ( gl( igtongl (ng) ) .le. gcutsigtmp ) ngmwcutoff = ng
  !enddo

!Choose whether is is Ex or Corr grid we are generating.
  ALLOCATE (nltmp(ngmtmp))

  if(corx.eq.1) then
     ALLOCATE (nlsco(ngmtmp))
    !ngmsco should be 4*ngmtmp
     ngmsco = ngmtmp
     ngmsig = ngmtmp
    else
     ALLOCATE (nlsex(ngmtmp))
     ngmsex = ngmtmp
  endif 
 !
 ! Now set nl with the correct fft correspondence
 !
  do ng = 1, ngmtmp
     ! n1 is going to be i+1, folded to positive when <= 0
     n1 = nint (g (1, ng) * at (1, 1) + g (2, ng) * at (2, 1) + g (3, ng) * at (3, 1) ) + 1
     if (n1.lt.1) n1 = n1 + nr1tmp

     ! n2 is going to be j+1, folded to positive when <= 0
     n2 = nint (g (1, ng) * at (1, 2) + g (2, ng) * at (2, 2) + g (3, ng) * at (3, 2) ) + 1
     if (n2.lt.1) n2 = n2 + nr2tmp

     ! n3 is going to be k+1, folded to positive when <= 0
     n3 = nint (g (1, ng) * at (1, 3) + g (2, ng) * at (2, 3) + g (3, ng) * at (3, 3) ) + 1
     if (n3.lt.1) n3 = n3 + nr3tmp

     !
     if (n1.le.nr1tmp.and.n2.le.nr2tmp.and.n3.le.nr3tmp) then
       nltmp (ng) = n1 + (n2 - 1) * nr1tmp + (n3 - 1) * nr1tmp * nr2tmp
     else
       !call error('ggens','Mesh too small?',ng)
        WRITE(6,'("ERROR: Mesh too small?")')
        STOP
     endif
  enddo

  nrtmp = nr1tmp * nr2tmp * nr3tmp 

  if (corx.eq.1) then
     nlsco = nltmp
     write(6,'(4x,"Correlation cutoffs: ")')
     write(6,'(5x,"ecutsco = ",f10.6)') ecuttmp
     write(6,'(5x,"ngmsco = ",i10)') ngmtmp
     write(6,'(5x,"nr1sco = ",i10)') nr1tmp
     write(6,'(5x,"nr2sco = ",i10)') nr2tmp
     write(6,'(5x,"nr3sco = ",i10)') nr3tmp
     write(6,'(5x,"nrsco  = ",i10)') nrtmp
     write(6,'(5x,"")') 
  
  else
     nlsex = nltmp
     write(6,'(4x,"Exchange cutoffs: ")')
     write(6,'(5x,"ecutsex = ",f10.7)') ecuttmp
     write(6,'(5x,"ngmsex = ",i10)') ngmtmp
     write(6,'(5x,"nr1sex = ",i10)') nr1tmp
     write(6,'(5x,"nr2sex = ",i10)') nr2tmp
     write(6,'(5x,"nr3sex = ",i10)') nr3tmp
     write(6,'(5x,"nrsex  = ",i10)') nrtmp
  endif

  DEALLOCATE(nltmp)
  END SUBROUTINE sig_fft_g
