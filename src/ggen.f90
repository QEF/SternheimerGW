  !
  !----------------------------------------------------------------
  subroutine ggen ( gcutm)
  !----------------------------------------------------------------
  !
  ! FG - adapted from PW/ggen.f90
  !
  !   This routine generates all the reciprocal lattice vectors
  !   contained in the sphere of radius gcutm. Furthermore it
  !   computes the indices nl which give the correspondence
  !   between the fft mesh points and the array of g vectors.
  !
  ! FG - first I use some temporary array with size ngmx,
  !      then I allocate those whihc will be used outside
  !      using the proper bound ngm
  !
  use parameters
  use gspace
  implicit none
  !
  ! intent in
  !
  real(kind=DP) :: gcutm
  !
  ! local variables
  !
  real(kind=DP) ::  t(3), tt, swap
  integer :: n1, n2, n3, i, j, k, ipol, ng, igl, iswap, indsw, ngmx 
  real(kind=DP), allocatable :: g2sort(:), gg(:)
  real(kind=DP), parameter :: eps8 = 1.0D-8
  real(kind=DP), allocatable :: g_tmp(:,:)
  !
  write(6,'(/4x,a)') repeat('-',67)
  write(6,'(4x,"Parameters for this run")') 
  write(6,'(4x,a/)') repeat('-',67)
  !
  ! bound to number of G-vectors
  !
  ngmx = (2*nr1+3)*(2*nr2+3)*(2*nr3+3)
  write(6,'(4x,"ngmx = ",i10)') ngmx
  allocate ( g2sort(ngmx), gg(ngmx), g_tmp(3, ngmx) )
  !
  g2sort(:) = 1.0d20
  gg(:) = gcutm + 1.d0
  !
  ! determine G-vectors within the cutoff
  !
  ngm = 0
  do i = -nr1-1, nr1+1
   do j = -nr2-1, nr2+1
    do k = -nr3-1, nr3+1
       !
       tt = 0.d0
       do ipol = 1, 3
          t (ipol) = i * bg (ipol, 1) + j * bg (ipol, 2) + k * bg (ipol, 3)
          tt = tt + t (ipol) * t (ipol)
       enddo
       if (tt <= gcutm) then
          ngm = ngm + 1
          if ( tt > eps8 ) then
             g2sort(ngm) = tt
          else
             g2sort(ngm) = 0.d0
          endif
          g_tmp (1:3, ngm) = t (1:3)
          gg (ngm) = tt
       end if
       !
    enddo
   enddo
  enddo
  !
  ! copy g_tmp into g with the proper dimensions
  !
  allocate ( igtongl(ngm), g(3, ngm), nl(ngm) )
  g (:, 1:ngm) = g_tmp (:, 1:ngm)
  deallocate ( g_tmp )
  !
  !   reorder the g's in order of increasing magnitude. On exit
  !   from hpsort g2sort is ordered, and nl contains the new order.
  !
  nl (1) = 0
  call hpsort_eps ( ngm, g2sort, nl, eps8 )
  deallocate ( g2sort )
  !
  !   reorder the g vectors, their norm, and nl
  !
  do ng = 1, ngm - 1
     indsw = nl (ng)
     do while (indsw.ne.ng)
        do ipol = 1, 3
           swap = g (ipol, indsw)
           g (ipol, indsw) = g (ipol, nl (indsw) )
           g (ipol, nl (indsw) ) = swap
        enddo
        swap = gg (indsw)
        gg (indsw) = gg (nl (indsw) )
        gg (nl (indsw) ) = swap
        iswap = nl (ng)
        nl (ng) = nl (indsw)
        nl (indsw) = iswap
        !
        indsw = nl (ng)
     enddo
  enddo
  !
  !     Now set nl with the correct fft correspondence
  !
  do ng = 1, ngm
     ! n1 is going to be i+1, folded to positive when <= 0
     n1 = nint (g (1, ng) * at (1, 1) + g (2, ng) * at (2, 1) + g (3, ng) * at (3, 1) ) + 1
     if (n1.lt.1) n1 = n1 + nr1
     ! n2 is going to be j+1, folded to positive when <= 0
     n2 = nint (g (1, ng) * at (1, 2) + g (2, ng) * at (2, 2) + g (3, ng) * at (3, 2) ) + 1
     if (n2.lt.1) n2 = n2 + nr2
     ! n3 is going to be k+1, folded to positive when <= 0
     n3 = nint (g (1, ng) * at (1, 3) + g (2, ng) * at (2, 3) + g (3, ng) * at (3, 3) ) + 1
     if (n3.lt.1) n3 = n3 + nr3
     !
     if (n1.le.nr1.and.n2.le.nr2.and.n3.le.nr3) then
       nl (ng) = n1 + (n2 - 1) * nr1 + (n3 - 1) * nr1 * nr2
     else
        call error('ggen','Mesh too small?',ng)
     endif
  enddo
  !
  ! calculate number of G shells: ngl
  !
  ! G vectors are grouped in shells with the same norm
  !
  ngl = 1
  igtongl (1) = 1
  do ng = 2, ngm
     if (gg (ng) > gg (ng - 1) + eps8) then
        ngl = ngl + 1
     endif
     igtongl (ng) = ngl
  enddo
  !
  allocate ( gl(ngl) )
  !
  gl (1) = gg (1)
  igl = 1
  do ng = 2, ngm
     if (gg (ng) > gg (ng - 1) + eps8) then
        igl = igl + 1
        gl (igl) = gg (ng)
     endif
  enddo
  !
  if (igl.ne.ngl) call error ('ggen', 'igl <> ngl', ngl)
  !
  !  size of sub-matrix used to initialize CG
  !  (I brutally use an Hamiltonian with a smaller cutoff)
  !
  do ng = 1, ngm
    if (gg(ng).le.ecut0/ecutwfc*gcutm) ngm0 = ng
  enddo  
  !
  write(6,'(4x,"ngm       = ",i5/4x,"ngm0      = ",i5)') ngm, ngm0
  !
  deallocate ( gg )
  !
  ! total number of real-space grid points
  !
  nr = nr1 * nr2 * nr3
  write(6,'(4x,"nr1  = ",i10)') nr1
  write(6,'(4x,"nr2  = ",i10)') nr2
  write(6,'(4x,"nr3  = ",i10)') nr3
  write(6,'(4x,"nr   = ",i10)') nr
  !
  end subroutine ggen
  !
