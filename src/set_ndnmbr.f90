  !
  !----------------------------------------------------------
  subroutine set_ndnmbr ( pool, proc, procp, npool, ndlab)
  !----------------------------------------------------------
  !
  !  create ndlab label from pool and proc numbers
  !
  !  The rule for deciding the node number is based on 
  !  the restriction set in startup.f90 that every pool 
  !  has the same number of procs.
  !
  !----------------------------------------------------------
  !
  implicit none
  character(len=3) :: ndlab
  integer :: pool, proc, procp, node, nprocs, npool 
  ! pool = 1,...,npool
  ! proc = 0,...,nproc_pool-1
  !
  nprocs = npool * procp
  !
  node = (pool-1)*procp + proc + 1
  !
  ndlab = '   '
  if ( nprocs < 10 ) then
    write (ndlab(1:1),'(i1)') node
  elseif ( nprocs < 100 ) then
    if ( node < 10 ) then
       ndlab = '0'
       write (ndlab(2:2),'(i1)') node
    else
       write (ndlab(1:2),'(i2)') node
    endif
  else
    if ( node < 10 ) then
       ndlab = '00'
       write (ndlab(3:3),'(i1)') node
    elseif ( node < 100 ) then
       ndlab = '0'
       write (ndlab(2:3),'(i2)') node
    else
       write (ndlab,'(i3)') node
    endif
  endif
  !
  end subroutine set_ndnmbr
  !----------------------------------------------------------
  !
