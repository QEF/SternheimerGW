  !-----------------------------------------------------------------------
  subroutine rw_haydock ( ik, ig1, ig2, a, b, a_term, b_term, norm, action)
  !-----------------------------------------------------------------------
  !
  use constants
  use parameters
  use gspace, only : ngms
  implicit none
  integer :: action
  real(DP) :: a(nstep+2,4), b(nstep+2,4), norm(4), a_term(4), b_term(4)
  real(DP) :: coeff_pack(8*nstep+28)
  integer :: i, ik, ig1, ig2, istep, ih, ios, rec0
  !
  rec0 = ig1 + (ig2-1)*ngms + (ik-1)*ngms*ngms
  !
  if (action.eq.1) then
    !
    i = 0
    do ih = 1, 4
      do istep = 1, nstep+2
         i = i + 1
         coeff_pack (i) = a(istep,ih)
      enddo
    enddo
    do ih = 1, 4
      do istep = 1, nstep+2
         i = i + 1
         coeff_pack (i) = b(istep,ih)
      enddo
    enddo
    do ih = 1, 4
      i = i + 1
      coeff_pack (i) = a_term(ih)
    enddo
    do ih = 1, 4
      i = i + 1
      coeff_pack (i) = b_term(ih)
    enddo
    do ih = 1, 4
      i = i + 1
      coeff_pack (i) = norm(ih)
    enddo
    !
    write ( iuncoeff, rec = rec0, iostat = ios) coeff_pack
    !
  elseif (action.eq.-1) then
    !
    read ( iuncoeff, rec = rec0, iostat = ios) coeff_pack
    !
    i = 0
    do ih = 1, 4
      do istep = 1, nstep+2
         i = i + 1
         a(istep,ih) = coeff_pack (i) 
      enddo
    enddo
    do ih = 1, 4
      do istep = 1, nstep+2
         i = i + 1
         b(istep,ih) = coeff_pack (i)
      enddo
    enddo
    do ih = 1, 4
      i = i + 1
      a_term(ih) = coeff_pack (i) 
    enddo
    do ih = 1, 4
      i = i + 1
      b_term(ih) = coeff_pack (i) 
    enddo
    do ih = 1, 4
      i = i + 1
      norm(ih) = coeff_pack (i) 
    enddo
    !
  else
    call error ('rw_haydock','wrong action specified',action)
  endif
  !
  return
  end subroutine rw_haydock
  !

