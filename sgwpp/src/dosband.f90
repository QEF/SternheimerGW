SUBROUTINE dosband()
  USE klist,                ONLY : nks, nkstot, ngauss, degauss, xk, wk, nelec
  USE wvfct,                ONLY : nbnd, et, npw, igk, npwx,  g2kin, ecutwfc
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : tpiba2, omega, at, alat, bg
  USE ener,                 ONLY : ef
  USE constants,            ONLY : pi, RYTOEV, e2, fpi, eps8
  USE control_coulmat,      ONlY : degaussfs, nbndmin, debye_e, do_lind, ngcoul, do_diag

  IMPLICIT NONE

  REAL(DP)    :: degaussw0, w0g1, w0g2
  real(kind=DP), external :: efermig, dos_ef, w0gauss, wgauss
  REAL(DP)    :: kcut
  REAL(DP)    :: dosnnp(2)
  REAL(DP) :: xkp_loc(3), xk_loc(3)
  INTEGER  :: ibnd, ik
  REAL(DP) :: enk 

  kcut = 0.33333
  dosnnp  = 0.0d0
  degaussw0 = degaussfs

  write(6,'("\sigma: ", f12.7, "\pi: ", f12.7)') dosnnp(1), dosnnp(2)
DO ik = 1, nks
  DO ibnd = 1, nbnd
     xk_loc(:)   = xk(:,ik)
     CALL cryst_to_cart(1, xk_loc(:), at, -1)
     enk = (et(ibnd, ik) - ef)
     w0g1 = w0gauss ( enk / degaussw0, 0) / degaussw0
     IF(sqrt(xk_loc(1)**2 + xk_loc(2)**2) .lt. kcut ) then
        dosnnp(1) = dosnnp(1) + w0g1*(wk(ik)/2.0)
     ELSE
        dosnnp(2) = dosnnp(2) + w0g1*(wk(ik)/2.0)
     ENDIF
  ENDDO
ENDDO

  write(6,*) 
  write(6,'("\sigma: ", f12.7, "\pi: ", f12.7)') dosnnp(1), dosnnp(2)
  write(6,*) 
  write(6,*) dosnnp(:)

END SUBROUTINE dosband
