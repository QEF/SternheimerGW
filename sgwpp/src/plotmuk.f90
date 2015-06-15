SUBROUTINE plotmuk()

  USE kinds,                ONLY : DP
  USE constants,            ONLY : pi, RYTOEV, e2, fpi, eps8
  USE klist,                ONLY : nks, nkstot, ngauss, degauss, xk, wk, nelec
  USE units_coulmat,        ONLY : iuncoulmat, lrcoulmat
  USE cell_base,            ONLY : tpiba2, omega, at, alat, bg
  USE mp_world,             ONLY : world_comm, mpime
  USE symm_base,            ONLY : nsym, s, time_reversal, t_rev, ftau, invs, nrot
  USE io_global,            ONLY : ionode, stdout

IMPLICIT NONE

REAL(DP) :: xk_loc(3)
INTEGER  :: ik, iq1, i, iq
INTEGER  :: imq, isq(48), nqstar, nqs
INTEGER  :: nsq(48)
REAL(DP) :: sxq(3,48), xqs(3,48)
REAL(DP) :: muk(nks)

!READ the mu(k) file:
if (ionode) CALL davcio (muk, lrcoulmat, iuncoulmat, 1, -1)
print*, nks
DO ik=1, nks
   xk_loc(:) = xk(:,ik)
   write(5000+mpime,'(7x, "k1", 3f12.7)') xk(:,ik)
   CALL star_q(xk(:,ik), at, bg, nsym, s, invs, nqs, sxq, isq, nsq, imq, .false. )
   DO iq = 1, nqs
      CALL cryst_to_cart(1, sxq(:,iq), at, -1)
   ENDDO
   WRITE( 5000+mpime, '(7x,4f14.9)') ((sxq(i,iq1), i=1,3), (muk(ik)), iq1=1,nqs)
ENDDO

END SUBROUTINE plotmuk

