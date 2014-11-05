SUBROUTINE sigma_exch(ik0)
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE disp,          ONLY : nqs, nq1, nq2, nq3, wq, x_q, xk_kpoints, num_k_pts
  USE control_gw,    ONLY : eta
  USE klist,         ONLY : wk, xk, nkstot, nks
  USE io_files,      ONLY : prefix, iunigk
  USE wvfct,         ONLY : nbnd, npw, npwx, igk, g2kin, et
  USE cell_base,     ONLY : omega, tpiba2, at, bg, tpiba, alat
  USE eqv,           ONLY : evq, eprec
  USE units_gw,      ONLY : iunsex, lrsex, lrwfc, iuwfc
  USE qpoint,        ONLY : xq, npwq, igkq, nksq
  USE gwsigma,       ONLY : sigma_x_st, nbnd_sig
  USE buffers,       ONLY : save_buffer, get_buffer
  USE io_global,     ONLY : stdout, ionode_id, ionode
  USE gvect,         ONLY : nl, ngm, g, nlm, gstart, gl, igtongl

IMPLICIT NONE

!ARRAYS to describe exchange operator.
  LOGICAL :: limit, lgamma
  REAL(DP) :: rcut, spal
  INTEGER :: ikmq, ik0, ik
  INTEGER :: ig, igp, ir, irp
  INTEGER :: iq, ipol, ibnd, jbnd, counter, ios
  REAL(DP) :: qg2, qg 
  COMPLEX(DP) :: ZDOTC
  COMPLEX(DP) :: czero, exch_element
!q-vector of coulomb potential xq_coul := k_{0} - xk(ik)
  REAL(DP) :: xq_coul(3)
!Arrays to handle case where nlsco does not contain all G vectors required for |k+G| < ecut
  INTEGER     :: igkq_ig(npwx) 
  INTEGER     :: igkq_tmp(npwx) 
  INTEGER     :: ikq
  COMPLEX(DP) :: sigma_band_ex(nbnd_sig, nbnd_sig)
  COMPLEX(DP), ALLOCATABLE :: sigma_ex(:,:)
  COMPLEX(DP), ALLOCATABLE :: greenf_na(:,:), greenf_nar(:,:)
  COMPLEX(DP), ALLOCATABLE :: barcoul(:,:), barcoulr(:,:)
  COMPLEX(DP), ALLOCATABLE :: sigma_g_ex(:,:)
! Self-Energy grid:
! iGv
  CALL start_clock('sigma_exch')
  ALLOCATE ( sigma_ex    (sigma_x_st%dfftt%nnr, sigma_x_st%dfftt%nnr) )
  if (nksq.gt.1) rewind (unit = iunigk)
!set appropriate weights for points in the brillouin zone. 
!weights of all the k-points are in odd positions in list.  

  wq(:) = 0.0d0
!check if we're looking at k = gamma.
  lgamma = (abs(xk_kpoints(1,1)).lt.1.0e-10_dp .and. &
            abs(xk_kpoints(2,1)) .lt. 1.0e-10_dp .and. &
            abs(xk_kpoints(3,1)) .lt. 1.0e-10_dp )

  do iq = 1, nksq
     if(lgamma) then
        wq(iq) = 0.5d0*wk(iq)
        print*, wq(iq)
     else
        wq(iq) = 0.5d0*wk(2*iq-1) 
     endif
  enddo

  if (lgamma) then
     ikq = ik0
  else
     ikq = 2*ik0
  endif
  write(6,'(4x,"Sigma exchange for k",i3, 3f12.7)') ik0, (xk_kpoints(ipol, ik0), ipol=1,3)
  czero = (0.0d0, 0.0d0)
  sigma_ex(:,:) = (0.0d0, 0.0d0)
  DO iq = 1, nksq
     write(6,'(4x,"q ",i3, 3f12.7)') iq, (xk(ipol, iq), ipol=1,3)
     if (lgamma) then
        ikq = iq
     else
        ikq = 2*iq
     endif
     call get_buffer (evq, lrwfc, iuwfc, ikq)
     if (nksq.gt.1) then
          read (iunigk, err = 100, iostat = ios) npw, igk
100       CALL errore ('sigma_exch', 'reading igk', abs (ios) )
     endif
     if (lgamma)  npwq = npw
     if (.not.lgamma.and.nksq.gt.1) then
           read (iunigk, err = 200, iostat = ios) npwq, igkq
200        call errore ('sigma_exch', 'reading igkq', abs (ios) )
     endif
!Need a loop to find all plane waves below ecutsco when igkq 
!takes us outside of this sphere.  
     counter  = 0
     igkq_tmp = 0
     igkq_ig  = 0
!Should have a data_type which keeps track of these indices...
     do ig = 1, npwq
        if((igkq(ig).le.sigma_x_st%ngmt).and.((igkq(ig)).gt.0)) then
           counter = counter + 1
           igkq_tmp (counter) = igkq(ig)
           igkq_ig  (counter) = ig
        endif
     enddo
     allocate ( greenf_na   (sigma_x_st%ngmt, sigma_x_st%ngmt) )
!    psi_{k+q}(r)psi^{*}_{k+q}(r')
     greenf_na = (0.0d0, 0.0d0)
     do ig = 1, counter
       do igp = 1, counter
         do ibnd = 1, nbnd
            greenf_na(igkq_tmp(ig),igkq_tmp(igp)) = greenf_na(igkq_tmp(ig), igkq_tmp(igp)) + &
                                            tpi * (0.0d0, 1.0d0) * (evq(igkq_ig(ig),ibnd))* &
                                            conjg((evq(igkq_ig(igp), ibnd)))
         enddo
       enddo
     enddo
!Fourier transform of green's function
     ALLOCATE ( greenf_nar  (sigma_x_st%dfftt%nnr, sigma_x_st%dfftt%nnr)  )
     greenf_nar(:,:) = czero
     call fft6(greenf_na(1,1), greenf_nar(1,1), sigma_x_st, 1)
     DEALLOCATE(greenf_na)
     ALLOCATE ( barcoul  (sigma_x_st%ngmt, sigma_x_st%ngmt) )

     rcut = (float(3)/float(4)/pi*omega*float(nq1*nq2*nq3))**(float(1)/float(3))
     barcoul(:,:) = (0.0d0,0.0d0)
     xq_coul(:) = xk_kpoints(:,ik0) - xk(:,ikq)

     do ig = 1, sigma_x_st%ngmt
        qg = sqrt((g(1,ig)  + xq_coul(1))**2.d0 + (g(2,ig) + xq_coul(2))**2.d0   &
                + (g(3,ig ) + xq_coul(3))**2.d0)

        qg2 =   (g(1,ig)   + xq_coul(1))**2.d0  + (g(2,ig) + xq_coul(2))**2.d0   &
                + ((g(3,ig)) + xq_coul(3))**2.d0
        limit = (qg.lt.eps8)
        if(.not.limit) then
            spal = 1.0d0 - cos (rcut * tpiba * qg)
            barcoul (ig, ig) = e2 * fpi / (tpiba2*qg2) * dcmplx(spal, 0.0d0)
        else
            barcoul(ig,ig)= (fpi*e2*(rcut**2))/2 
        endif
     enddo
     ALLOCATE (barcoulr    (sigma_x_st%dfftt%nnr,  sigma_x_st%dfftt%nnr))
     barcoulr(:,:) = (0.0d0, 0.0d0)
     call fft6(barcoul(1,1), barcoulr(1,1), sigma_x_st, 1)
     DEALLOCATE(barcoul)
     sigma_ex = sigma_ex + wq(iq)* (0.0d0,1.0d0) / tpi *  greenf_nar * barcoulr
     DEALLOCATE(barcoulr)
     DEALLOCATE(greenf_nar)
   ENDDO

   ALLOCATE ( sigma_g_ex  (sigma_x_st%ngmt, sigma_x_st%ngmt) )
   sigma_g_ex(:,:) = (0.0d0,0.0d0)
   call fft6(sigma_g_ex, sigma_ex, sigma_x_st, -1)
   CALL davcio(sigma_g_ex, lrsex, iunsex, ik0, 1)
   DEALLOCATE(sigma_g_ex)
   DEALLOCATE (sigma_ex)
   CALL stop_clock('sigma_exch')
END SUBROUTINE sigma_exch
