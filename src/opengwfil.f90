SUBROUTINE opengwfil()
  USE units_gw,        ONLY : iuncoul, iungreen, lrgrn, lrcoul, iunsigma, &
                              lrsigma, lrsex, iunsex, iunresid, lrresid, & 
                              iunalphabeta, lralphabeta, iunsigext, lrsigext
  USE io_files,        ONLY : tmp_dir, diropn, seqopn
  USE freq_gw,         ONLY : nfs, nwsigma
  USE gwsigma,         ONLY : sigma_x_st, sigma_c_st
  USE control_gw,      ONLY : multishift, do_serial
  USE wvfct,           ONLY : nbnd,npwx

IMPLICIT  NONE
  LOGICAL :: exst

  
    print*, tmp_dir

    iuncoul = 28
    lrcoul = 2 * sigma_c_st%ngmt * sigma_c_st%ngmt * nfs
    CALL diropn (iuncoul, 'coul', lrcoul, exst)

    iungreen = 31
    lrgrn  = 2 * sigma_c_st%ngmt * sigma_c_st%ngmt
    CALL diropn (iungreen, 'green', lrgrn, exst)

!\Sigma^{c}(\G,\G';\omega)
    iunsigma = 32
    if(do_serial) then
       lrsigma = 2 * sigma_c_st%ngmt * sigma_c_st%ngmt
    else  
       lrsigma = 2 * sigma_c_st%ngmt * sigma_c_st%ngmt * nwsigma
    endif
    CALL diropn(iunsigma, 'sigma', lrsigma, exst)

!\Sigma^{x}(\G,\G';\omega)
    iunsex = 33
    lrsex = 2 * sigma_x_st%ngmt * sigma_x_st%ngmt
    CALL diropn(iunsex, 'sigma_ex', lrsex, exst)

    if(multishift) then
       iunresid = 34
       lrresid  = 2*npwx
       iunalphabeta = 35
       lralphabeta  = 4
    endif

END SUBROUTINE opengwfil
