SUBROUTINE sgw_opening_message()

    USE kinds, ONLY: DP
    USE io_global,  ONLY :  stdout, ionode, ionode_id
    USE sgw_version, ONLY: sgw_version_number, sgw_svn_revision

    IMPLICIT NONE

    WRITE( stdout, '(/5X, "You are using SGW version: ", A)') sgw_version_number
    WRITE( stdout, '(/5X, "With svn revision number: ", A)') sgw_svn_revision

    WRITE( stdout, '(/5X,"This program is part of the open-source Quantum ", &
         &  "ESPRESSO suite", &
         &/9X," Please also cite H. Lambert and F. Giustino, Phys. Rev. B ",&
         &    " 88. 075117 (2013);", &
         &/9X," URL http://www.sternheimer.org.uk")' )
    RETURN
END SUBROUTINE sgw_opening_message
