!     This program obtains Vxc(r) and calculates Vxc(r)*|psi_n(r)|^2,
!     to be integrated over to get <Vxc>

      PROGRAM vxc_sandwich

       IMPLICIT NONE

       CHARACTER(LEN=35)             :: statelist, inputfile,outputfile
       INTEGER                       :: nstates, istate, iat, nat
       INTEGER                       :: selected_k, selected_i
       INTEGER                       :: nr1sx, nr2sx, nr3sx, ix, iy, iz
       DOUBLE PRECISION, ALLOCATABLE :: states(:,:,:,:), v_xc(:,:,:)
       DOUBLE PRECISION              :: temp(6), sandwich, checkcharge
       DOUBLE PRECISION              :: a1(3), a2(3), a3(3), vox_vol


!  We first need V_xc
!  Can't obtain this directly from pp; need (V_bare+V_har+V_xc) - (V_bare+V_xc)

 
       OPEN(41,FILE='tot_potential')
       OPEN(42,FILE='vh_vbare_potential')
       OPEN(50,FILE='v_xc.dat')

! READ THE CUBE FILE HEADER
! First two lines are comments
          READ(41,*)
          READ(41,*)
! Third line gives number of atoms and origin
          READ(41,*) nat
! Then, number of voxels and spacing:
          READ(41,*)nr1sx, a1
          READ(41,*)nr2sx, a2
          READ(41,*)nr3sx, a3

          CALL tripleproduct(a1,a2,a3,vox_vol)
! Then positions       
          DO iat = 1, nat
             READ(41,*)
          END DO
          
! Allocate array for v_xc
        
          ALLOCATE(v_xc(nr1sx,nr2sx,nr3sx))

           DO iat = 1, nat + 6

              READ(42,*)

           END DO

! We're now at right place in both files

           DO ix = 1, nr1sx

             DO iy = 1, nr2sx

                DO iz = 1, nr3sx, 6
                   IF((iz+5).LT.(nr3sx+1)) THEN
                   READ(41,*) v_xc(ix,iy,iz:iz+5)
                   READ(42,*) temp
                   v_xc(ix,iy,iz:iz+5) = v_xc(ix,iy,iz:iz+5)  - temp
                   WRITE(50,*)v_xc(ix,iy,iz:iz+5)
                   ELSE
                   READ(41,*) v_xc(ix,iy,iz:nr3sx)
                   READ(42,*) temp(1:nr3sx+1-iz)
                   v_xc(ix,iy,iz:nr3sx) = v_xc(ix,iy,iz:nr3sx)  - temp(1:nr3sx+1-iz)
                   WRITE(50,*)v_xc(ix,iy,iz:nr3sx)
                   END IF
                END DO

              END DO

           END DO
     
      CLOSE(41)
      CLOSE(42)
      CLOSE(50)
! Choose which states you want to calculate the product for
!  input contains number of states, and then their filenames
       statelist = 'input.txt'



       OPEN(40,FILE=statelist)
       READ(40,*) nstates

 !      ALLOCATE(states(nstates,nr1sx,nr2sx,nr3sx))
       DO istate = 1, nstates
          checkcharge = 0.0
          sandwich = 0.0

          READ(40,*)inputfile, selected_k, selected_i
!          WRITE(outputfile,'(I3,A4)')istate,'.dat'
 !         outputfile = ADJUSTL(outputfile)
          OPEN(30,FILE=inputfile)

        !  OPEN(51,FILE=TRIM(outputfile))


          DO iat = 1, nat + 6

             READ(30,*)

          END DO


!  We're now at the right place at the cube file

!  Read all the densities

          DO ix = 1, nr1sx

             DO iy = 1, nr2sx

                DO iz = 1, nr3sx, 6

                   IF((iz+5).LT.(nr3sx+1)) THEN
                   READ(30,*) temp
                   checkcharge = checkcharge + SUM(temp)
                   temp = temp * v_xc(ix,iy,iz:iz+5)
                   sandwich = sandwich + SUM(temp)
                   ELSE
                   READ(30,*) temp(1:nr3sx+1-iz)
                   checkcharge = checkcharge + SUM(temp(1:nr3sx+1-iz))
                   temp(1:nr3sx+1-iz) = temp(1:nr3sx+1-iz) * v_xc(ix,iy,iz:nr3sx)
                   sandwich = sandwich + SUM(temp(1:nr3sx+1-iz))
                   END IF

                   END DO

              END DO

           END DO
 
           CLOSE(30)

       WRITE(*,*)selected_k, selected_i,sandwich*vox_vol, checkcharge*vox_vol
       !    CLOSE(51)

       END DO



       CLOSE(40)




 !      DEALLOCATE(states,v_xc)
 !      DEALLOCATE(v_xc)

      END PROGRAM vxc_sandwich


  SUBROUTINE tripleproduct(vone,vtwo,vthree,volume)
   DOUBLE PRECISION vone(3), vtwo(3), vthree(3)
   DOUBLE PRECISION volume
   DOUBLE PRECISION temp(3)
   temp(1) = vone(2)*vtwo(3) - vtwo(2)*vone(3)
   temp(2) = vone(3)*vtwo(1) - vtwo(3)*vone(1)
   temp(3) = vone(1)*vtwo(2) - vtwo(1)*vone(2)
 
   volume =  DOT_PRODUCT(temp,vthree)
  END SUBROUTINE tripleproduct

