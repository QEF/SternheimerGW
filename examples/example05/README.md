
COULOMB PSEUDOPOTENTIAL

Calculation of W(\r,\r';\omega=0), \mu, \mu* for MgB2.
This is a (relatively) quick calculation. The ground
state is calculated on a 10x10x10 grid. W_{\q} is then 
calculated for 84 q-points.

  mpirun -np 8 ~/PW/pw.x -npool 8 < mgb2.scf.in > mgb2.scf.out
  mpirun -np NPROC  ~/SGW/gw.x -nimage N -npool M < mgb2.coul.in > mgb2.coul.out

where N*M = np.

Once W has been calculated mv./tmp/_gw0/mgb2.coul1 to ./tmp/mgb2.coul1
(i.e. up one directory). Then you can run the post-processing code
to perform the Fermi surface matrix element:

  mpirun -np N  ~/SGW/sgwpp/coulmat.x -nimage N < mgb2.sgwpp.in > mgb2.sgwpp.out

!!!!NOTA BENE: the post processing code only accepts image level parallelism (do not pool)!!!
!!!! Also  np must = nimage                                                               !!!

To run the GW code on an MgB2 system:  
mpirun -np NPROC  ~/SGW/gw.x -nimage N -npool M < mgb2.gw.in > mgb2.gw.out

If the calculation fails to finish calculating W. w_of_q_start can be set
to the last completed q-point.


