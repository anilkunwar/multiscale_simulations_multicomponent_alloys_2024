#!/bin/bash
# This script executes sequentially a series of LAMMPS simulations at different temperatures.
# https://icme.hpc.msstate.edu/mediawiki/index.php/LAMMPS_Stacking_Fault_Energy.html

#lammps="../../lammps/src/lmp_serial" # Path to LAMMPS executable.
#lammps="lmp_serial" # conda env lmp_serial
lammps="mpirun -np 7 lmp_mpi"  # conda env lmp_parallel

# Setup list of parameters to loop over.
#T=(   100    400    700   1000   1300   1600)
#a=(2.8841 2.9115 2.9315 2.9484 2.9637 2.9782)
#k=( 5.787  4.866  4.073  3.373  2.799  2.443)

mkdir -p data # Create directory structure for data output.

${lammps} -in  lmp_gsfe_alloy.in \
          -log data/lammps_sim.log 
# Run job.
#for n in $(seq 0 5)
#do
#  printf "Running T = ${T[n]}K simulation.\n"
#  ${lammps} -in  in.lmp                   \
#            -log data/lammps_${T[n]}K.log \
#            -screen none                  \
#            -var RANDOM ${RANDOM}         \
#            -var T ${T[n]}                \
#            -var a ${a[n]}                \
#            -var k ${k[n]}
#done
