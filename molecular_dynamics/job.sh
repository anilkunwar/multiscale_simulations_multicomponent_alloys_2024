#!/bin/bash

#lammps="../../lammps/src/lmp_serial" # Path to LAMMPS executable.
#lammps="lmp_serial" # conda env lmp_serial
lammps="mpirun -np 7 lmp_mpi"  # conda env lmp_parallel


mkdir -p data # Create directory structure for data output.

${lammps} -in  lmp_gsfe_alloy.in \
          -log data/lammps_sim.log 

