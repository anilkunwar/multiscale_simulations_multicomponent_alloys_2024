# Codes and Data 

DESCRIPTION OF THE CODES AND DATA HERE

# How to Perform Numerical Simulations using FEM

The finite element analysis has been performed using Elmer multiphysics software.

The FORTRAN 90 files have to be compiled using the following command:

$ elmerf90 -o filename.so filename.F90

The filename.F90 is the user developed fortan file.

The user developed codes made available at this repository can be assembled in a directory together with an input file (yet to be made available). The input file can be run with the the following command at the terminal:

$ ElmerSolver inputfilename.sif


The MWE for the input file (sif) file will be made available later.


# How to Perform Molecular Dynamics Simulations

The MD simulations have been performed using LAMMPS software. 

The following command at the terminal can be employed to run the MD simulations:

$ mpirun -np 7 lmp_mpi -in  filename.in -log data/lammps_sim.log 

The filename.in is the LAMMPS input file. 






