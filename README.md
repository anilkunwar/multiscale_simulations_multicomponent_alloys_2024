# Codes and Data 

DESCRIPTION OF THE CODES AND DATA HERE

# How to Perform Numerical Simulations using FEM

The finite element analysis has been performed using Elmer multiphysics software. The source codes of the Elmer software are available at https://github.com/elmercsc/elmerfem.

HOW TO INSTALL ELMER SOFTWARE IN UBUNTU OS:

In Ubuntu OS, the most easiest way to perform the installation of this software is by using the launchpad

$ sudo apt-add-repository ppa:elmer-csc-ubuntu/elmer-csc-ppa

$ sudo apt-get update

$ sudo apt-get install elmerfem-csc

HOW TO RUN THE SIMULATION:

Once the installation of the Elmer software is complete, the computer simulation can be performed:
The FORTRAN 90 files have to be compiled using the following command:

$ elmerf90 -o filename.so filename.F90

The user developed codes made available at this repository can be assembled in a directory together with an input file (yet to be made available). The input file can be run with the the following command at the terminal:

$ ElmerSolver inputfilename.sif


The MWE for the input file (sif) file will be made available later.


# How to Perform Molecular Dynamics Simulations

The MD simulations have been performed using LAMMPS software. 

There are many ways to install LAMMPS software in Linux computer. 
It is recommended to install LAMMPS within a given conda environment. 

The following command at the terminal can be employed to run the MD simulations:

$ mpirun -np 7 lmp_mpi -in  filename.in -log data/lammps_sim.log 

The filename.in is the LAMMPS input file. 






