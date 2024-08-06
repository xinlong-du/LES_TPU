#!/bin/bash
#----------------------------------------------------
# Sample Slurm job script
#   for TACC Stampede3 SKX nodes
#
#   *** MPI Job in SKX Queue ***
# 
# Last revised: 23 April 2024
#
# Notes:
#
#   -- Launch this script by executing
#      "sbatch skx.mpi.slurm" on Stampede3 login node.
#
#   -- Use ibrun to launch MPI codes on TACC systems.
#      Do not use mpirun or mpiexec.
#
#   -- Max recommended MPI ranks per SKX node: 48
#      (start small, increase gradually).
#
#   -- If you're running out of memory, try running
#      fewer tasks per node to give each task more memory.
#
#----------------------------------------------------

#SBATCH -J myjob           # Job name
#SBATCH -o myjob.o%j       # Name of stdout output file
#SBATCH -e myjob.e%j       # Name of stderr error file
#SBATCH -p skx-dev         # Queue (partition) name
#SBATCH -N 4               # Total # of nodes 
#SBATCH -n 192             # Total # of mpi tasks
#SBATCH -t 2:00:00         # Run time (hh:mm:ss)
#SBATCH --mail-user=du.xinl@northeastern.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A DesignSafe-SimCenter       # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...

ml intel/24.0
source /work2/05660/abiy/stampede3/OpenFOAM/OpenFOAM-10/etc/bashrc

# Launch MPI code... 

ibrun ./myprogram         # Use ibrun instead of mpirun or mpiexec