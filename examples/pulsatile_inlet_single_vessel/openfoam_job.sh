#!/bin/bash
#PBS -l ncpus=8
#PBS -l mem=20GB
#PBS -l jobfs=100MB
#PBS -l walltime=10:00:00
#PBS -l software=openFOAM
#PBS -l wd

# Unload modules.
module rm openmpi intel-cc intel-fc
module load openmpi/4.0.2

# Load modules and set the values.
. /scratch/m45/cm5094/OpenFOAM/OpenFOAM-2.1.1/etc/bashrc

# decompose the case into the desired number of cores
decomposePar
mpirun ericNonNewtonianImplicitFoam -parallel > $PBS_JOBID.log
# reconstruct the case
reconstructPar
# now remove the processor folders to save space
rm -rf processor*
# calculate the wall shear stress
ericWallTractionShearStress
# calculate get the vtk data
foamToVTK