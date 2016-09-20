#!/usr/bin/env sh
#PBS -l walltime=40:00:00
#PBS -l nodes=1:ppn=40
#PBS -A mnh5174_collab
#PBS -j oe
#PBS -M michael.hallquist@psu.edu
#PBS -m abe

env
cd $PBS_O_WORKDIR

matlab_cpus=40

module load matlab/R2014b

export matlab_cpus

matlab -nodisplay -r TC_identifiability_vetting

