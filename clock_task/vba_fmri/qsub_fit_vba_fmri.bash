#!/usr/bin/env sh

#PBS -l nodes=1:ppn=37
#PBS -A mnh5174_collab
#PBS -j oe
#PBS -M michael.hallquist@psu.edu
#PBS -m abe

env
cd $PBS_O_WORKDIR

module load matlab
export matlab_cpus=37

matlab -nodisplay -r sceptic_fit_group_vba_fmri

