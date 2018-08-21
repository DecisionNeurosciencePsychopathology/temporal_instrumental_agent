#!/usr/bin/env sh

#PBS -l nodes=1:ppn=37
#PBS -l walltime=08:00:00
#PBS -A mnh5174_a_g_hc_default
#PBS -j oe
#PBS -M michael.hallquist@psu.edu
#PBS -m abe

env
cd $PBS_O_WORKDIR

#module load matlab/R2014b
export matlab_cpus=37
module load matlab/R2017b

# MMClock MEG
export sceptic_dataset=mmclock_meg
export sceptic_model=decay
#matlab -nodisplay -r sceptic_fit_group_vba_fmri
matlab -nodisplay -r sceptic_fit_group_vba_fmri_mfx

export sceptic_dataset=mmclock_meg
export sceptic_model=fixed
#matlab -nodisplay -r sceptic_fit_group_vba_fmri
matlab -nodisplay -r sceptic_fit_group_vba_fmri_mfx

# MMClock fMRI
export sceptic_dataset=mmclock_fmri
export sceptic_model=decay
#matlab -nodisplay -r sceptic_fit_group_vba_fmri
matlab -nodisplay -r sceptic_fit_group_vba_fmri_mfx

export sceptic_dataset=mmclock_fmri
export sceptic_model=fixed
#matlab -nodisplay -r sceptic_fit_group_vba_fmri
matlab -nodisplay -r sceptic_fit_group_vba_fmri_mfx
