#!/usr/bin/env sh

#PBS -l nodes=1:ppn=37
#PBS -l walltime=14:00:00
#PBS -A wff3_a_g_hc_default
#PBS -j oe
#PBS -M michael.hallquist@psu.edu
#PBS -m abe

# mnh5174_a_g_hc_default

# env
cd $PBS_O_WORKDIR

export matlab_cpus=37
module load matlab/R2017b

#datasets="mmclock_meg mmclock_fmri"
datasets="mmclock_fmri"
models="decay fixed decay_factorize decay_uniform decay_ps_equate decay_uniform_ps_equate" 

for d in $datasets; do
    for m in $models; do
	export sceptic_dataset=$d
	export sceptic_model=$m

	if [ ! -f "ofiles/sceptic_fit_mfx_${d}_${m}.out" ]; then
	    matlab -nodisplay -r sceptic_fit_group_vba_ffx 2>&1 | tee ofiles/sceptic_fit_ffx_${d}_${m}.out
	    matlab -nodisplay -r sceptic_fit_group_vba_mfx 2>&1 | tee ofiles/sceptic_fit_mfx_${d}_${m}.out
	else
	    echo "Skipping run because sceptic_fit_ffx_${d}_${m}.out exists."
	fi	
    done
done

R CMD BATCH --no-save --no-restore compile_trial_level_dataframes.R
