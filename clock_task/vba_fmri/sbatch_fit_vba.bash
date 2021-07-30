#!/usr/bin/env sh

#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=6G
#SBATCH -p general

# env
cd $SLURM_SUBMIT_DIR

export matlab_cpus=16
module use /proj/mnhallqlab/sw/modules

module load matlab/2021a
module load r/4.0.3_depend

#datasets="mmclock_meg mmclock_fmri"
#datasets="mmclock_fmri"
datasets="mmclock_meg"
#models="decay" # fixed decay_factorize decay_uniform decay_ps_equate decay_uniform_ps_equate"
#models="fixed_uv_ureset fixed_uv_baked_ureset"
#models="fixed_uv_ureset_fixedparams_fmri"
#models="fixed_uv_ureset fixed_uv_baked_ureset"
models="fixed"

#factorize="separate factorize"
#decay="uniform selective"
#ps="psnarrow psequate"

#settings for group fixed pararams in fmri
#export sceptic_dataset="mmclock_fmri"
#export sceptic_model="decay_factorize_selective_psequate_fixedparams_fmri"

#matlab -nodisplay -r sceptic_fit_group_vba_ffx 2>&1 | tee ofiles/sceptic_fit_ffx_mmclock_fmri_${sceptic_model}.out

#settings for group fixed pararams in meg
#export sceptic_dataset="mmclock_meg"
#export sceptic_model="decay_factorize_selective_psequate_fixedparams_meg"

#matlab -nodisplay -r sceptic_fit_group_vba_ffx 2>&1 | tee ofiles/sceptic_fit_ffx_${sceptic_dataset}_${sceptic_model}.out

for d in $datasets; do
    for m in $models; do
	export sceptic_dataset=$d
	export sceptic_model="${m}"

	if [ ! -f "ofiles/sceptic_fit_mfx_${d}_${sceptic_model}.out" ]; then
	    #matlab -nodisplay -r sceptic_fit_group_vba_ffx 2>&1 | tee ofiles/sceptic_fit_ffx_${d}_${sceptic_model}.out
	    matlab -nodisplay -r sceptic_fit_group_vba_mfx 2>&1 | tee ofiles/sceptic_fit_mfx_${d}_${sceptic_model}.out
	else
	    echo "Skipping run because sceptic_fit_mfx_${d}_${sceptic_model}.out exists."
	fi
    done
done

#for full examination of factorization, prop spread equality, and uniform versus selective decay
# for d in $datasets; do
#     for m in $models; do
# 	for f in $factorize; do
# 	    for u in $decay; do
# 		for p in $ps; do
# 		    export sceptic_dataset=$d
# 		    export sceptic_model="${m}_${f}_${u}_${p}"

# 		    if [ ! -f "ofiles/sceptic_fit_mfx_${d}_${sceptic_model}.out" ]; then
# 			matlab -nodisplay -r sceptic_fit_group_vba_ffx 2>&1 | tee ofiles/sceptic_fit_ffx_${d}_${sceptic_model}.out
# 			matlab -nodisplay -r sceptic_fit_group_vba_mfx 2>&1 | tee ofiles/sceptic_fit_mfx_${d}_${sceptic_model}.out
# 		    else
# 			echo "Skipping run because sceptic_fit_ffx_${d}_${sceptic_model}.out exists."
# 		    fi
# 		done
# 	    done
# 	done
#     done
# done

R CMD BATCH --no-save --no-restore compile_trial_level_dataframes.R
