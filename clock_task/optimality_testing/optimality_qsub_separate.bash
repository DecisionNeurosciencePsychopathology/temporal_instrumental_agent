#!/usr/bin/env sh

#master script to submit one qsub job per agent
nagents=13
matlab_cpus=20 #number of CPUs per qsub job

for (( c=1; c<=$nagents; c++ )); do
    qsub optimality_qsub_worker.bash -l walltime=72:00:00 -l nodes=1:ppn=$matlab_cpus -v matlab_cpus=$matlab_cpus,which_agent=$c
    #matlab parpool starts seem to be failing sometimes even when using separate profiles. Maybe collision due to multitude of launches? Add a bit of sleep
    sleep 20
done
