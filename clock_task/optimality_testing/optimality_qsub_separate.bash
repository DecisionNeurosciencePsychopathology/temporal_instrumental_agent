#!/usr/bin/env sh

#master script to submit one qsub job per agent
nagents=12
matlab_cpus=25 #number of CPUs per qsub job

for (( c=1; c<=$nagents; c++ )); do
#for c in 5 6 7 8; do #models that didn't finish
    #echo "qsub optimality_qsub_worker.bash -l walltime=72:00:00 -l nodes=1:ppn=$matlab_cpus -v matlab_cpus=$matlab_cpus,which_agent=$c,runspercond=10,whichopt=sinusoid,fixbeta=1.0"
    echo "qsub optimality_qsub_worker.bash -l walltime=72:00:00 -l nodes=1:ppn=$matlab_cpus -v matlab_cpus=$matlab_cpus,which_agent=$c,whichopt=sinusoid"
    qsub optimality_qsub_worker.bash -l walltime=72:00:00 -l nodes=1:ppn=$matlab_cpus -v matlab_cpus=$matlab_cpus,which_agent=$c,whichopt=sinusoid

    #echo "qsub optimality_qsub_worker.bash -l walltime=60:00:00 -l nodes=1:ppn=$matlab_cpus -v matlab_cpus=$matlab_cpus,which_agent=$c,whichopt=allequate,runspercond=5"
    #qsub optimality_qsub_worker.bash -l walltime=60:00:00 -l nodes=1:ppn=$matlab_cpus -v matlab_cpus=$matlab_cpus,which_agent=$c,whichopt=allequate,runspercond=5

    #matlab parpool starts seem to be failing sometimes even when using separate profiles. Maybe collision due to multitude of launches? Add a bit of sleep
    sleep 40
done
