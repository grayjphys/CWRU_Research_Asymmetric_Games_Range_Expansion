#!/bin/bash

numruns=$1;
timesteps=$2;

./replace_vals.sh 0.1 $timesteps $numruns;
cd WAVE_rw_0.1/;
job_submit=$(sbatch job_init_wave.slurm 0.1);
job_array_id=$(echo $job_submit | cut -d' ' -f4);
cd ..;

./replace_vals.sh 0.5 $timesteps $numruns;
cd WAVE_rw_0.5/;
job_submit=$(sbatch job_init_wave.slurm 0.5);
job_array_id=$(echo $job_submit | cut -d' ' -f4);
cd ..;

./replace_vals.sh 0.9 $timesteps $numruns;
cd WAVE_rw_0.9/;
job_submit=$(sbatch job_init_wave.slurm 0.9);
job_array_id=$(echo $job_submit | cut -d' ' -f4);
cd ..;
