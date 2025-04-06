#!/bin/bash

numruns=$1;
timesteps=$2;
rwc=$3
rmc=$4

./replace_vals.sh -1.0 1.0 $timesteps $numruns $rwc $rmc;
cd Pwmm_-1.0_Pwmw_1.0_Pmwm_1.0_Pmww_-1.0/;
job_submit_0=$(sbatch job_surf_probs.slurm -1.0 1.0 1.0 -1.0 $rwc $rmc);
job_array_id0=$(echo $job_submit_0 | cut -d' ' -f4);
cd ..;

./replace_vals.sh 0.5 -0.25 $timesteps $numruns $rwc $rmc;
cd Pwmm_0.5_Pwmw_-0.5_Pmwm_-0.25_Pmww_0.25/;
job_submit_1=$(sbatch --dependency=aftercorr:$job_array_id0 job_surf_probs.slurm 0.5 -0.5 -0.25 0.25 $rwc $rmc);
job_array_id1=$(echo $job_submit_1 | cut -d' ' -f4);
cd ..;

./replace_vals.sh 0.5 0.75 $timesteps $numruns $rwc $rmc;
cd Pwmm_0.5_Pwmw_-0.5_Pmwm_0.75_Pmww_-0.75/;
job_submit_2=$(sbatch --dependency=aftercorr:$job_array_id1 job_surf_probs.slurm 0.5 -0.5 0.75 -0.75 $rwc $rmc);
job_array_id2=$(echo $job_submit_2 | cut -d' ' -f4);
cd ..;

./replace_vals.sh -0.25 0.75 $timesteps $numruns $rwc $rmc;
cd Pwmm_-0.25_Pwmw_0.25_Pmwm_0.75_Pmww_-0.75/;
job_submit_3=$(sbatch --dependency=aftercorr:$job_array_id2 job_surf_probs.slurm -0.25 0.25 0.75 -0.75 $rwc $rmc);
job_array_id3=$(echo $job_submit_3 | cut -d' ' -f4);
cd ..;
