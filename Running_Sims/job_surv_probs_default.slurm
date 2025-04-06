#!/bin/bash
#SBATCH --time=13-00:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6000MB
#SBATCH --array=1-64

omrw=$(
echo -n "0";
echo "1.0-$5" | bc;)

omrm=$(
echo -n "0";
echo "1.0-$6" | bc;)

echo -e "Pwmm\tPwmw\tPmwm\tPmww\n$1\t$2\t$3\t$4" > "params_$SLURM_ARRAY_JOB_ID.txt";
echo -e "$SLURM_CPUS_ON_NODE $SLURM_NNODES $SLURM_NTASKS $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_MAX NUMITERS" > "SLURM_params.txt";

mkdir Games-Surfing-rw_$5-rm_$6-Pwmm_$1-Pmwm_$3-T_3000000-K_100-M_400;

cp ../progress.sh .;
cp ../iprog.sh .;

/home/jmg367/JULIA/julia-1.8.0/bin/julia -t $SLURM_CPUS_ON_NODE $PWD/surv_probs_re_gill_games.jl \
$SLURM_NNODES \
$SLURM_CPUS_ON_NODE \
$SLURM_ARRAY_TASK_ID \
$SLURM_ARRAY_JOB_ID \
$SLURM_ARRAY_TASK_MAX \
NUMITERS;
