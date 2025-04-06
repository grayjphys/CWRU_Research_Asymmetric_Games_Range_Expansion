#!/bin/bash
#SBATCH --time=13-00:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6000MB
#SBATCH --array=1-64

omrw=$(
echo -n "0";
echo "1.0-$1" | bc;)

echo -e "rw\t$1" > "params_$SLURM_ARRAY_JOB_ID.txt";
echo -e "$SLURM_CPUS_ON_NODE $SLURM_NNODES $SLURM_NTASKS $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_MAX NUMITERS" > "SLURM_params.txt";

mkdir Init-Wave-Pwvw_1.0-Pwvv_-1.0-Pvwv_$omrw-Pvww_-$omrw-T_20000000-K_100-M_200-Zone_150;
cp ../progress.sh .;
cp ../iprog.sh .;

/home/jmg367/JULIA/julia-1.8.0/bin/julia -t $SLURM_CPUS_ON_NODE $PWD/init_wave_re_gill_games.jl \
$SLURM_NNODES \
$SLURM_CPUS_ON_NODE \
$SLURM_ARRAY_TASK_ID \
$SLURM_ARRAY_JOB_ID \
$SLURM_ARRAY_TASK_MAX \
NUMITERS;

#/home/jmg367/JULIA/julia-1.8.0/bin/julia -t $SLURM_NTASKS $PWD/surf_probs_re_gill_games.jl \
#$SLURM_NNODES $SLURM_CPUS_ON_NODE $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_MAX;
