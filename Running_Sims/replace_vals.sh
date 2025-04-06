
#!/bin/bash
cp surv_probs_re_gill_games_default.jl surv_probs_re_gill_games.jl;
cp job_surv_probs_default.slurm job_surv_probs.slurm
sed -i "s/Pwmm=-1.0/Pwmm=$1/g" surv_probs_re_gill_games.jl ;
sed -i "s/Pmwm=1.0/Pmwm=$2/g" surv_probs_re_gill_games.jl;
sed -i "s/r_w=0.1/r_w=$5/g" surv_probs_re_gill_games.jl ;
sed -i "s/r_m=0.1/r_m=$6/g" surv_probs_re_gill_games.jl ;
sed -i "s/T=3_000_000/T=$3/g" surv_probs_re_gill_games.jl;
sed -i "s/T_3000000/T_$3/g" job_surv_probs.slurm;
sed -i "s/NUMITERS/$4/g" job_surv_probs.slurm;

n1=$1;
n2=$2;
nn1=$(julia -e "println(-($n1))");
nn2=$(julia -e "println(-($n2))");

echo -e "Pwmm=$n1";
echo -e "Pwmw=$nn1";
echo -e "Pmwm=$n2";
echo -e "Pmww=$nn2\n";

dirname="Pwmm_${n1}_Pwmw_${nn1}_Pmwm_${n2}_Pmww_${nn2}";
mkdir $dirname;
cd $dirname;
cp ../job_surv_probs.slurm .;
cp ../surv_probs_re_gill_games.jl .;
cp ../continuous_init_wave-rw_$5-K_100-L_302.txt .
cp ../get_iter_dist.jl .
cp ../iter_dist.sh .
#sbatch job_surf_probs.slurm $n1 $nn1 $n2 $nn2;
