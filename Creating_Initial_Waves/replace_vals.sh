
#!/bin/bash
cp init_wave_re_gill_games_default.jl init_wave_re_gill_games.jl;
cp job_init_wave_default.slurm job_init_wave.slurm
sed -i "s/r_w=0.1/r_w=$1/g" init_wave_re_gill_games.jl ;
sed -i "s/time_steps=20000000/time_steps=$2/g" init_wave_re_gill_games.jl;
sed -i "s/T_20000000/T_$2/g" job_init_wave.slurm;
sed -i "s/NUMITERS/$3/g" job_init_wave.slurm;

n1=$1;

echo -e "rw=$n1";

dirname="WAVE_rw_${n1}";
mkdir $dirname;
cd $dirname;
cp ../job_init_wave.slurm .;
cp ../init_wave_re_gill_games.jl .;
cp ../continuous_init_wave-rw_0.1-K_100-L_302.txt .
cp ../get_iter_dist.jl .
cp ../iter_dist.sh .
