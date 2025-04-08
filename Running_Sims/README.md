This folder contains all of the files needed to run 1D range expansion simulations with asymmetric game 
interactions as outlined in:
https://www.biorxiv.org/content/10.1101/2024.12.14.628506v1.

###########################################################################################################
To change the interaction parameter sets one wants to test, one must modify the "runningall.sh" file. 
The current file looks at:

Pwmm=-1.0, Pmwm=1.0 and correspondingly Pwmw=-Pwmm, Pmww=-Pmwm

Pwmm=0.5, Pmwm=-0.25 and correspondingly Pwmw=-Pwmm, Pmww=-Pmwm

Pwmm=0.5, Pmwm=0.75 and correspondingly Pwmw=-Pwmm, Pmww=-Pmwm

Pwmm=-0.25, Pmwm=0.75 and correspondingly Pwmw=-Pwmm, Pmww=-Pmwm

Each successive chunk of code starting at line 14 is dependent on the last. This means that (for example) 
simulations of job_id=1234 and task_id=23 corresponding to Pwmm=0.5, Pmwm=-0.25 won't run until job_id=4321 
and task_id=23 corresponding to Pwmm=-1.0, Pmwm=1.0 are finished (notice that the task_id is the same). 
The task_id is the id of one of the tasks in the task array in slurm.

For these parameter sets one needs to make the directories:

Pwmm_-1.0_Pwmw_1.0_Pmwm_1.0_Pmww_-1.0/

Pwmm_0.5_Pwmw_-0.5_Pmwm_-0.25_Pmww_0.25/

Pwmm_0.5_Pwmw_-0.5_Pmwm_0.75_Pmww_-0.75/

Pwmm_-0.25_Pwmw_0.25_Pmwm_0.75_Pmww_-0.75/

In order for the files to work properly one needs to run the following commands:

chmod +x runningall.sh

chmod +x replace_vals.sh

cp job_surv_probs_default.slurm job_surv_probs.slurm 

cp surv_probs_re_gill_games_default.jl surv_probs_re_gill_games.jl 

for d in */; do cp *surv_probs* continuous_init_wave-rw_0.* $d; done

###########################################################################################################

In order to simulate arbitrary wild-type growth rates (rw) one needs to create a new initial wave profile 
file (e.g. continuous_init_wave-rw_0.5-K_100-L_302.txt). The files corresponding to rw=0.1, 0.5, or 0.9 are 
available already. 

***IMPORTANT: You must first update the job_surv_probs_default.slurm file to account for changes in where 
your julia executable is, and specific cluster requirements for memory and other architectural things if 
not using the CWRU pioneer cluster. 

For a given number of total trajectories per initiation site (for example 1024), total number of simulation 
time steps (for example 10,000,000), an intrinsic wild-type growth rate (for example rw=0.5), and an 
intrinsic mutant growth rate (for example rm=0.9), one would run the following command:

./runningall.sh 1024 10000000 0.5 0.9

###########################################################################################################

To track the progress of your running sims you can run:

./progress 64 1024 150

Where 64 is the size of the array in job_surv_probs_default.slurm, 1024 is the total number of sims, and 
150 is the number of mutant initiation sites.
