This folder contains all of the files needed to run 1D range expansion simulations to get the initial wild-type wave profile as outlined in https://www.biorxiv.org/content/10.1101/2024.12.14.628506v1.

########################################################################################################### To change the initial wild-type growth rate sets one wants to test, one must modify the "runningall.sh" file. The current file looks at:

rw = 0.1, rw = 0.5, and rw = 0.9
###########################################################################################################

***IMPORTANT: You must first update the job_init_wave_default.slurm file to account for changes in where your julia executable is, and specific cluster requirements for memory and other architectural things if not using the CWRU pioneer cluster.

For a given number of total trajectories per initiation site (for example 1024), and a total number of simulation time steps (for example 10,000,000), one would run the following command:

./runningall.sh 1024 10000000

###########################################################################################################
