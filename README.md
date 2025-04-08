This repository contains all of the files needed to run all simulations and get figures corresponding to
the paper "Asymmetric Interactions Shape Survival During Population Range Expansions":
https://www.biorxiv.org/content/10.1101/2024.12.14.628506v1.

########################################################################################################### 

The "Creating_Initial_Waves" folder has all of the code necessary for running the simulations of a 1D
range expansion of wild-type cells, and allows one to get an average initial wave profile for input into
range expansion simulations involving a mutant. 

The "Running_Sims" folder has all of the code necessary for running simulations of 1D range expansions
with an initial wild-type wave profile (some wave profiles are already given), and mutants initiated at 
sites within the simulation box, where there are asymmetric game interactions between wild-types and 
mutants. The output of the code is the survival probabilities, surfing probabilities, abiding probabilities, 
and average surfing times of a mutant at each initiation site. 

The "Getting_Figures" directory contains the code which analyzes and plots the output of the previous 
simulations from the Running_Sims folder. One may create figures which shows the surfing/abiing
probabilities, the total survival probabilities along with a numerical solution to the ODE in the paper, 
and average surfing times. 

