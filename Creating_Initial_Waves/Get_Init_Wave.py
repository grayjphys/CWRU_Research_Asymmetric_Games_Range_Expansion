import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# rw=float(sys.argv[1])
rw = 0.1
tiploc=75

states_dir=os.fsencode("/home/jmg367/GAMES-SURFING/WAVE_rw_%s/Init-Wave-Pwvw_1.0-Pwvv_-1.0-Pvwv_%s-Pvww_%s-T_%s-K_100-M_200-Zone_150"%(rw,1-rw,-1+rw,1000000))

states_data = []
for file in os.listdir(states_dir):
    if b"states_matrix" in file:
        states_data.append(np.loadtxt(os.path.join(states_dir, file).decode("utf-8")))

states_data=np.array(states_data)
states_data_array = []
times = []
for state_data in states_data:
    if len(state_data.shape) == 1:
        states_data_array.append(state_data[:-1])
        times.append(state_data[-1])
    else:
        for state in state_data:
            states_data_array.append(state[:-1])
            times.append(state[-1])

states_data_array=np.array(states_data_array)
loc_50=np.zeros(len(states_data_array),dtype=int)
for i,state in enumerate(states_data_array):
    loc_50[i] = np.argmin(np.abs(state - 50))

min_loc_50 = np.min(loc_50)

for i,state in enumerate(states_data_array):
    states_data_array[i] = np.roll(state,min_loc_50-loc_50[i])
    if min_loc_50-loc_50[i] != 0:
        states_data_array[i][min_loc_50-loc_50[i]:]=0.0

ave_wave = np.mean(states_data_array,axis=0)
minarg = np.argmin(ave_wave[np.nonzero(ave_wave)])
rolled_wave=np.roll(ave_wave,tiploc-minarg)
rolled_wave[tiploc-minarg:]=0.0

for i in range(302-len(rolled_wave)):
    rolled_wave=np.append(rolled_wave,0.0)

np.savetxt("/home/jmg367/GAMES-SURFING/WAVE-K_100-M_200/WAVE_rw_%s/continuous_init_wave-rw_%s-K_100-L_302.txt"%(rw,rw), rolled_wave, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)

plt.plot(rolled_wave)
plt.savefig("/home/jmg367/GAMES-SURFING/WAVE-K_100-M_200/WAVE_rw_%s/continuous_init_wave-rw_%s-K_100-L_302.png"%(rw,rw), bbox_inches='tight')
