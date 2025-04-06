##########################################################
# Sets the current directory to where the file is run
cd(@__DIR__)
ENV["GKSwstype"] = "nul"
##########################################################
# If running on a cluster (might have to modify based on the specific cluster),
# this section finds the most even distribution of the total number of simulations.
# Based on how many tasks there are (how big your array is on the cluster),
# it calculates how many of the simulations to run on each task.
# Based on the task id, the code then grabs the number it should use for the particular
# instance of the code running (the task).

using Distributed
using ClusterManagers

nn_passed = 0
np_passed = 0
job_task_passed = 0
job_name_passed = 0
arr_size_passed = 0
tot_num_iterations_passed = 0

nn = 0
np = 0
job_task = 0
job_name = 0
arr_size = 0
tot_num_iterations = 0

iterations_vec = zeros(Int64,arr_size)
num_iterations = 1

if length(ARGS) > 0
    nn_passed = parse(Int64,ARGS[1])
    np_passed = parse(Int64,ARGS[2])
    job_task_passed = parse(Int64,ARGS[3])
    job_name_passed = parse(Int64,ARGS[4])
    arr_size_passed = parse(Int64,ARGS[5])
    tot_num_iterations_passed = parse(Int64,ARGS[6])

    nn = nn_passed
    np = np_passed
    job_task = job_task_passed
    job_name = job_name_passed
    arr_size = arr_size_passed
    tot_num_iterations = tot_num_iterations_passed

    iterations_vec = zeros(Int64,arr_size)
    iterations_vec[1]=Int64(ceil(tot_num_iterations/arr_size))

    if arr_size > 1
        for i in 2:arr_size
            N = (tot_num_iterations-sum(iterations_vec))*1.0
            iterations_vec[i] = Int64(ceil(N/(arr_size-i+1)))
        end
    end

    num_iterations = iterations_vec[job_task]
end

#########################################################

using DelimitedFiles
using Dates
using StatsBase

K = Int16(100)
T=3_000_000 # default value of the number of timesteps in a given simulation. *Replaced by replace_vals.sh.
snaps=3_000
M=400 # Box size minus boundaries.
L=(M+2)÷2-1 # Number of locations where mutants are initiated.
EM = 18*M # Number of possible events (swapping, birth, death, etc.) times the number of sites gives the total number of unique events. 

r_m=0.1 # default value of the number of mutant growth rate in a given simulation. *Replaced by replace_vals.sh
Pmvm=1.0 
Pvmm=-(1.0-r_m)
Pwmm=-1.0 # default value of the P_{wm}^m asymmetric payoff tensor component in a given simulation. *Replaced by replace_vals.sh
Pmwm=1.0 # default value of the P_{mw}^m asymmetric payoff tensor component in a given simulation. *Replaced by replace_vals.sh

r_w=0.1 # default value of the number of wild-type growth rate in a given simulation. *Replaced by replace_vals.sh
Pwvw=1.0
Pvww=-(1.0-r_w)
Pwmw=-Pwmm # this is required in order to have a finite carrying capacity.
Pmww=-Pmwm # this is required in order to have a finite carrying capacity.

Pwvv=-Pwvw # this is required in order to have a finite carrying capacity.
Pvwv=-Pvww # this is required in order to have a finite carrying capacity.
Pmvv=-Pmvm # this is required in order to have a finite carrying capacity.
Pvmv=-Pvmm # this is required in order to have a finite carrying capacity.
halfD = 0.5 # factor which sets the relative timescale of interactions.

post = Int64(ceil((M+2)÷2*1.8)) # The location of the barrier, where when the population crosses it, moves the co-moving reference frame.

# initial wave profile is initially 300 sites long, this just adds more sites to the end.
init_wild_wave = convert(Vector{Int16},[ceil.(readdlm("continuous_init_wave-rw_$(r_w)-K_$(K)-L_302.txt").*K./100
                                             )[1:min(M+2,302)]; zeros(max(0,M-300))])
state_init = [init_wild_wave'; zeros(Int16,M+2)'; K .- init_wild_wave' ] # the initial state of the system, a 3x402 matrix.

# creates the directory for the output files.
dirname_gill = "./Games-Surfing-rw_$(r_w)-rm_$(r_m)-Pwmm_$(Pwmm)-Pmwm_$(Pmwm)-T_$(T)-K_$(K)-M_$(M)"
mkpath(dirname_gill)

# the following code builds the stoichiometric matrix (tensor).
change_matrix_left =   [[0,0,0];; # bw                                
                        [0,0,0];; # dw
                        [0,0,0];; # bm
                        [0,0,0];; # dm
                        [0,0,0];; # wm
                        [0,0,0];; # mw
                        [0,0,0];; # wv+           
                        [1,0,-1];;# wv-
                        [0,0,0];; # vw+
                        [-1,0,1];;# vw-
                        [0,0,0];; # mv+
                        [0,1,-1];;# mv-
                        [0,0,0];; # vm+
                        [0,-1,1];;# vm-
                        [0,0,0];; # wm+
                        [1,-1,0];;# wm-
                        [0,0,0];; # mw+
                        [-1,1,0]] # mw-

change_matrix_center = [[Int16(sign(Pwvw)),0,Int16(sign(Pwvv))];;# bw   
                        [Int16(sign(Pvww)),0,Int16(sign(Pvwv))];;# dw
                        [0,Int16(sign(Pmvm)),Int16(sign(Pmvv))];;# bm
                        [0,Int16(sign(Pvmm)),Int16(sign(Pvmv))];;# dm
                        [Int16(sign(Pwmw)),Int16(sign(Pwmm)),0];;# wm
                        [Int16(sign(Pmww)),Int16(sign(Pmwm)),0];;# mw
                        [-1,0,1];;# wv+  
                        [-1,0,1];;# wv-
                        [1,0,-1];;# vw+
                        [1,0,-1];;# vw-
                        [0,-1,1];;# mv+
                        [0,-1,1];;# mv-
                        [0,1,-1];;# vm+
                        [0,1,-1];;# vm-
                        [-1,1,0];;# wm+
                        [-1,1,0];;# wm-
                        [1,-1,0];;# mw+
                        [1,-1,0]] # mw-

change_matrix_right =  [[0,0,0];; # bw
                        [0,0,0];; # dw
                        [0,0,0];; # bm
                        [0,0,0];; # dm
                        [0,0,0];; # wm
                        [0,0,0];; # mw
                        [1,0,-1];;# wv+  
                        [0,0,0];; # wv-
                        [-1,0,1];;# vw+
                        [0,0,0];; # vw-
                        [0,1,-1];;# mv+
                        [0,0,0];; # mv-
                        [0,-1,1];;# vm+
                        [0,0,0];; # vm-
                        [1,-1,0];;# wm+
                        [0,0,0];; # wm-
                        [-1,1,0];;# mw+
                        [0,0,0]]  # mw-

change = convert(Array{Int16,3},permutedims([change_matrix_left;;; 
                                             change_matrix_center;;; 
                                             change_matrix_right],(1,3,2))) # the stoichiometric matrix

function sim(M::Int64,
             L::Int64,
             state_init::Matrix{Int16},
             change::Array{Int16,3},
             EM::Int64,
             T::Int64,
             snaps::Int64,
             K::Int16,
             r_m::Float64,
             Pmvm::Float64,
             Pvmm::Float64,
             Pwmm::Float64,
             Pmwm::Float64,
             r_w::Float64,
             Pwvw::Float64,
             Pvww::Float64,
             Pwmw::Float64,
             Pmww::Float64,
             Pwvv::Float64,
             Pvwv::Float64,
             Pmvv::Float64,
             Pvmv::Float64,
             halfD::Float64,
             post::Int64
             )
    # The absolute values of the payoff tensor components. Required for the probability matrix. 

    Amvm=abs(Pmvm)
    Avmm=abs(Pvmm)
    Awmm=abs(Pwmm)
    Amwm=abs(Pmwm)
    Awvw=abs(Pwvw)
    Avww=abs(Pvww)
    Awmw=abs(Pwmw)
    Amww=abs(Pmww)
    Awvv=abs(Pwvv)
    Avwv=abs(Pvwv)
    Amvv=abs(Pmvv)
    Avmv=abs(Pvmv)

    prob_vec = zeros(EM) # A matrix (vector of vectors) of probabilities of all possible events, initiated to zeros.
    args = collect(1:EM) # Get arguments of these events. Later will be used to get events with non-zero probability. 
    R=0.0 # The sum of all propensities, initiated zero.
    r=0.0 # a random number between 0 and 1, initiated to zero. 
    arg=0 # The argument of the actual event which happens, initiated to zero.
    event=0 # The number of the type of event that occurs (1-18), initiated to zero. 
    loc=0 # The location that the event takes place, initiated to zero. 
    delta=zeros(Int16,(2,3)) # The change that happens to the system (each event modifies a maximum of two sites and three species), initiated to zeros.
    tau = 0.0 # The increase in Gillespie time that occurs during the event, initiated to zero.
    time = 0.0 # The Gillespie time.
    sim_time = 0.0 # The time given by the for loop over time steps.
    shift_dist = Int64(floor((M+2)÷4)) # The distance to shift the co-moving reference frame. 
    tip_loc = 0 # The location of the tip of the population wave, initiated to zero. 
    new_tip_loc = post - shift_dist # The location that the population wave moves after passing the barrier.

    state = zeros(Int16,(3,M+2)) # The state of the system to be updated.
    data_matrix = zeros((5,L)) # The data to be collected, described later.

    @inbounds for l in 1:L # a new simulation for every initiation site.

        state .= state_init # sets the state to the initial wild-type wave profile

        state[2,l] = 1 # places a mutant at the initiation site
        if state[1,l] > 0 # makes sure that the carrying capacity is maintained with the proper number of vacancies
            state[1,l] -= 1
        else
            state[3,l] -= 1
        end

        prob_vec = zeros(EM)
        args = collect(1:EM)
        R=0.0
        r=0.0
        arg=0
        event=0
        loc=0
        delta=zeros(Int16,(2,3))
        tau = 0.0
        time = 0.0
        sim_time = 0.0
        shift_dist = Int64(floor((M+2)÷4))
        tip_loc = 0
        new_tip_loc = post - shift_dist

        @inbounds for t in 1:T # updates the state matrix for each time step

            # sets the probability matrix as the propensities times the rates
            prob_vec .= [Awvw .* @views(state[1,2:end-1]) .* @views(state[3,2:end-1]); # bw
                        Avww .* @views(state[1,2:end-1]) .* @views(state[3,2:end-1]); # dw
                        Amvm .* @views(state[2,2:end-1]) .* @views(state[3,2:end-1]); # bm
                        Avmv .* @views(state[2,2:end-1]) .* @views(state[3,2:end-1]); # dm
                        Awmw .* @views(state[1,2:end-1]) .* @views(state[2,2:end-1]);  # wm
                        Amwm .* @views(state[1,2:end-1]) .* @views(state[2,2:end-1]); # mw
                        halfD .* @views(state[1,2:end-1]) .* @views(state[3,3:end]); # wv+
                        halfD .* @views(state[1,2:end-1]) .* @views(state[3,1:end-2]); # wv-
                        halfD .* @views(state[3,2:end-1]) .* @views(state[1,3:end]); # vw+
                        halfD .* @views(state[3,2:end-1]) .* @views(state[1,1:end-2]); # vw-
                        halfD .* @views(state[2,2:end-1]) .* @views(state[3,3:end]); # mv+
                        halfD .* @views(state[2,2:end-1]) .* @views(state[3,1:end-2]); # mv-
                        halfD .* @views(state[3,2:end-1]) .* @views(state[2,3:end]); # vm+
                        halfD .* @views(state[3,2:end-1]) .* @views(state[2,1:end-2]); # vm-
                        halfD .* @views(state[1,2:end-1]) .* @views(state[2,3:end]); # wm+
                        halfD .* @views(state[1,2:end-1]) .* @views(state[2,1:end-2]); # wm-
                        halfD .* @views(state[2,2:end-1]) .* @views(state[1,3:end]); # mw+
                        halfD .* @views(state[2,2:end-1]) .* @views(state[1,1:end-2]) # mw-
            ]

            R = sum(prob_vec) # sum of all probabilities of all events.
            prob_vec ./= R # normalizes the probability matrix.
            try 
                arg = sample(args, Weights(prob_vec)) # gets the argument of the event that occurs from prob_vec.    
                event=((arg-1)÷M)+1 # gets the type of event that occurs.
                loc=(arg-1)%M+2 # gets the location of the event that occured.
                delta = @views(change[:,:,event]) # gets the change to the state matrix.
                state[:,loc-1:loc+1] .+= delta # changes the state matrix based on the event that happened.
                tau = -R*log(rand()) # gets the change in time.
                time += tau # updated the time.
                if t%(T÷snaps)==0 # after some computational time steps checks whether the barrier has been reached or surpassed.
                    if sum(@views(state[2,:])) < 1 || sum(@views(state[1,:])) < 1 || (sum(@views(state[1,1:L])) < 1 && sum(@views(state[2,L+2:end])) < 1)
                        break
                    end
                    if sum(@views(state[1:2,post+1])) > 0 
                        state[:,L+1] .= 0 # if the barrier is passed, the stationary and comoving reference frame are decoupled.
                        for j in M+2:-1:1 # gets the location of the population wave-front
                            if state[3,j] < K
                                tip_loc = j
                                break
                            end
                        end 
                        # updates the comoving reference frame by shifting sites cyclically and removing cells.
                        state[:,L+2:end] .= circshift(@views(state[:,L+2:end]),(0,-shift_dist-tip_loc+post))
                        state[:,new_tip_loc+1:end] .= @views(state_init[:,new_tip_loc+1:end])
                    end
                end
            catch
                break
            end
        end

        if sum(@views(state[2,:])) > 0 && sum(@views(state[1,L+2:end])) > 0
            data_matrix[1,l] += 1.0 # survival without surfing
        end

        if sum(@views(state[2,:])) > 0 && sum(@views(state[1,L+2:end])) < 1
            data_matrix[2,l] += 1.0 # surfing counter
            data_matrix[3,l] += time # surfing time
            data_matrix[4,l] += sim_time # surfing sim_time
        end
	
	if sum(@views(state[2,:])) < 1
            data_matrix[5,l] += 1.0 # no survival
        end
    end

    return data_matrix
end


function mc_sim(M::Int64,
                L::Int64,
                state_init::Matrix{Int16},
                change::Array{Int16,3},
                EM::Int64,
                T::Int64,
                snaps::Int64,
                K::Int16,
                r_m::Float64,
                Pmvm::Float64,
                Pvmm::Float64,
                Pwmm::Float64,
                Pmwm::Float64,
                r_w::Float64,
                Pwvw::Float64,
                Pvww::Float64,
                Pwmw::Float64,
                Pmww::Float64,
                Pwvv::Float64,
                Pvwv::Float64,
                Pmvv::Float64,
                Pvmv::Float64,
                halfD::Float64,
                post::Int64,
                data_aggregate::Matrix{Float64}
                )

    @inbounds for n in 1:num_iterations # num iterations calculated based on id and total number of sims. simulates a set of 1-L sims per iteration.
            data_aggregate .+= sim(M,
                                       L,
                                       state_init,
                                       change,
                                       EM,
                                       T,
                                       snaps,
                                       K,
                                       r_m,
                                       Pmvm,
                                       Pvmm,
                                       Pwmm,
                                       Pmwm,
                                       r_w,
                                       Pwvw,
                                       Pvww,
                                       Pwmw,
                                       Pmww,
                                       Pwvv,
                                       Pvwv,
                                       Pmvv,
                                       Pvmv,
                                       halfD,
                                       post
                                       )
    end

    return data_aggregate # all data for a task of the slurm array on the cluster.
end

# Timing the simulations: Nice to have for if things go wrong, or to figure out how to partition your sims on the cluster.
println(
    "---------------------------------------------------------------------------",
)
t_start = now()
println("Started: ", t_start)

data_aggregate = zeros((5,L)) # Initiate the data matrix to zero.

data_aggregate .= mc_sim(M,
                             L,
                             state_init,
                             change,
                             EM,
                             T,
                             snaps,
                             K,
                             r_m,
                             Pmvm,
                             Pvmm,
                             Pwmm,
                             Pmwm,
                             r_w,
                             Pwvw,
                             Pvww,
                             Pwmw,
                             Pmww,
                             Pwvv,
                             Pvwv,
                             Pmvv,
                             Pvmv,
                             halfD,
                             post,
                             data_aggregate
                            )
                            
# Output data file to "./Games-Surfing-rw_$(r_w)-rm_$(r_m)-Pwmm_$(Pwmm)-Pmwm_$(Pmwm)-T_$(T)-K_$(K)-M_$(M)"
writedlm("$(dirname_gill)/fixation_matrix_r_w$(r_w)-r_m$(r_m)-Pwmm$(Pwmm)-Pmwm$(Pmwm)-T$(T)-K$(K)-M$(M)-tot_num_interations$(tot_num_iterations)-num_interations$(num_iterations)_$(job_name)_$(job_task).txt",data_aggregate)

# A progress file, to be read by progress.sh to estimate how far along your sims are. 
writedlm("$(dirname_gill)/prog_num_$(job_name)_$(job_task).txt",[num_iterations])

# Show Total Runtime in terms of hours, minutes, etc.
t_end = now()
println("Ended: ", t_end)

t_diff = (t_end - t_start).value
println(
    t_diff ÷ 3600000,
    " hours ",
    Int(floor((t_diff / 3600000 - t_diff ÷ 3600000) * 60)),
    " minutes ",
    Int(
        floor(
            (
                (t_diff / 3600000 - t_diff ÷ 3600000) * 60 -
                Int(floor((t_diff / 3600000 - t_diff ÷ 3600000) * 60))
            ) * 60,
        ),
    ),
    " seconds",
)
println("---------------------------------------------------------------------------")



