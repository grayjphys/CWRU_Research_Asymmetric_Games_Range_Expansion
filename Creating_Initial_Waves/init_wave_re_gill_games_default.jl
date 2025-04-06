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

println(
    "---------------------------------------------------------------------------",
)
t_start = now()
println("Started: ", t_start)



time_steps=20000000 # default value of the number of timesteps in a given simulation. *Replaced by replace_vals.sh.
num_snaps = time_steps÷10
snap = time_steps÷num_snaps
M = 200 # box size minus boundaries.
K = 100 # carrying capacity
shift_zone=150 
barrier_zone=191
r_w=0.1  # default value of the number of wild-type growth rate in a given simulation. *Replaced by replace_vals.sh
Pwvw= 1.0 # birth wild regular value 1.0
Pwvv=-1.0 # regular value -1.0
Pvwv=1-r_w # death wild regular value 1.0-r_w
Pvww=-(1-r_w) # regular value -(1.0-r_w)

# initial wave profile is initially 300 sites long, this just takes the front 200.

init_wild_wave = convert(Vector{Int8},ceil.(readdlm("continuous_init_wave-rw_0.1-K_$(Int64(K))-L_302.txt")[1:M+2]))
state_init = [init_wild_wave Int8(K) .- init_wild_wave]

state = deepcopy(state_init) # sets the state to the initial wild-type wave profile
states = zeros((num_iterations,M+2+1)) # the state at the end of every trajectory
dims = size(state) # the dimensions of the state matrix
state_plus = zeros(Int8,dims) # the state shifted to the left by one
state_minus = zeros(Int8,dims) # the state shifted to the right by one
prop_mat = zeros((M+2,6)) # propensity matrix
prob_mat = zeros((M+2,6)) # probability matrix
time = 0.0 # the initial Gillespie time
sim_time = 0.0 # the initial simulation time
R = 0.0 # the sum of probabilities
r = 0.0 # a random number between 0 and 1
loc = 0 # the location of a randomly chosen event
event = 0 # the type of a randomly chosen event

# the location of output data files
dirname_gill = "./Init-Wave-Pwvw_$(Pwvw)-Pwvv_$(Pwvv)-Pvwv_$(Pvwv)-Pvww_$(Pvww)-T_$(time_steps)-K_$(Int64(K))-M_$(M)-Zone_$(shift_zone)"
mkpath(dirname_gill)

# the stoichiometric matrix
#                       bw                                   dw                                    wv+       wv-     vw+       vw-    
change_matrix_left = [ [0,0]                               [0,0]                               [0,0]  [1,-1] [0,0]  [-1,1]]
change_matrix_bulk = [ [Int8(sign(Pwvw)),Int8(sign(Pwvv))] [Int8(sign(Pvww)),Int8(sign(Pvwv))] [-1,1] [-1,1] [1,-1] [1,-1]]
change_matrix_right = [[0,0]                               [0,0]                               [1,-1] [0,0]  [-1,1] [0,0]]

global change_tensor_define = [change_matrix_left[:,1];;change_matrix_bulk[:,1];;change_matrix_right[:,1]]'
for i in 2:6
    global change_tensor_define = [change_tensor_define;;;[change_matrix_left[:,i];;change_matrix_bulk[:,i];;change_matrix_right[:,i]]'][:,:,:]
end

change_tensor = zeros(Int8,size(change_tensor_define))
change_tensor .= change_tensor_define 
change_tensor_define = nothing

# the simulation function
function gill(
    r_w::Float64,
    M::Int64,
    K::Int64,
    time_steps::Int64,
    change_tensor::Array{Int8,3},
    state_init::Matrix{Int8},
    dirname_gill::String,
    state::Matrix{Int8},
    dims::Tuple{Int64, Int64},
    state_plus::Matrix{Int8},
    state_minus::Matrix{Int8},
    prop_mat::Matrix{Float64},
    prob_mat::Matrix{Float64},
    time::Float64,
    sim_time::Float64,
    shift_zone::Int64,
    barrier_zone::Int64,
    R::Float64,
    r::Float64,
    loc::Int64,
    event::Int64
)   
          
    state .= state_init # sets the state to the initial wild-type wave profile

    state_plus .= zeros(Int8,dims)
    state_minus .= zeros(Int8,dims)
    prop_mat .= zeros((M+2,6))
    prob_mat .= zeros((M+2,6))
    time = 0.0
    sim_time = 0.0
    R = 0.0
    r = 0.0
    loc = 0
    event = 0

    @inbounds for t in 1:time_steps # updates the state matrix for each time step
        state_plus .= circshift(@views(state), -1)
        state_plus[end,:] .= [0,0]
        state_minus .= circshift(@views(state), 1)
        state_minus[1,:] .= [0,0]

        # sets the probability matrix as the propensities times the rates
        prop_mat .= [abs(Pwvw).* @views(state[:,1]) .* @views(state[:,2]);; # bw
                        abs(Pvww).*@views(state[:,1]) .* @views(state[:,2]);; # dw
                        0.5 .* @views(state[:,1]) .* @views(state_plus[:,2]);; # wv+
                        0.5 .* @views(state[:,1]) .* @views(state_minus[:,2]);; # wv-
                        0.5 .* @views(state[:,2]) .* @views(state_plus[:,1]);; # vw+
                        0.5 .* @views(state[:,2]) .* @views(state_minus[:,1])]# vw -

        R = sum(prop_mat) # sum of all probabilities of all events.
        tau = -R*log(rand()) # gets the change in time.
        time += tau # updates the time.
        sim_time += 1.0 # updated the sim time.
        prob_mat .= prop_mat ./ R # normalizes the probability matrix.
        event_arguments = findall(>(0.0), prob_mat) # gets all nonzero probability events
        prob_event_vec = prob_mat[event_arguments] # creates vector of event probabilities
        num_events = length(event_arguments) # number of events to choose from
        seq_prob_vec = zeros(num_events+1) # cumulative probabilities

        for i in 2:num_events+1
            seq_prob_vec[i] = seq_prob_vec[i-1]+prob_event_vec[i-1]
        end

        r = rand() # gets a random number between 0 and 1
        loc = 0
        event = 0
        @inbounds for i in 2:num_events+1 # selects event to occur
            if r >= seq_prob_vec[i-1] && r < seq_prob_vec[i]
                loc = event_arguments[i-1][1]
                event = event_arguments[i-1][2]
            end
        end
        
        event_arguments = nothing
        prob_event_vec = nothing
        num_events = nothing
        seq_prob_vec = nothing

        state[loc-1:loc+1,:] .+= change_tensor[:,:,event] # changes the state based on the event that occurs

        if mod(t,snap) == 0 # shifts the state if it goes past the barrier
            arguments=findall(x -> x != 0,state[:,1])
            tip_loc = maximum([arguments[i,1][1] for i in 1:size(arguments)[1]])
            if tip_loc > barrier_zone 
                state .= circshift(@views(state),(M+2-tip_loc+shift_zone-1,0))
                state[shift_zone+1:end,1] .= 0
                state[shift_zone+1:end,2] .= K
            end                    
        end
    end

    return [state[:,1]; time]
end

function parallel_aggregate_states(r_w::Float64,
                                     M::Int64,
                                     K::Int64,
                                     time_steps::Int64,
                                     change_tensor::Array{Int8,3},
                                     state_init::Matrix{Int8},
                                     dirname_gill::String,
                                     state::Matrix{Int8},
                                     dims::Tuple{Int64, Int64},
                                     state_plus::Matrix{Int8},
                                     state_minus::Matrix{Int8},
                                     prop_mat::Matrix{Float64},
                                     prob_mat::Matrix{Float64},
                                     time::Float64,
                                     sim_time::Float64,
                                     shift_zone::Int64,
                                     barrier_zone::Int64,
                                     R::Float64,
                                     r::Float64,
                                     loc::Int64,
                                     event::Int64)
    states = zeros((num_iterations,M+2+1))
    @inbounds for n in 1:num_iterations
        states[n,:] .= gill(r_w,
            M,
            K,
            time_steps,
            change_tensor,
            state_init,
            dirname_gill,
            state,
            dims,
            state_plus,
            state_minus,
            prop_mat,
            prob_mat,
            time,
            sim_time,
            shift_zone,
            barrier_zone,
            R,
            r,
            loc,
            event)
    end
    return states
end

states .= parallel_aggregate_states(r_w,
                                    M,
                                    K,
                                    time_steps,
                                    change_tensor,
                                    state_init,
                                    dirname_gill,
                                    state,
                                    dims,
                                    state_plus,
                                    state_minus,
                                    prop_mat,
                                    prob_mat,
                                    time,
                                    sim_time,
                                    shift_zone,
                                    barrier_zone,
                                    R,
                                    r,
                                    loc,
                                    event)

# writes all of the data files to "./Init-Wave-Pwvw_$(Pwvw)-Pwvv_$(Pwvv)-Pvwv_$(Pvwv)-Pvww_$(Pvww)-T_$(time_steps)-K_$(Int64(K))-M_$(M)-Zone_$(shift_zone)"
writedlm("$(dirname_gill)/states_matrix_r_w$(r_w)-T$(time_steps)-K$(Int64(K))-M$(M)-tot_num_interations$(tot_num_iterations)-num_interations$(num_iterations)_$(job_name)_$(job_task).txt",states)

# Show Total Runtime
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
