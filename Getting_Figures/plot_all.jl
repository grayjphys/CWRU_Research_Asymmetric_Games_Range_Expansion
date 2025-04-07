cd(@__DIR__)
ENV["GKSwstype"] = "nul"

using Plots
using Measures
using LaTeXStrings
using DelimitedFiles
using Interpolations
using StatsBase
using LsqFit

# The possible plots you can make with this code correspond to the surfing/abiding probabilities, survival probabilities, and fixation times
# of a mutant in a 1D population range expansion with asymmetric games. The keywords for each of these are "surfabide", "surv", and "times".
# Admittedly, the average surfing time plots are a bit specific to my data used in the paper, 
# so you might want to write your own code if you are interested in looking at these for your own parameter sets.

to_plot = "surfabide"
save_fig=false # default is to display the figures, make this true to save the figures.

# sets where you want figures saved
curdir=pwd()
fig_dest=curdir*"\\Figs\\"

# If you have the data for the range expansion simulations, the parameters below tell you which data to plot.
rws = [0.1,0.5,0.9] # plots to be made for wild-type growth rates of 0.1, 0.5, and 0.9.
games = [[-1.0 1.0]; [-0.25 0.75]; [0.5 0.75]; [0.5 -0.25]] # the games to be plotted, with [Pwmm Pmwm] for each game.

M=200 # Half of the system size of your sims
add_to_bulk=0 # If you need to shift the wave form to get the numerics to work better, you can add a buffer to the front. I didn't...
K=100 # The carrying capacity
D = 1.0 # The diffusion constant
max_time = 0.0 # Helps with the fixation time figures in scaling the axes, initialized to zero. 

function ReLu(x) # ϕ^+ in the paper
    if x > 0 
        return x
    else return eps(0.0)
    end
end

# γ in the paper, corresponds to the ratio of backwards to forwards transitions at the boundaries, takes in Pwmm and Pmwm, or Pvmm and Pmvm in any order.
function γ(a,b) 
    if a > 0.0
        if b > 0.0
            return 0.0
        elseif b < 0.0
            if a ≈ -b
                return 1.0
            else
                return -b/a
            end
        elseif b ≈ 0.0
            return 0.0
        end
    elseif a < 0.0
        if b > 0.0
            if a ≈ -b
                return 1.0
            else
                return -a/b
            end
        elseif b < 0.0
            return prevfloat(typemax(Float64))
        elseif b ≈ 0.0
            return prevfloat(typemax(Float64))
        end
    elseif a ≈ 0.0
        if b > 0.0
            return 0.0
        elseif b < 0.0
            return prevfloat(typemax(Float64))
        elseif b ≈ 0.0
            return 0.5
        end
    end
end

# Computes the boundary condition for  Pwmm and Pmwm, or Pvmm and Pmvm in any order as a and b. 
# c is a large number, BC terms are approximating an infinite system. 
function BC(a,b,c)
    if a > 0.0
        if b > 0.0
            return 1.0
        elseif b < 0.0
            if a ≈ -b
                return 0.0
            else
                return (1.0 - γ(a,b))/(1-γ(a,b)^c)
            end
        elseif b ≈ 0.0
            return 1.0
        end
    elseif a < 0.0
        if b > 0.0
            if a ≈ -b
                return 0.0
            else
                return (1.0 - γ(a,b))/(1-γ(a,b)^c)
            end
        elseif b < 0.0
            return 0.0
        elseif b ≈ 0.0
            return 0.0
        end
    elseif a ≈ 0.0
        if b > 0.0
            return 1.0
        elseif b < 0.0
            return 0.0
        elseif b ≈ 0.0
            return 0.0
        end
    end
end

figcounter = 0
for gnum in 1:4

    # get game parameters for a specific game
    Pwmm=games[gnum,1] 
    Pmwm=games[gnum,2]

    # gets all of the directories for all of the rws and the specific game, which contain the data for analyzing/plotting.
    # saves it as gamedirs.
    for rw in rws
        figcounter += 1
        gamedirs = []
        for (root, dirs, files) in walkdir(curdir)
            for dir in dirs
                if occursin("Games-Surfing-rw_$(rw)",joinpath(root, dir)) && occursin("Pwmm_$(Pwmm)-Pmwm_$(Pmwm)",joinpath(root, dir)) 
                    if length(readdir(joinpath(root, dir))) > 0 && !occursin("TEST",joinpath(root, dir)) 
                        gamedirs = [gamedirs; joinpath(root, dir)] # path to directories
                    end
                end
            end
        end
        
        # The first dimension of gamedirs corresponds to how many rms per rw per game you did.
        # Assuming you tested the same number of rms for each. 
        num_rms = length(gamedirs) 
        rms = zeros(num_rms) # a list of the rms, to be filled in later.
        endtimes = zeros(num_rms) # not the apocalypse :). gets the times when a mutant surfed, if it surfed.
        colors = [:pink,:violet,:lightblue] # colors for each rm scatter plot or line plot.
        strokecolors = [:red,:purple,:blue] # colors for the outline of each rm scatter plot point.
        markershapes = [:circle,:rect,:utriangle] # shapes for surfing, abiding, and survival scatter plots.

        # The initial wave profile to be read in from the appropriate directory. Takes in the value of rw.
        init_wave = [K .* ones(add_to_bulk); readdlm("./"*readdir(".")[occursin.("continuous_init_wave-rw_$(rw)", readdir("."))][1])[1:min(M,300)]; zeros(max(0,M-300))]./K

        # Below I create an approximate version of the initial wave profile to put into the numerical solver of the survival ode
        # I use a tanh function to approximate it.
        numxs = 30_000 # number of positions 
        xs = LinRange(1.0,M+add_to_bulk,numxs) # the positions
        h = xs[2]-xs[1] # the difference between adjacent position
        w_coeffs0 = [-10.0,1/6] # initial parameters for the fit
        winit_eq(x,coeffs) = ( 0.5 .- 0.5 .*tanh.(coeffs[1].+coeffs[2].*x)) # the fitting function
        fit = curve_fit(winit_eq, LinRange(1.0,M+add_to_bulk,M+add_to_bulk), init_wave, w_coeffs0) # fits the function
        w_coeffs = fit.param # the new parameters of the fit
        ws=winit_eq(xs,w_coeffs) # creates an array of values of the new function for every x in xs. 
        vs = 1.0 .- ws  # creates an array of values for the initial vacancy distribution for every x in xs. 

        p = plot()  # p is for the numerical + scatter survival probability plots
        q = plot()  # q is for the scatter average surfing time plots
        counter = 1 # helps to avoid redundant plotting
        max_time = 0.0 # helps with scaling axes for the average surfing time plots
        min_time = 1e32 # helps with scaling axes for the average surfing time plots

        for (dirnum,dir) in enumerate(gamedirs) # for each rw, rm, and game 

            # Below are calculations of the wave speeds of the initial wave profile 
            # You can optionally calculate wave speeds for the mutant in the bulk and at the front. 

            noise_correction_w_front = (1.0-pi^2/(2*log(K*sqrt(rw))^2)) 
            vw_front=2*sqrt(rw)*noise_correction_w_front
            
            # noise_correction_m_bulk = (1.0-pi^2/(2*log(K*sqrt(Pwmm+Pmwm))^2)) 
            # vm_bulk=2*sqrt(Pwmm+Pmwm)*noise_correction_m_bulk

            # gathers the rms from the file names
            rm=parse(Float64,split(split(dir,"rm_")[2],"-")[1])
            rms[dirnum] = rm

            # noise_correction_m_front = (1.0-pi^2/(2*log(K*sqrt(rm))^2)) 
            # vm_front=2*sqrt(rm)*noise_correction_m_front

            # calculates Pvmm and Pmvm from rm
            Pvmm=-(1.0-rm)
            Pmvm=1.0

            # gets all files from your data directories
            files = [dir*"\\"*file for file in readdir(dir)[occursin.("fixation",readdir(dir))]]
            progs = [dir*"\\"*file for file in readdir(dir)[occursin.("prog",readdir(dir))]] 
            num_trajs = sum([readdlm(prog) for prog in progs])[:][1] # use this to check and see if all of your sims converged

            data = readdlm(files[1])
            for i in 2:length(files)
                data .+= readdlm(files[i])
            end

            num_surfs = data[2,:] # number of trajectories where a mutant surfed for a given rm, rw, and game
            data[1,:] ./= num_trajs # the abiding probability across positions
            data[2,:] ./= num_trajs # the surfing probability across positions
            data[3,:] ./= num_surfs # the average surfing Gillespie-time across positions
            data[4,:] ./= num_surfs # the average surfing simulation-time across positions
            data[5,:] ./= num_trajs # the probability that the mutant didn't survive at all across positions

            if to_plot == "surfabide"
                #### BEGIN PLOT SURF AND ABIDING ####
                
                z=plot(grid=false,
                border=:box,
                size=(900,600),
                xlims=(0,150),
                yticks=LinRange(0.0,1.0,11),
                ytickfontsize=16,
                ylabel=L"\textrm{Probability}",
                ylabelfontsize=18,
                leftmargin=5mm,
                xlabel=L"\textrm{Initiation}\;\textrm{Site}\;(x)",
                xlabelfontsize=18,
                xtickfontsize=16,
                margin=5mm)

                if rm < 0.6
                    plot!(z,legend=:topright,legendfontsize=18)
                else
                    plot!(z,legend=:right,legendfontsize=18)
                end

                plot!(z,init_wave[add_to_bulk+1:end],label="",color=:black,linewidth=2)
                
                half_loc = argmin(abs.(init_wave.-0.5))

                scatter!(z,data[2,:],
                         label=L"s_\textrm{sim}(x)",
                         markercolor=:pink,
                         markerstrokecolor=:red,
                         markeralpha=0.5,
                         markersize=7.0,
                         markershape=:circle)
                scatter!(z,data[1,:],
                         label=L"a_\textrm{sim}(x)",
                         markercolor=:yellow,
                         markerstrokecolor=:gold,
                         markeralpha=0.5,
                         markersize=7.0,
                         markershape=:rect)
                scatter!(z,data[1,:].+data[2,:],
                         label=L"u_\textrm{sim}(x)",
                         markercolor=:lightblue,
                         markerstrokecolor=:blue,
                         markeralpha=0.5,
                         markersize=7.0,
                         markershape=:utriangle)
                annotate!(z,[half_loc + 15],[0.55],text(L"\langle w \rangle_{\textrm{init}}(x)",:black,20))
                annotate!(z,[half_loc + 15],[0.45],text(L"r_w="*latexstring("$(rw)"),:black,20))
                        
                if save_fig == true
                    savefig(z,fig_dest*"Surf_vs_Abide_rw_$(rw)_rm_$(rm)_Pwmm_$(Pwmm)_Pmwm_$(Pmwm).png")
                end
                display(z)

                #### END PLOT SURF AND ABIDING ####
            elseif to_plot == "surv"
                #### BEGIN PLOT SURV AS NUMERICAL ####

                if counter == 1
                    plot!(p,grid=false,
                    border=:box,
                    # legend=:outerright,
                    size=(900,600),
                    xlims=(0,150),
                    yticks=LinRange(0.0,1.0,11),
                    ytickfontsize=16,
                    ylabel=L"\textrm{Survival}\;\textrm{Probability}",
                    ylabelfontsize=18,
                    leftmargin=5mm,
                    xlabel=L"\textrm{Initiation}\;\textrm{Site}\;(x)",
                    xlabelfontsize=18,
                    xtickfontsize=16,
                    margin=5mm, 
                    legend=true, 
                    legendfontsize=18)

                    plot!(p, xs, ws,
                          label="",
                          color=:black,
                          linewidth=2)
                    
                    annotate!(p,[36-(0.9/rw)^1.2],[0.9],text(L"\langle w \rangle_{\textrm{init}}(x)",:black,20))
                    annotate!(p,[36-(0.9/rw)^1.2],[0.8],text(L"r_w="*latexstring("$(rw)"),:black,20))

                end

                # Get BOTH boundary conditions
                Bulk_BC = BC(Pwmm,Pmwm,1e6)
                Front_BC = BC(Pvmm,Pmvm,1e6)
                println(Bulk_BC," ",Front_BC)

                # Below, the numerical calculations of the survival ode are performed.
                # If the Bulk BC is small, I perform the numerical integration from right to left, 
                # because 0 Bulk BC gives a constant solution of zero. The mathematical formulas
                # for what I am doing here are in the paper: https://www.biorxiv.org/content/10.1101/2024.12.14.628506v1

                F = (Pmvm + Pvmm) .* vs .+ (Pmwm + Pwmm) .* ws 
                J = ((ReLu(Pmvm) + ReLu(Pvmm)) .* vs .+ (ReLu(Pmwm) + ReLu(Pwmm)) .* ws) 

                Us = zeros(length(xs))
                Us[1] = Bulk_BC
                Us[2] = Bulk_BC
                Us[3] = Bulk_BC
                Us[4] = Bulk_BC
                Us[end] = Front_BC
                Us[end-1] = Front_BC
                Us[end-2] = Front_BC
                Us[end-3] = Front_BC

                if Bulk_BC >= 1.0/K

                    g = 1.0 ./(2.0 .*D .+h .* vw_front)
                    @inbounds for j in 3:length(xs)
                        Us[j] =  g[j-1]*(
                                    4.0*D*Us[j-1] 
                                + 2.0*F[j-1]*(h^2)*Us[j-1] 
                                - 2.0*(h^2)*J[j-1]*(Us[j-1]^2) 
                                - 2.0*D*Us[j-2] 
                                + h*Us[j-2]*vw_front[j-1] )
                    end  

                elseif Bulk_BC < 1.0/K

                    g = 1.0 ./(2.0 .*D .+h .* vw_front)
                    @inbounds for j in length(xs)-2:-1:1
                        Us[j] =  g[j+1]*(
                                    4.0*D*Us[j+1] 
                                + 2.0*F[j+1]*(h^2)*Us[j+1] 
                                - 2.0*(h^2)*J[j+1]*(Us[j+1]^2) 
                                - 2.0*D*Us[j+2] 
                                + h*Us[j+2]*vw_front[j+1] )
                    end   
                end

                # Total survival is the sum of the surfing and abiding probabilities, alternatively one could write 1.0 .- data[5,:]
                surv_data = data[1,:].+data[2,:]
                scatter!(p,surv_data,
                         label=L"r_m="*latexstring("$(rm)"),
                         markercolor=colors[counter],
                         markerstrokecolor=strokecolors[counter],
                         markeralpha=0.5,
                         markersize=7.0,
                         markershape=markershapes[counter])

                plot!(p, xs, Us,
                      label="", 
                      color=strokecolors[counter], 
                      legend=:topright)        
                
                counter += 1 
                #### END PLOT SURV AS NUMERICAL ####
            elseif to_plot == "times"
                #### BEGIN PLOT SURF TIMES ####

                scatter!(q,data[3,:],
                         label=L"r_m="*latexstring("$(rm)"),
                         markercolor=colors[counter],
                         markerstrokecolor=strokecolors[counter],
                         markeralpha=0.5,
                         markersize=5.0,
                         markershape=markershapes[counter])

                endtimes[dirnum]=data[3,end]

                if maximum(data[3,.!isnan.(data[3,:])]) > max_time
                    max_time = maximum(data[3,.!isnan.(data[3,:])])
                end
                if minimum(data[3,.!isnan.(data[3,:])]) < min_time
                    min_time = minimum(data[3,.!isnan.(data[3,:])])
                end
                counter += 1 
                #### END PLOT SURF TIMES ####
            end
        end

        if to_plot == "surfabide"
            println(to_plot)
        elseif to_plot == "surv"
            println(to_plot)

            plot!(p,ylims=(-0.01,1.01))
            display(p)
            if save_fig == true
                savefig(p,"Figs/paper_rw_$(rw)_Pwmm_$(Pwmm)_Pmwm_$(Pmwm)_no_num_psurv_multiple_rms.png")
            end
        elseif to_plot == "times"
            println(to_plot)

            plot!(q,grid=true,
            border=:box,
            # legend=false,
            size=(900,600),
            ytickfontsize=16,
            ylabel=L"\textrm{Surfing}\;\textrm{Time}",
            ylabelfontsize=18,
            leftmargin=5mm,
            xlabel=L"\textrm{Initiation}\;\textrm{Site}\;(x)",
            xlabelfontsize=18,
            xtickfontsize=16,
            bottommargin=5mm,
            rightmargin=5mm, 
            ylims=(min_time*0.99,max_time*1.01),
            legend=true,
            legendfontsize=18)

            plot!(q, xs, (max_time-min_time) .* ws .+ min_time,
                  label="",
                  color=:black,
                  linewidth=2)

            if figcounter in [2 3 4 6 7 9 10 11 12]
                annotate!(q,[36-(0.9/rw)^1.2],[0.55*(max_time-min_time)+min_time],text(L"\tau_{max}\langle w \rangle_{\textrm{init}}(x)",:black,20))
                annotate!(q,[36-(0.9/rw)^1.2],[0.45*(max_time-min_time)+min_time],text(L"r_w="*latexstring("$(rw)"),:black,20))
            elseif figcounter == 5
                annotate!(q,[100-(0.9/rw)^1.2],[0.55*(max_time-min_time)+min_time],text(L"\tau_{max}\langle w \rangle_{\textrm{init}}(x)",:black,20))
                annotate!(q,[100-(0.9/rw)^1.2],[0.45*(max_time-min_time)+min_time],text(L"r_w="*latexstring("$(rw)"),:black,20))
            elseif figcounter in [1 8]
                annotate!(q,[36.5-(0.9/rw)^1.2],[0.57*(max_time-min_time)+min_time],text(L"\tau_{max}\langle w \rangle_{\textrm{init}}(x)",:black,20))
                annotate!(q,[36.5-(0.9/rw)^1.2],[0.47*(max_time-min_time)+min_time],text(L"r_w="*latexstring("$(rw)"),:black,20))
            end

            display(q)
            if save_fig == true
                savefig(q,"Figs/paper_rw_$(rw)_Pwmm_$(Pwmm)_Pmwm_$(Pmwm)_times_multiple_rms.png")
            end
        end
    end

end