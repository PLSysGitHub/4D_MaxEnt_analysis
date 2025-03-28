"""
This script was used to check how the fraction of small loops varies in time.

We load data from 1D simulations, and analyze the loop size distributions. 

Additionally, the analytical predictions for the fraction of small loops are calculated.
"""


include("src/Plots.jl")
include("src/statistics_functions.jl")

using OrdinaryDiffEq

M_at_t(t, M0, vr, N) = M0 + M0 * 2 * vr * t / N

#Numerically solve the ODE that describes the fraction of small loops
function solve_nsmall(n0, D_vals, times, τ, M0, vr, N, Δx, k_grow)    
    function dn_small_dt(n, p, t)
        M=M_at_t(t, M0, vr, N)
        D_t = D_vals[round(Int, vr*t*N/4050)+1] #nearest lattice point
        term1 = (1 - n) / τ
        term2 = 2*vr*M0/N/M
        term3 = (D_t *vr)*(1-n) / Δx / M
        term4 = -n * k_grow
        
        return term1 + term2 + term3 + term4
    end
    
    tspan = (minimum(times), maximum(times))
    prob = ODEProblem(dn_small_dt, n0, tspan)
    sol = solve(prob, Tsit5(), saveat=times)
    
    return sol.u
end

#Numerically solve the ODE that describes the fraction of small loops, without loop-extruder unbinding at the forks
function solve_nsmall_no_fork_unbind(n0, times, τ, M0, vr, N, k_grow)    
    function dn_small_dt_no_fork(n, p, t)
        M=M_at_t(t, M0, vr, N)
        term1 = (1 - n) / τ
        term2 = 2*vr*M0/N/M
        term4 = -n * k_grow
        return term1 + term2 + term4
    end
    
    tspan = (minimum(times), maximum(times))
    prob = ODEProblem(dn_small_dt_no_fork, n0, tspan)
    sol = solve(prob, Tsit5(), saveat=times)
    
    return sol.u
end

#Find the indices of the correct time-points
function fetch_time_indices(ts, t_targets)
    t_indices=[findlast(ts.<=t) for t in t_targets]
    t_indices[1]=1
    t_indices[end]=length(ts)
    return t_indices
end

plots_dir="Plots/1D_tests_MolecularDynamics/"
mkpath(plots_dir)

simulation_data_1D_parent="Results_1D/"
times=[0, 10, 30, 45, 60, 75, 80] #the 80 time-point is the last one; actually corresponds to initial distribution of loop-extruders due to inheritance
loop_size_bins=0:10:4050/2
nt=length(times)


#Loop over all except the generational simulations; those were used to check how many generations are needed before the initial loop-extruder positions are sampled correctly
simulation_data_1d_dirs=filter(!contains("generational"), readdir(simulation_data_1D_parent, join=true))
simulation_data_1d_dirs=filter(isdir, simulation_data_1d_dirs)

#Prefixes for saving plots
file_prefixes=basename.(simulation_data_1d_dirs)
file_prefixes=replace.(file_prefixes, "_GPU_5"=>"")

#Loop over all simulations
for (simulation_data_1d_dir, file_prefix) in zip(simulation_data_1d_dirs, file_prefixes)
    files=filter(contains("/simulation_"), readdir(simulation_data_1d_dir, join=true))
    N=parse_after(simulation_data_1d_dir, "N_")
    num_sims=0
    chipseq=zeros((N,nt))
    loop_size_hist=zeros(length(loop_size_bins)-1, nt)

    for f in files
        try
            h5open(f) do file
                ts=read(file, "ts")
                t_indices=fetch_time_indices(ts, times)

                l_sites=read(file, "l_sites")[:,t_indices].+1
                r_sites=read(file, "r_sites")[:,t_indices].+1

                for t in 1:nt
                    smcs=hcat(l_sites[:,t], r_sites[:,t])
                    bound=smcs[:,1].>0

                    smcs=smcs[bound,:]
                    smcs=permutedims(smcs)

                    chipseq_t=fetch_chipseq(smcs, N)
                    chipseq[:,t].+=chipseq_t

                    loop_size_counts=fetch_loop_size_distribution(smcs, loop_size_bins,N)
                    loop_size_hist[:,t].+=loop_size_counts
                end
                num_sims+=1
            end
        catch
            println("couldnt open $f")
        end
    end

    #normalize
    chipseq./=num_sims
    for i in 1:nt
        loop_size_hist[:,i]./=sum(loop_size_hist[:,i])
    end

    #do a plot over time
    plot(xlabel="Genomic position [Mb]", ylabel="ChIP-seq signal", xticks=xticks_Mb(N))
    for t in 1:nt-1
        #smoothen the curve a bit to get rid of noise
        chip_seq_curve=chipseq[:,t]

        plot!(circshift(chip_seq_curve, Int(N/2)), label="t $(times[t])", alpha=1-0.05*(t-1)/10, color=t)
        vline!([N/2+23*N/4050*times[t], N/2-23*times[t]*N/4050], color=t)
    end
    savefig(plots_dir*"chipseq_time_$(file_prefix).pdf")

    plot_loop_size_histograms(loop_size_bins, loop_size_hist, times)
    savefig(plots_dir*"loop_size_distributions_$(file_prefix).pdf")

    #Calculate the relative density of small loops, together with analytical estimate
    M0=parse_after(file_prefix, "_M_")
    n0=sum(loop_size_hist[1:4,1])

    vr=4050/2/85 #in kb/min
    Δx=4050/N #in kb
    le_speed=18 #kb/min

    #loop growth rates were roughly fitted for simulations with by-passing of the forks
    if contains(file_prefix, "t_stall_0.0")
        #no stalling; faster loop growth
        k_grow=1.95*le_speed/Δx/17 # in 1/min
    elseif contains(file_prefix, "ori_aff_14580")
        #90% chance load at ori; stall upon collision
        #more collisions, so slower loop growth
        k_grow=0.91*le_speed/Δx/17 # in 1/min
    else
        #50% chance load at ori; stall upon collision
        k_grow=le_speed/Δx/17 # in 1/min
    end
    τ=parse_after(file_prefix, "_loopsize_")/2*4050/N/le_speed
    ts=0:2:times[end-1]

    println("parameters are: n0=$n0, M0=$M0, vr=$vr, Δx=$Δx, k_grow=$k_grow, τ=$τ")

    if contains(file_prefix, "fork_bypass")
        n_small=solve_nsmall_no_fork_unbind(n0, ts, τ, M0, vr, N, k_grow)
    else
        D_vals=chipseq[:,end] #last time point inherited by next generation; so this is the initial condition
        D_vals.*=2*M0/sum(D_vals) #normalize; total number should be 2*M0, 2 for two legs
        n_small=solve_nsmall(n0, D_vals, ts, τ, M0, vr, N, Δx, k_grow)
    end

    n_small./=n_small[1]

    plot_dens_small_loops(loop_size_hist[:,1:end-1], loop_size_bins, times[1:end-1], 40, label="Simulation", size=(500,400), legend=:topright)
    plot!(ts, n_small, label="Analytical", color=:red, ylims=(0.9,3))
    savefig(plots_dir*"rel_dens_small_loops_$(file_prefix).pdf")

    #Whisker plots of loop sizes
    plot_loop_size_distributions(loop_size_hist[:,1:end-1], times[1:end-1])
    savefig(plots_dir*"whiskers_loops_$(file_prefix).pdf")
end
