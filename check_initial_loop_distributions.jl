"""
This script was used to check how many generations of 1D loop-extruder simulations 
where the final loop-extruder positions are inherited by the next generation are 
needed to reach a self-consistent distribution.

In other words, each generation inherits the final loop-extruder positions of the last,
and we require a large enough generation count that the final distributions of loop-extruders
are the same for consecutive generations.

For actual simulations, we exclude these first few generations before a consistent
state has been reached.
"""


include("src/Plots.jl")
include("src/statistics_functions.jl")


function moving_mean(data, window_size)
    smoothed_data = similar(data)
    half_window = div(window_size, 2)
    
    for i in 1:length(data)
        start_index = max(1, i - half_window)
        end_index = min(length(data), i + half_window)
        smoothed_data[i] = mean(data[start_index:end_index])
    end
    
    return smoothed_data
end

plots_dir="Plots/1D_tests_MolecularDynamics/"
simulation_data_1d_dir="Results_1D/Replicating_caulobacter_generational_N_1620_M_30_loopsize_480.00000000000006_ter_size_40_ter_strength_100_GPU_0/"

files=filter(contains("/gen_"), readdir(simulation_data_1d_dir, join=true))
N=parse_after(simulation_data_1d_dir, "N_")
num_gens=80
num_sims=0
initial_chipseq=zeros((N, num_gens))
loop_size_bins=1:25:2025
initial_loop_size_hist=zeros(length(loop_size_bins)-1, num_gens)

for f in files
    h5open(f) do file
        l_sites=read(file, "l_sites")[:,1].+1 #1 based indexing in Julia
        r_sites=read(file, "r_sites")[:,1].+1

        #exclude unbound
        smcs=hcat(l_sites[:,t], r_sites[:,t])
        bound=smcs[:,1].>0
        smcs=smcs[bound,:]
        smcs=permutedims(smcs)
        
        gen_num=parse_after(f, "gen_")+1
        sim_num=parse_after(f, "simulation_")+1
        chipseq=fetch_chipseq(smcs, N)
        initial_chipseq[:,gen_num].+=chipseq

        loop_size_counts=fetch_loop_size_distribution(smcs, loop_size_bins,N)
        initial_loop_size_hist[:,gen_num].+=loop_size_counts

        global num_sims=max(num_sims, sim_num)
    end
end

#normalize
initial_chipseq./=num_sims
initial_loop_size_hist./=num_sims

#do a plot over time
plot(xlabel="Genomic position [Mb]", ylabel="Initial ChIP-seq signal", xticks=xticks_Mb(N))
for gen in 1:10:num_gens
    #smoothen the curve a bit to get rid of noise
    chip_seq_curve=moving_mean(initial_chipseq[:,gen],10)

    plot!(circshift(chip_seq_curve, Int(N/2)), label="Gen $gen", alpha=1-0.05*(gen-1)/10)
end
savefig(plots_dir*"initial_chipseq_curves_different_generations.pdf")

plot(xlabel="Loop size [kb]", ylabel="Frequency")
for gen in 1:10:num_gens
    smoothened_distr=moving_mean(initial_loop_size_hist[:,gen], 10)
    plot!(loop_size_bins[1:end-1], smoothened_distr, label="Gen $gen", alpha=1-0.05*(gen-1)/10)
end
savefig(plots_dir*"initial_loop_size_distributions_different_generations.pdf")
