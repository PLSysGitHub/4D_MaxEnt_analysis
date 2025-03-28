"""
This file was used to check whether the diffusion of the chromosome segments in the simulations
occurrs at the expected rate.

The diffusion constant in the presence of loop-extruders was of the same order of magnitude
as that found by Weber et.al. 2012.
"""

include("src/analyze_md.jl")
include("src/Plots.jl")

fig_type=".pdf" #can also save png or other formats
nt=500 #total 500 s, sample every 1 s

data_dir="MD_Results_3D/Short_times_start_top_monomer_0/"
out_dir="Stats/MD_Short_times/"
plots_dir="Plots/MolecularDynamics/Check_time_scales/"
skip_done=false #whether or not existing stat files should be overwritten

mkpath(out_dir)
mkpath(plots_dir)

#For each directory, calculate the mean squared distance of the origin and all monomers
for subdir in subdirs(data_dir)
    println(subdir)
    N=parse_after(subdir, "N_")
    tracked_loci=[1, floor(Int,N/4), Int(N/2), floor(Int,3*N/4)]
    monomer_size=parse_after(subdir, "_b_")/1000
    initial_stage_files=readdir(subdir, join=true)
    initial_stage_files=filter(contains("blocks"), initial_stage_files)
    L_0=parse_after(subdir, "_L_")*monomer_size
    stats_file=replace(subdir, data_dir=>out_dir)*".h5"

    if !skip_done || !isfile(stats_file)
        mean_rsq, error_rsq=fetch_msd_trajectories(initial_stage_files,monomer_size, tracked_loci)
        h5open(stats_file, "w") do f
            @write f mean_rsq
            @write f error_rsq
        end
    else
        mean_rsq=h5read(stats_file, "mean_rsq")
        error_rsq=h5read(stats_file, "error_rsq")
    end

    ts=findall(mean_rsq[1,:].>0)
    mean_rsq=mean_rsq[:,ts]

    out_name=replace(subdir,data_dir=>plots_dir)*"diffusion_early_all_loci"*fig_type

    plot_mean_rsq(ts, mean_rsq)
    savefig(out_name)
end
