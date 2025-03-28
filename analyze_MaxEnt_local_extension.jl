"""
This script was used to analyze the local extension ahead/behind the replication forks,
as well as averaged over different genomic regions.

We fit exponential functions to the averaged local extensions, and find exponential decay
after replication.

The script should be run **after** process_MaxEnt.jl, since that's where the local extension profiles are calculated.
"""


include("src/Analyze_ME.jl")
include("src/Plots.jl")

data_file="Stats/MaxEnt/Full_MaxEnt_correct_oridist_fast.h5"
stat_dir="Stats/MaxEnt_extension_time_scales/"
plot_directory="Plots/MaxEnt/Strain_timescales/"
mkpath(plot_directory)
mkpath(stat_dir)

#Load data
replication_speed=2.6 #bins/min

sampled_ts=[0.,10.,30.,45.,60.,75.]
region_widths=5:5:25

extensions_abs=h5read(data_file, "mean_local_extensions_varying_n")[:,:,1]
extensions=relative_extensions(extensions_abs) #relative to t=0
forks=Int.(h5read(data_file, "mean_forks"))

plot_local_extensions_ahead_of_forks(forks, extensions, sampled_ts, plot_directory, region_widths)

plot_extension_decay_after_replication(sampled_ts,extensions, stat_dir,plot_directory, region_widths)

make_latex_table_timescales(stat_dir, plot_directory)
