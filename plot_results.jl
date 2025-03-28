"""
Main plotting script. For both MD and 4D-MaxEnt simulations, plot various statistics.

Also does the validation plots.
"""



include("src/Plots.jl")

times=[0, 10, 30, 45, 60, 75] #timepoints in minutes
skip_done=false #skip a directory if plots already made

#First do molecular dynamics plots
stats_directory="Stats/MolecularDynamics/"
plots_directory="Plots/MolecularDynamics/"
make_fast_plots(times, stats_directory, plots_directory, skip_done=skip_done)
make_spring_length_histograms(stats_directory, plots_directory)
make_hic_plots(times, stats_directory, plots_directory; skip_done=skip_done, fig_type=".png")

#Then do MaxEnt plots
stats_directory="Stats/MaxEnt/"
plots_directory="Plots/MaxEnt/"
make_spring_length_histograms(stats_directory, plots_directory)
make_fast_plots(times, stats_directory, plots_directory, skip_done=skip_done)

#make plots of the effective energies
energies_parent="Converged_energies"

for directory_in in subdirs(energies_parent)
   if !contains(directory_in, "Ori constraint only")
       directory_out="Plots/MaxEnt/Energies/".*basename(directory_in).*"/" #where to save the output
       mkpath(directory_out)

       visualize_energies(directory_in, directory_out)
   end
end

#Validation plots
validation_dir="Validation_stats/"
experimental_dir="Validation_stats/Experimental_validation_data/"
validation_plot_directory="Plots/Validation/"

oclocks=[1,2,3,4,6,8,9,11] #for reading in data
locus_positions = [0.43,0.60,0.96,1.50,2.03,2.67,3.03,3.65] #genomic position in Mb
replication_speed = 0.026
genomic_length = 4.017 #in Mb
locus_replication_times = [min(x,genomic_length-x)/replication_speed for x in locus_positions] #expected time of replication

make_validation_plots(times,oclocks, locus_replication_times, validation_dir, experimental_dir, validation_plot_directory)
