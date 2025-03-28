"""
This script processes 4D-MaxEnt simulation data and calculates statistics.

Raw data directory should contain subdirectories corresponding to different models.

For each subdirectory, we first save configurations in .h5 files, then calculate statistics over these configurations.
"""

include("src/preprocess_ME_data.jl")
include("src/Analyze_ME.jl")
include("src/validation.jl")


raw_data_directory="Raw_ME_configs/"
processed_data_directory="Processed_ME_configs/"
stats_directory="Stats/MaxEnt/"
validation_out_dir="Validation_stats/" #Where to output means and stds for validation data
validation_exp_dir="Experimental_validation_data/" #should contain csv files "[oclock]oclock_ori.csv" with experimental data for locus positions

#Options
skip_done=false #skip a directory if statistics have already been calculated

# MaxEnt simulation parameters
times=[0, 10, 30, 45, 60, 75] #minutes
Ls_sim=[18, 19, 22, 24, 28, 31].+2*3.7 #simulation units
N=405 #number of Hi-C bins or MaxEnt monomers
monomer_size=0.088 #in μm
Ls=Ls_sim.*monomer_size #convert to μm

#####################################
#Run the pipeline
#compress_MaxEnt_files(raw_data_directory,processed_data_directory) #skips if files already exist

calc_fast_ME_stats(times, Ls_sim, processed_data_directory, stats_directory, N=N, skip_done=skip_done) # skips if files already exist
calc_spring_lengths(times, processed_data_directory,stats_directory; skip_done=skip_done, bin_width=44)

calc_hic(times, processed_data_directory, stats_directory, N=N, skip_done=true) # skips if files already exist


#####################################
#Validation against experimental data

oclocks=[1,2,3,4,6,8,9,11] #for reading in data
locus_positions = [0.43,0.60,0.96,1.50,2.03,2.67,3.03,3.65] #genomic position in Mb
replication_speed = 0.026
genomic_length = 4.017 #in Mb
locus_degrees=[round(Int,x*360/genomic_length) for x in locus_positions] #genomic position in degrees
locus_replication_times = [min(x,genomic_length-x)/replication_speed for x in locus_positions] #expected time of replication
exp_replication_times=[30,30,45,60,200,60,45,30] #the times when loci were seem to be replicated in the experiment
threshold_for_conditioning=0.320/monomer_size #if separation below this, don't count. 320 nm is the resolution of the experiment

calc_validation_stats_experiment(oclocks,exp_replication_times,validation_exp_dir,validation_out_dir)
calc_validation_stats_MaxEnt(times, Ls_sim, oclocks, locus_positions, locus_replication_times, genomic_length,N,threshold_for_conditioning, processed_data_directory, validation_out_dir)
