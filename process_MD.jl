"""
This script calculates statistics for molecular dynamics simulations.

Raw data directory should contain subdirectories corresponding to different simulations

For each subdirectory, we calculate averaged statistics over all blocks_*.h5 files, which should contain polymer configurations over time.

"""


include("src/analyze_md.jl")
include("src/validation.jl")

#Directories
raw_data_directory="MD_Results_3D/Combined_replicating/"
stats_directory="Stats/MolecularDynamics/"
validation_out_dir="Validation_stats/" #Where to output means and stds for validation data

#Options
skip_done=false #skip a directory if statistics have already been calculated
ME_orient=true #orient so that strand 1 has ori further from ter
fork_orient=false #orient so that fork closer to pole 1 is fork 1
times=[0, 10, 30, 45, 60, 75] #times in minutes to calculate statistics for

#Basic statistics
calc_fast_stats(times, raw_data_directory, stats_directory, skip_done=skip_done, ME_orient=ME_orient, fork_orient=fork_orient)
calc_hic(times, raw_data_directory, stats_directory, ME_orient=ME_orient, skip_done=skip_done)


#Spring lengths
calc_spring_lengths(raw_data_directory, stats_directory, skip_done=skip_done)

