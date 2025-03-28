"""
Given raw Hi-C counts from different time-points, estimate the fork positions by analyzing the relative counts compared to t=0.
Make plots that illustrate the steps, as well as the inferred fork positions.

This pipeline follows Le et al. 2012
"""

include("src/Plots.jl")
include("src/Analyze_ME.jl")

N=405
hic_folder="Input_sims/crescentus_Hi-c/Data_Le_2013/"
hic_files=filter(contains("before_normalization"),readdir(hic_folder, join=true))

plot_folder="Plots/Fork_infer/"
mkpath(plot_folder)

hic0=readdlm(hic_files[findfirst(contains.(hic_files, "_0min_overlap_before_normalization"))])
colsums=sum(hic0, dims=2)
colsums.*=N/sum(colsums)

#Parse the hic file names to get the time points
times=map(x->parse(Int, split(split(x,"min")[1], '_')[end]), hic_files)

#Set the filter width for the gaussian kernel
filter_width=50;
σ=35 #standard deviation of the gaussian kernel

inferred_forks_all=zeros(2,length(times))

for (index,t) in enumerate(times)
    hic_ind=findfirst(contains.(hic_files, "_$(t)min_overlap_before_normalization"))
    hic=readdlm(hic_files[hic_ind])

    rel_col_sums=vec(sum(hic, dims=2)) #sum counts for each column
    rel_col_sums.*=N/sum(rel_col_sums) #normalize
    rel_col_sums./=colsums #divide by the 0min counts

    ders, smoothed=smooth_and_differentiate_periodic(rel_col_sums, filter_width,σ)
    found_forks= t==0 ? [1,N] : find_fork_indices(ders)
    inferred_forks_all[:,index].=found_forks

    #Save plots of (smoothened) curves
    plot_relative_counts(rel_col_sums, smoothed, found_forks, t, plot_folder)
    plot_ders(ders, found_forks, t, plot_folder)
end

#Do a linear fit
inferred_forks_all./=100 #Mb
used_inds=[3,4,5] #the indices of the time points used for the fit

fit_right=linear_fit(times[used_inds], inferred_forks_all[1,used_inds])
fit_left=linear_fit(times[used_inds], N/100 .-inferred_forks_all[2,used_inds]) #left fork is going "backwards"

plot_inferred_forks(times, inferred_forks_all, fit_right, fit_left, plot_folder)
