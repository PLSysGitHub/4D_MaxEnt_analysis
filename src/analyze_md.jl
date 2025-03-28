"""
This file runs a pipeline that calculates statistics given a folder of loop-extrusion simulation results

Takes a path to input data as an argument.

Calculates:

- Mean z-positions and histograms
- Mean fork positions
- Mean local extensions
- Hi-C maps
- Mean chipseq and loop size distributions [if there's smcs in the sims]

averaged over all samples, at roughly the times requested in 'times'

The statistics are saved in .h5 files
"""

include("statistics_functions.jl")
using Base.Threads

"""
    calc_fast_stats(times, data_parent_dir, out_parent_dir; skip_done=true, ME_orient=false, fork_orient=false, ns=[1,2,3,4,5,6], total_kb=4050, loop_bin_width=10)

Calculates statistics for a folder of loop-extrusion simulation results

# Arguments
- `times::Array{Float64}`: Times at which to calculate statistics
- `data_parent_dir::String`: Path to the folder containing the simulation results
- `out_parent_dir::String`: Path to the folder where the statistics will be saved
- `skip_done::Bool=true`: If true, skips directories for which statistics have already been calculated
- `ME_orient::Bool=false`: If true, the ME orientation is used; strand 1 is ter-distal one
- `fork_orient::Bool=false`: If true, the fork closer to the ter-distal pole is labeled as fork 1
- `ns::Array{Int}=[1,2,3,4,5,6]`: Local extension lengths to calculate
- `total_kb::Int=4050`: Total length of the chromosome in kb
- `loop_bin_width::Int=10`: Width of the loop size bins
"""
function calc_fast_stats(times, data_parent_dir,out_parent_dir; skip_done=true, ME_orient=false, fork_orient=false, ns=[1,2,3,4,5,6], total_kb=4050, loop_bin_width=10)

    subdirs=readdir(data_parent_dir, join=true)
    subdirs=subdirs[isdir.(subdirs)]
    mkpath(out_parent_dir)

    loop_size_bins=0:loop_bin_width:total_kb/2

    for dir in subdirs
        out_file_name=replace(dir, data_parent_dir=>out_parent_dir)*".h5"
        if contains(dir, "GPU") || (skip_done && isfile(out_file_name))
            println("Found previous results for $dir, skipping...")
        else
            println("Calculating statistics for $dir")
            num_samps=zeros(Int, length(times)) #should be same, but if a file has been deleted/corrupted, might slightly differ

            N=parse_after(dir, "N_")
            monomer_size=parse_after(dir, "_b_")/1000

            #For comparing to experiments, we want distance corresponding to 10 kb
            incr=round(Int, N/405)

            #Initialize output arrays
            mean_zs=zeros(2*N, length(times))
            mean_zs_resolution_lim=zeros(2*N, length(times))
            counts_zs_resolution_limit=zeros(Int, 2*N, length(times))
            mean_radial_pos=zeros(2*N,length(times))
            mean_forks=zeros(2, length(times))
            max_h=time_to_height(times[end])/monomer_size #in monomer lengths
            hist_bins=-max_h/2:0.088/monomer_size:max_h/2 #the step matches MaxEnt model
            z_histograms=zeros(length(hist_bins)-1, 3, length(times))

            if !contains(dir, "No_smcs") && !contains(dir, "nosmcs")
                loop_size_distr=zeros(length(loop_size_bins)-1,length(times))
                mean_chipseq=zeros(N, length(times))
            end
            mean_local_extensions_different_n=zeros(2*N,length(times), length(ns))
            #Loop over files, each containing a single time trajectory
            for filename in filter(contains(".h5"), readdir(dir,join=true))
                h5open(filename, "r") do file
                    list_inds=get_indices_at_times(file,times)

                    #Get all statistics at required times
                    for t in 1:length(times)
                        for i in list_inds[t]
                            num_samps[t]+=1
                            pos, fork=fetch_pos_fork_at_ind(file, i, ME_orient=ME_orient, fork_orient=fork_orient)
                            mean_forks[:,t].+=fork
                            mean_zs[:,t].+=fetch_zs(pos,fork)
                            mean_radial_pos[:,t].+=fetch_radial_positions(pos, fork)
                            zs,resolved=fetch_zs_resolution_limited(pos, fork)
                            mean_zs_resolution_lim[resolved,t].+=zs[resolved]
                            counts_zs_resolution_limit[resolved,t].+=1

                            if !contains(dir, "No_smcs") && !contains(dir, "nosmcs")
                                SMCs=fetch_smcs_at_ind(file, i)
                                mean_chipseq[:,t].+=fetch_chipseq(SMCs,N)
                                loop_size_distr[:,t].+=fetch_loop_size_distribution(SMCs,loop_size_bins,N)
                            end

                            for (j,n) in enumerate(ns)
                                mean_local_extensions_different_n[:,t,j].+=fetch_local_extensions(pos,fork,n*incr)
                            end

                            z1,z2,z3=fetch_z_histograms(pos, fork, hist_bins)
                            z_histograms[:,1,t].+=z1
                            z_histograms[:,2,t].+=z2
                            z_histograms[:,3,t].+=z3
                        end
                    end
                end
            end

            #divide by number of samples to get average. Note that the number might vary per time point
            for t in 1:length(times)
                mean_zs[:,t]./=num_samps[t]
                mean_radial_pos[:,t]./=num_samps[t]
                mean_forks[:,t]./=num_samps[t]
                z_histograms[:,1,t]./=num_samps[t]
                z_histograms[:,2,t]./=num_samps[t]
                z_histograms[:,3,t]./=num_samps[t]
            end

            mean_zs_resolution_lim./=counts_zs_resolution_limit
            #save all the statistics to a .h5 file
            hist_bins=collect(hist_bins) #need these to be vectors for h5
            loop_size_bins=collect(loop_size_bins)

            out_file_name=replace(dir, data_parent_dir=>out_parent_dir)*"_fast.h5"
            h5open(out_file_name, "w") do file
                @write file num_samps
                @write file mean_zs
                @write file mean_zs_resolution_lim
                @write file mean_forks
                @write file mean_radial_pos
                if !contains(dir, "No_smcs") && !contains(dir, "nosmcs")
                    #normalize the matrices
                    for t in 1:length(times)
                        mean_chipseq[:,t]./=num_samps[t]
                        norm_const=sum(loop_size_distr[:,t])
                        loop_size_distr[:,t]./=norm_const  #probability density
                    end

                    @write file loop_size_distr
                    @write file loop_size_bins
                    @write file mean_chipseq

                    #also calculate the mean loop sizes
                    mean_loop_sizes=[sum(loop_size_distr[:,t].*loop_size_bins[1:end-1])/sum(loop_size_distr[:,t]) for t in 1:length(times)]
                    @write file mean_loop_sizes
                end
                @write file mean_local_extensions_different_n
                @write file ns
                @write file hist_bins
                @write file z_histograms
            end
        end
    end
end


"""
    calc_spring_lengths(data_parent_dir, out_parent_dir; skip_done=false, bin_width=44)

Calculates spring lengths statistics for a folder of loop-extrusion simulation results

# Arguments
- `data_parent_dir::String`: Path to the folder containing the simulation results
- `out_parent_dir::String`: Path to the folder where the statistics will be saved
- `skip_done::Bool=false`: If true, skips directories for which statistics have already been calculated
- `bin_width::Int=44`: Width of the spring length bins
"""
function calc_spring_lengths(data_parent_dir,out_parent_dir; skip_done=false, bin_width=44)
    subdirs=readdir(data_parent_dir, join=true)
    subdirs=subdirs[isdir.(subdirs)]
    mkpath(out_parent_dir)

    for dir in subdirs 
        monomer_size=parse_after(dir, "_b_")
        N=parse_after(dir, "N_")

        #For comparing to experiments, we want distance corresponding to 10 kb
        incr=round(Int, N/405)

        out_file_name=replace(dir, data_parent_dir=>out_parent_dir)*".h5"
        if (skip_done && isfile(out_file_name))
            println("Found previous results for $dir, skipping...")
        else
            println("Calculating spring length statistics for $dir")
            #Initialize output arrays
            d_bins=collect(0.0:bin_width:500) #in the simulation units where monomer size is 1
            spring_lengths=zeros(Int, length(d_bins))
            num_samples=0

            #Loop over files, each containing a single time trajectory
            for filename in filter(contains("blocks"), readdir(dir,join=true))
                h5open(filename, "r") do file
                    try
                        list_inds=keys(file)
                        for list_ind in list_inds[end-5:end]
                            if !isnothing(list_ind)
                                position=read(file["$list_ind/pos"]).*monomer_size
                                bin_indices=fetch_all_spring_lengths(position[:,1:N], d_bins, incr)
                                for bin_ind in bin_indices #spring_lengths[bin_indices] .+=1 doesn't work for repeated indices
                                    spring_lengths[bin_ind]+=1
                                end
                                num_samples+=N
                            end
                        end
                    catch
                        println("couldnt get inds for $filename")
                    end
                end
                @assert num_samples<=typemax(Int)-1 "Overflow in number of samples"
                @assert maximum(spring_lengths)<typemax(Int)-1 "Overflow in spring length statistics!"
            end
            #normalize
            spring_lengths = spring_lengths ./num_samples

            #save all the statistics to a .h5 file
            out_file_name=replace(dir, data_parent_dir=>out_parent_dir)*"_10_kb_lengths.h5"
            h5open(out_file_name, "w") do file
                @write file spring_lengths
                @write file d_bins
            end
        end
    end
end

"""
    calc_hic(times, data_parent_dir, out_parent_dir; skip_done=true, ME_orient=false)

Calculates Hi-C maps for a folder of loop-extrusion simulation results

# Arguments
- `times::Array{Float64}`: Times at which to calculate statistics
- `data_parent_dir::String`: Path to the folder containing the simulation results
- `out_parent_dir::String`: Path to the folder where the statistics will be saved
- `skip_done::Bool=true`: If true, skips directories for which statistics have already been calculated
- `ME_orient::Bool=false`: If true, the ME orientation is used; strand 1 is ter-distal one
"""
function calc_hic(times, data_parent_dir,out_parent_dir; skip_done=true, ME_orient=false, hic_bins=405)

    subdirs=readdir(data_parent_dir, join=true)
    subdirs=subdirs[isdir.(subdirs)]
    mkpath(out_parent_dir)

    Threads.@threads for dir in subdirs 
        out_file_name=replace(dir, data_parent_dir=>out_parent_dir)*".h5"
        if contains(dir, "GPU") || (skip_done && isfile(out_file_name))
            println("Found previous results for $dir, skipping...")
        else
            println("Calculating contact map for $dir")
            num_samp=0

            N=parse_after(dir, "N_")

            if N<hic_bins
                N_hic=N
                incr=1
            else
                N_hic=hic_bins
                incr=Int(N/hic_bins)
            end

            #Initialize output arrays
            mean_forks=zeros(2, length(times))
            mean_hic=zeros(2*N_hic, 2*N_hic, length(times)) #the 4 is the dimension for all/old/new/inter

            #Loop over files, each containing a single time trajectory
            for filename in filter(contains(".h5"), readdir(dir,join=true))
                h5open(filename, "r") do file
                    list_inds=get_indices_at_times(file,times)

                    #Get all statistics at required times
                    for t in 1:length(times)
                        for i in list_inds[t]
                            try
                                pos, fork=fetch_pos_fork_at_ind(file, i, ME_orient=ME_orient)
                                pos=pos[:,1:incr:end]
                                fork[1]=floor(Int,fork[1]/incr)
                                fork[2]=ceil(Int, fork[2]/incr)

                                d_map=fetch_distance_matrix(pos,fork)
                                mean_forks[:,t].+=fork

                                c_array=fetch_contacts(d_map, fork, 4)
                                mean_hic[:,:,t].+=c_array
                                if t==1
                                    num_samp+=1
                                end
                            catch
                                println("Failed to read $i from $filename")
                            end
                        end
                    end
                end
            end

            #divide by number of samples to get average
            for array in [mean_forks, mean_hic]
                array ./= num_samp
            end
            #save all the statistics to a .h5 file
            if N_hic%5==0
                chunk_size=(Int(N_hic/5), Int(N_hic/5),1)
            elseif N%4==0
                chunk_size=(Int(N_hic/4), Int(N_hic/4),1)
            else
                chunk_size=(Int(N_hic/2), Int(N_hic/2),1)
            end

            out_file_name=replace(dir, data_parent_dir=>out_parent_dir)*"_hic_distances.h5"
			h5open(out_file_name, "w") do file
				@write file num_samp
				@write file mean_forks
                dset_hic = create_dataset(file, "mean_hic", datatype(mean_hic), size(mean_hic);  chunk=chunk_size, compress=5)
                write(dset_hic, mean_hic)
			end
        end
    end
end

"""
    fetch_msd_trajectories(initial_stage_files, monomer_size, tracked_loci=[1, 1012, 2025, 3037])

Fetches the mean squared displacements for a folder of loop-extrusion simulation results. Needs to be run for simulations with a higher time-resolution.

# Arguments
- `initial_stage_files::Array{String}`: Paths to the files containing the simulation results
- `monomer_size::Float64`: Monomer size in nm
- `tracked_loci::Array{Int}=[1, 1012, 2025, 3037]`: Loci to track
"""
function fetch_msd_trajectories(initial_stage_files, monomer_size, tracked_loci=[1, 1012, 2025, 3037])
    #Make arrays and set parameters
    r_squared=Dict{Int, Array{Float64}}() #for each sim, store an array x,y,t
    n_loci=length(tracked_loci)

    for (index, f) in enumerate(initial_stage_files)
        h5open(f) do file
            for i in keys(file)
                block=read_attribute(file[i], "block")+1
                sim = read_attribute(file[i], "simulation_run")
                pos=read(file[i], "pos")[2:3,tracked_loci].*monomer_size #2D diffusion
                if haskey(r_squared, sim)
                    r_squared[sim][:,:,block].=pos
                else
                    r_squared[sim]=zeros(2,n_loci, nt)
                    r_squared[sim][:, :, block].=pos
                end
            end
        end
    end
    #For each simulation, calculate r^2(t)
    for sim in keys(r_squared)
        ori_pos=r_squared[sim]
        ori_pos.-=ori_pos[:,:,1] #subtract first time point
        ori_pos.=ori_pos.^2
        r_squared[sim]=(ori_pos[1,:,:].+ori_pos[2,:,:]) #sum y and z
    end

    mean_rsq=zeros(n_loci, nt)
    error_rsq=zeros(n_loci, nt)
    for sim in keys(r_squared)
        mean_rsq.+=r_squared[sim]
        error_rsq.+=r_squared[sim].^2
    end
    mean_rsq./=length(keys(r_squared))
    error_rsq./=length(keys(r_squared))
    error_rsq.=sqrt.(error_rsq.-mean_rsq.^2)
    error_rsq./=sqrt(length(keys(r_squared)))

    return mean_rsq, error_rsq
end
