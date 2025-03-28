"""
For each folder in the ME_configs folder, analyze the data and save to .h5 files.

Also contains functions for inferring fork positions from Hi-C data, and for analyzing the decay of extension behind replication forks.
"""

include("statistics_functions.jl")
using Base.Threads

"""
    smooth_and_differentiate_periodic(data::Vector{Float64}, threshold::Int, σ=1)

Smooths and differentiates a given periodic vector.

# Arguments
    data: A vector of size N.
    threshold: An integer representing the smoothing threshold.
    σ: The standard deviation of the gaussian kernel.

# Returns
    smoothed_der: A vector of size N containing the smoothed derivative.
    smoothed_data: A vector of size N containing the smoothed data.
"""
function smooth_and_differentiate_periodic(data::Vector{Float64}, threshold::Int, σ=1)
    # Check if threshold is valid
    if threshold <= 0
        error("Threshold must be a positive integer.")
    end

    # Get the length of the input data
    N = length(data)

    # Use gaussian weights
    gaussian_weights=Kernel.gaussian([σ],[2*threshold+1])

    # Smooth the data with moving average for periodic data
    smoothed_data = [
        mean([data[mod(i - j - 1, N) + 1] for j in -threshold:threshold], weights(gaussian_weights))
        for i in 1:N
    ]

    # Calculate the numerical derivative
    derivative_curve = diff(smoothed_data) ./ diff(1:N)

    # Adjust the last derivative value for periodic data
    push!(derivative_curve, (smoothed_data[1] - smoothed_data[end]) / (N - (N - 1)))

    smoothed_der = [
        mean([derivative_curve[mod(i - j - 1, N) + 1] for j in -threshold:threshold], weights(gaussian_weights))
        for i in 1:N
    ]

    return vec(smoothed_der), vec(smoothed_data)
end

"""
    find_global_extrema_indices(smoothed_vector::Vector)

Takes a smoothened vector and returns the indices of its global extrema, i.e., the global maximum
and minimum indices.

# Arguments
- `smoothed_vector::Vector`: The smoothened input vector.

# Returns
- vector with indices corresponding to the inferred fork positions
"""
function find_fork_indices(smoothed_vector::Vector)
    N=length(smoothed_vector)
    ter=floor(Int,N/2)
    # Find the index of the global minimum
    (_, min_index) = findmin(smoothed_vector[1:ter])

    # Find the index of the global maximum
    (_, max_index) = findmax(smoothed_vector[ter:end])

    # Return a vector with indices
    return [min_index, max_index+ter]
end

"""
    calc_fast_ME_stats(times, lengths, file_dir,out_parent_dir; N=405, skip_done=true, ns=[1,2,3,4,5,6], monomer_size=0.088)

Calculate the statistics for each folder in the file_dir directory, and save to .h5 files in the out_parent_dir directory

# Arguments
- `times`: times to analyze
- `lengths`: lengths of the confinement
- `file_dir`: directory containing the input configurations
- `out_parent_dir`: directory to save the output .h5 files
- `N`: number of monomers
- `skip_done`: whether to skip files that already exist
- `ns`: list of n values to calculate local extensions
- `monomer_size`: size of the monomer
"""
function calc_fast_ME_stats(times, lengths, file_dir,out_parent_dir; N=405, skip_done=true, ns=[1,2,3,4,5,6], monomer_size=0.088)
    if !isdir(out_parent_dir)
        mkpath(out_parent_dir)
    end

    for subdir in subdirs(file_dir)
        out_file_name=replace(subdir, file_dir=>out_parent_dir)*"_fast.h5"

        if !isfile(out_file_name) || !skip_done
            mean_seg_frac=zeros(length(times))
            mean_zs=zeros(2*N,length(times))
            mean_zs_resolution_lim=zeros(2*N, length(times))
            counts_zs_resolution_limit=zeros(Int, 2*N, length(times))
            mean_local_extensions_different_n=zeros(2*N,length(times), length(ns))
            mean_sister_distances=zeros(N,length(times))
            mean_radial_pos=zeros(2*N,length(times))
            mean_forks=zeros(2,length(times))
            energies=zeros(N,N,length(times))
            max_h=time_to_height(times[end])/monomer_size #in monomer lengths
            hist_bins=-max_h/2:1:max_h/2
            z_histograms=zeros(length(hist_bins)-1, 3, length(times))
            num_samps=zeros(length(times))

            for (index,t) in enumerate(times)
                lin_file=h5open(subdir*"/lin_t$t.h5", "r")
                ring_file=h5open(subdir*"/ring_t$t.h5", "r")
                fork=read(ring_file, "fork")
                L=lengths[index]
                if haskey(ring_file,"energies")
                    energies[:,:,index].=read(ring_file,"energies")
                end
                num_samps[index]=length(keys(lin_file))
                for k in keys(lin_file)
                    pos=cat(read(ring_file,k), read(lin_file,k), dims=2)
                    mean_seg_frac[index]+=fetch_seg_fraction(pos, fork)./num_samps[index]
                    mean_zs[:,index].+=fetch_zs(pos,fork)./num_samps[index]
                    zs,resolved=fetch_zs_resolution_limited(pos, fork, N)
                    mean_zs_resolution_lim[resolved,index].+=zs[resolved]
                    counts_zs_resolution_limit[resolved,index].+=1
                    mean_radial_pos[:,index].+=fetch_radial_positions(pos, fork)./num_samps[index]
                    mean_sister_distances[:,index].+=fetch_sister_distances(pos,fork)./num_samps[index]
                    for (j,n) in enumerate(ns)
                        mean_local_extensions_different_n[:,index,j].+=fetch_local_extensions(pos,fork,n)./num_samps[index]
                    end

                    z1,z2,z3=fetch_z_histograms(pos, fork, hist_bins)
                    z_histograms[:,1,index].+=z1./num_samps[index]
                    z_histograms[:,2,index].+=z2./num_samps[index]
                    z_histograms[:,3,index].+=z3./num_samps[index]
                end
                mean_forks[:,index].=fork
                close(lin_file)
                close(ring_file)
            end
            mean_zs_resolution_lim./=counts_zs_resolution_limit
            hist_bins=collect(hist_bins)

            h5open(out_file_name, "w") do file
                @write file mean_zs
                @write file num_samps
                @write file mean_seg_frac
                @write file mean_sister_distances
                @write file mean_zs_resolution_lim
                @write file mean_forks
                @write file mean_radial_pos
                @write file energies
                @write file ns
                @write file mean_local_extensions_different_n
                @write file hist_bins
                @write file z_histograms
            end
        else
            println("File $out_file_name already exists, skipping...")
        end
    end
end

"""
    calc_hic(times=[0,10,30,45,60,75],file_dir="ME_configs/Compressed_configs/",out_parent_dir="Stats/"; N=405, skip_done=true)

Calculate the Hi-C statistics for each folder in the file_dir directory, and save to .h5 files in the out_parent_dir directory

# Arguments
- `times`: times to analyze
- `file_dir`: directory containing the input configurations
- `out_parent_dir`: directory to save the output .h5 files
- `N`: number of monomers
- `skip_done`: whether to skip files that already exist
"""
function calc_hic(times=[0,10,30,45,60,75],file_dir="ME_configs/Compressed_configs/",out_parent_dir="Stats/"; N=405, skip_done=true)
    if !isdir(out_parent_dir)
        mkpath(out_parent_dir)
    end

    for subdir in subdirs(file_dir)
        out_file_name=replace(subdir, file_dir=>out_parent_dir)*"_hic_distances.h5"

        if !isfile(out_file_name) || !skip_done
            mean_hic=zeros(2*N,2*N,length(times))
            mean_forks=zeros(2,length(times))
            for (index,t) in enumerate(times)
                lin_file=h5open(subdir*"/lin_t$t.h5", "r")
                ring_file=h5open(subdir*"/ring_t$t.h5", "r")
                fork=read(ring_file, "fork")

                num_samps=length(keys(lin_file))
                for k in keys(lin_file)
                    pos=cat(read(ring_file,k), read(lin_file,k), dims=2)
                    d_map=fetch_distance_matrix(pos, fork)
                    mean_hic[:,:,index].+=fetch_contacts(d_map,fork, 0.5)./num_samps
                end
                mean_forks[:,index].=fork
                close(lin_file)
                close(ring_file)
            end
            #save all the statistics to a .h5 file
            if N%5==0
                chunk_size=(Int(N/5), Int(N/5),1)
            elseif N%4==0
                chunk_size=(Int(N/4), Int(N/4),1)
            else
                chunk_size=(Int(N/2), Int(N/2),1)
            end
            h5open(out_file_name, "w") do file
                @write file mean_forks
                dset_hic = create_dataset(file, "mean_hic", datatype(mean_hic), size(mean_hic);  chunk=chunk_size, compress=5)
                write(dset_hic, mean_hic)
            end
        else
            println("File $out_file_name already exists, skipping...")
        end
    end
end

"""
    calc_spring_lengths(times, file_dir,out_parent_dir; skip_done=false, bin_width=0.5, monomer_size=88, N=405)

Calculate the spring lengths statistics for each folder in the file_dir directory, and save to .h5 files in the out_parent_dir directory

# Arguments
- `times`: times to analyze
- `file_dir`: directory containing the input configurations
- `out_parent_dir`: directory to save the output .h5 files
- `skip_done`: whether to skip files that already exist
- `bin_width`: width of the bins
- `monomer_size`: size of the monomer
- `N`: number of monomers
"""
function calc_spring_lengths(times, file_dir,out_parent_dir; skip_done=false, bin_width=0.5, monomer_size=88, N=405)
    if !isdir(out_parent_dir)
        mkpath(out_parent_dir)
    end

    for subdir in subdirs(file_dir)
        out_file_name=replace(subdir, file_dir=>out_parent_dir)*"_10_kb_lengths.h5"

        if !isfile(out_file_name) || !skip_done
            d_bins=collect(0.0:bin_width:500) #in the simulation units where monomer size is 1
            spring_lengths=zeros(Int, length(d_bins))
            num_samples=0
            for (index,t) in enumerate(times)
                ring_file=h5open(subdir*"/ring_t$t.h5", "r")
                for k in keys(ring_file)
                    pos=read(ring_file,k).*monomer_size
                    bin_indices=fetch_all_spring_lengths(pos, d_bins, 1)
                    for bin_ind in bin_indices[bin_indices.<=length(d_bins)] #spring_lengths[bin_indices] .+=1 doesn't work for repeated indices
                        spring_lengths[bin_ind]+=1
                    end
                    num_samples+=N
                end
                close(ring_file)
            end
            #normalize
            spring_lengths = spring_lengths ./num_samples

            #save all the statistics to a .h5 file
            h5open(out_file_name, "w") do file
                @write file spring_lengths
                @write file d_bins
            end
        else
            println("File $out_file_name already exists, skipping...")
        end
    end
end

"""
    load_extension_data(filename)

Helper function for loading the extension decay data from a file

# Arguments
- `filename`: name of the file
"""
function load_extension_data(filename)
    data = []
    open(filename, "r") do file
        for line in eachline(file)
            push!(data, split(line))
        end
    end
    return data
end

"""
    fetch_all_exponents(stat_directory)

Fetch all the exponents for extension decay from the files in the stat_directory

# Arguments
- `stat_directory`: directory containing the statistics files
"""
function fetch_all_exponents(stat_directory)
    filenames = filter(contains("strain_timescales_region_width_"), readdir(stat_directory, join=true))

    # Initialize an array to store all data
    all_data = []

    # Load data from each file
    for filename in filenames
        data = load_extension_data(filename)
        region_width = parse_after(filename, "region_width_")*10
        for row in data
            center_pos = round(Int,parse(Float64,row[1])*10) #in kb
            exponent = parse(Float64,row[2])
            if !isnan(exponent)
                push!(all_data, (center_pos, region_width, exponent))
            end
        end
    end

    # Sort data by region width and center position
    sort!(all_data, by = x -> (x[2], x[1]))
    return all_data
end

"""
    make_latex_table_timescales(stat_directory::String, plot_directory::String)

Generate a LaTeX table of the timescales for post-replication extension decay.

# Arguments
- `stat_directory::String`: The directory containing the statistical data files.
- `plot_directory::String`: The directory where the LaTeX table will be saved.

# Description
This function fetches all exponents from the given `stat_directory` and generates a LaTeX table
with the timescales for post-replication extension decay. The table is saved as `exponents_table.tex`
in the specified `plot_directory`.

# Output
The generated LaTeX table includes columns for center position (in kb), region width (in kb), and fit timescale (in minutes).
Rows with different region widths are highlighted with a gray background.
"""
# Function to make a LaTeX table of the timescales for post-replication extension decay
function make_latex_table_timescales(stat_directory,plot_directory)
    all_data=fetch_all_exponents(stat_directory)

    # Generate LaTeX table
    # File to output the LaTeX table
    output_file = joinpath(plot_directory, "exponents_table.tex")

    open(output_file, "w") do file
        println(file, "\\begin{table}[ht]")
        println(file, "\\centering")
        println(file, "\\begin{tabular}{|c|c|c|}")
        println(file, "\\hline")
        println(file, "\\textbf{Center position} \$i\$ [kb] & \\textbf{Region width \$n\$ [kb]} & \\textbf{Fit timescale [min]}\\\\")
        println(file, "\\hline")

        current_width = ""
        for row in all_data
            center_pos, region_width, exponent = row
            if region_width != current_width
                println(file, "\\rowcolor{gray!20}")
                current_width = region_width
            end
            println(file, "$center_pos & $region_width & $exponent \\\\")
            println(file, "\\hline")
        end

        println(file, "\\end{tabular}")
        println(file, "\\caption{\\textbf{Timescales of extension relaxation.} The inverse of the fit exponents for the relative extension after replication.}")
        println(file, "\\label{table:exponents}")
        println(file, "\\end{table}")
    end
end
