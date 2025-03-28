"""
Given that statistics have been calculated, this file can be used to make the plots
"""

using Plots, DelimitedFiles, HDF5, LaTeXStrings, Plots.PlotMeasures, StatsBase, CurveFit
include("helpers.jl")
include("statistics_functions.jl")

pythonplot(grid=false,label="",framestyle=:box,colorbar=true, linewidth=2,
    guidefontsize=15, tickfontsize=15,colorbar_tickfontsize=12,legend=:outerright, markersize=7,
    colorbar_titlefontsize=15, legendfontsize=15, palette=:tab10, topmargin=30px)

##################
# Helper functions

function xticks_Mb(N; replicated=false)
    increment=N/404*100
    if !replicated
        return (1:increment:N, string.([2, 3, 0, 1, 2]))
    else
        return (vcat(1:increment:3*increment+1, 5*increment+1:increment:8*increment+1), ["2", "3", "0", "1", "3'", "0'", "1'", "2'"])
    end
end

function yticks_Mb(N; replicated=false)
    increment=N/404*100

    if !replicated
        return (1+increment:increment:N, string.([3, 0, 1, 2]))
    else
        return (vcat(1+increment:increment:3*increment+1, 5*increment+1:increment:8*increment+1), ["3", "0", "1", "3'", "0'", "1'", "2'"])
    end
end

"""
Either pass a Hi-C maps for a single chromosome, or a replicated one (set replicated=true)

Plots the Hi-C map with monomer mid at the center of the plot,
if R is not NaN, draw dashed lines to indicate fork positions
"""
function plot_hic(M; fork=NaN, mid=1)
    N=size(M)[1]
    cticks=0:0.05:0.2

    M_shifted=shifted_map(M,mid)

    pl=heatmap(M_shifted,clims=(cticks[1],cticks[end]),color=cgrad(:BuPu), size=(500,450),
            xlims=(1,N),ylims=(1,N), xlabel="Genomic position [Mb]",
            colorbar_ticks=(cticks,cticks), aspect_ratio=1., 
            ylabel="Genomic position [Mb]", colorbartitle="Contact frequency",
            yticks=yticks_Mb(N), xticks=xticks_Mb(N))
    if !any(isnan.(fork))
        pl=vline!(pl,fork.+[N/2,-N/2], color=:red, linestyle=:dash)
        pl=hline!(pl,fork.+[N/2,-N/2], color=:red, linestyle=:dash)
    end
    return pl
end

"""
Given inter-chromosomal contacts, plots the trans-contact map with monomer mid at the center of the plot,
if R is not NaN, draw dashed lines to indicate fork positions
"""
function plot_inter_hic(M; fork=NaN, mid=1)
    N=size(M)[1]
    cticks=0:0.02:0.1

    M_shifted=shifted_map(M,mid)

    pl=heatmap(M_shifted,clims=(0,0.1),color=cgrad(:BuPu), size=(500,450),
            xlims=(1,N),ylims=(1,N), xlabel="Genomic position, Strand 1 [Mb]",
            colorbar_ticks=(cticks,cticks), aspect_ratio=1.,
            ylabel="Genomic position, Strand 2 [Mb]", colorbartitle="Contact frequency",
            xticks=xticks_Mb(N), yticks=yticks_Mb(N))
    if !any(isnan.(fork))
        pl=vline!(pl,fork.+[N/2,-N/2], color=:red, linestyle=:dash)
        pl=hline!(pl,fork.+[N/2,-N/2], color=:red, linestyle=:dash)
    end
    return pl
end

# Hi-C counts
function plot_counts(v, fork;mid=1)
    N=length(v)
    pl=plot(1:N,v,color=pu,xlims=[1,N], legend=false,grid=false, ylims=(0.6,1.4),
            ylabel="Hi-C count\nvs t=0 (rescaled)",xlabel="Genomic position (Mb)",
            xticks=xticks_Mb(N), xminorticks=4, tick_direction=:out, size=[500,300])
    plot!(pl, [fork[1],fork[1]],[0.6,1.4], linestyle=:dash, color=re)
    plot!(pl, [fork[2],fork[2]],[0.6,1.4], linestyle=:dash, color=re)

    return pl
end

#Plot mean squared displacements over time
function plot_mean_rsq(ts, mean_rsq)
    scatter(ts, mean_rsq[2,:], label=L"ter", color=:black, markerstrokewidth=0, alpha=0.5, xlabel="Simulation time [s]", ylabel="MSD [μm²]", size=(800, 400), scale=:log10)
    scatter!(ts, mean_rsq[1,:], label=L"ori", color=:red, markerstrokewidth=0, alpha=0.5)
    scatter!(ts, mean_rsq[3,:], label=L"90^{\circ}", color=3, markerstrokewidth=0, alpha=0.5)
    scatter!(ts, mean_rsq[4,:], label=L"-90^{\circ}", color=5, markerstrokewidth=0, alpha=0.5)

    weber_bottom=10^(-3)*ts.^0.4
    weber_mean=10^(-2)*ts.^0.4
    weber_top=5*10^(-2)*ts.^0.4

    plot!(ts,weber_mean, ribbon=(weber_mean.-weber_bottom, weber_top.-weber_mean), linestyle=:dash, color=:gray, alpha=0.5, fillcolor=:gray, fillalpha=0.1, label="Weber et al., PRL, 2010")
end

# Compute weighted quantiles for a discrete, sorted support.
function weighted_quantiles(bins::AbstractVector, weights::AbstractVector, qs::AbstractVector)
    total = sum(weights)
    cumw = cumsum(weights)
    quantiles = similar(qs)
    for (i,q) in pairs(qs)
        # find first index where cumulative weight exceeds q*total
        idx = findfirst(x -> x ≥ q * total, cumw)
        quantiles[i] = bins[idx]
    end
    return quantiles
end

#Plot loop size distributions over time; meadian and 25th and 75th percentiles
function plot_loop_size_distributions(data::AbstractMatrix, times::AbstractVector, kb_scale=10)
    N, nt = size(data)
    medians = zeros(nt)
    q25s = zeros(nt)
    q75s = zeros(nt)
    mins = zeros(nt)
    maxs = zeros(nt)
    bins = collect(1:N)
    for t in 1:nt
        weights = data[:, t]
        total = sum(weights)
        if total == 0
            medians[t] = NaN
            q25s[t]   = NaN
            q75s[t]   = NaN
            mins[t]   = NaN
            maxs[t]   = NaN
            continue
        end
        # weighted quantiles at 25%, 50%, 75%
        q = weighted_quantiles(bins, weights, [0.25, 0.5, 0.75])
        q25s[t], medians[t], q75s[t] = q[1], q[2], q[3]
        # determine min and max bins with nonzero weight
        min_bin = findfirst(>(0), weights)
        max_bin = findlast(>(0), weights)
        mins[t] = bins[min_bin]
        maxs[t] = bins[max_bin]
    end

    # Box error bars: difference from median to q25 and q75
    box_lower = medians .- q25s
    box_upper = q75s .- medians

    # Whisker error bars: difference from median to min and max
    whisker_lower = medians .- mins
    whisker_upper = maxs .- medians

    #Adjust the scale for the y-axis
    for v in [medians, box_lower, box_upper, whisker_lower, whisker_upper]
        v .*= kb_scale
    end

    plt = scatter(times, medians, yerr=(box_lower, box_upper), markersize=5, color=:black, size=(500,400), ylims=(0,800), xlabel="Time [min]", ylabel="Median loop size [kb]")
    return plt
end

#Plot the relative fraction of loops that are below a certain size over time
function plot_dens_small_loops(loop_size_weights, loop_size_bins, times, crit_loop_size=40; kwargs...)
    nt=length(times)
    max_bin=findfirst(loop_size_bins.>crit_loop_size)-1

    loops_below_n=map(t->sum(loop_size_weights[1:max_bin,t]), 1:nt) #fraction of loops below n
    loops_below_n./=loops_below_n[1]

    scatter(times, loops_below_n, xlabel="Time [min]", ylabel="Rel. dens. small loops", color=:black, ylims=(0.9,2.2); kwargs...)
    hline!([1], color=:grey, linestyle=:dash, alpha=0.5)
end

#Plot the loop size histograms over time
function plot_loop_size_histograms(loop_size_bins, loop_size_weights, times, N_kb=4050)
    N=length(loop_size_bins)
    nt=length(times)
    plot(xlabel="Loop size [kb]", ylabel="Probability", xlims=(0,N_kb/2), ylims=(0, 0.03))
    for i in 1:nt
        alpha=1-(i-1)/(nt+1)
        plot!(loop_size_bins, loop_size_weights[:,i], label="$(times[i]) min", color=i, alpha=alpha)
    end
end

"""
Plot densities across cell
"""
function plot_z_histograms(z_hist, bins, height, label="";kwargs...)
    plot([-height/2, -height/2, height/2, height/2], [0, 1, 1,0], color=:turquoise, fillalpha=0.2, ribbon=(1,0), xticks=-1.5:1.0:1.5;kwargs...)
    plot!(bins,z_hist[:,1], xlabel="Long axis position [μm]", ylabel="Mean frac. monomers", color=1, ylims=(0, 0.055), label=label;kwargs...)
    plot!(bins, z_hist[:,2], color=2;kwargs...)
    plot!(bins, z_hist[:,3], color=:black;kwargs...)
end
    
"""
Plot densities across cell
"""
function plot_z_histograms!(z_hist, bins, height, label="LE+ori-pull";kwargs...)
    plot!(bins,z_hist[:,1], color=1;kwargs...)
    plot!(bins, z_hist[:,2], color=2, label=label;kwargs...)
    plot!(bins, z_hist[:,3], color=:black;kwargs...)
end

#############
# Mean z-position plots
"""
Given average long axis positions from simulations with and without smcs,
plot them on the same figure, with an indication of the cell size and fork positions
"""
function plot_av_z(av1, fork, L)
    N=Int(length(av1)/2)
    ter=ceil(Int,N/2)

    unrepl1=fork[1]:ter-1
    unrepl2=ter:min(N,fork[2]+1)

    repl=vcat(min(N,fork[2]+1):N, 1:fork[1])
    repl_inds= repl.-(ter)
    repl_inds[repl_inds.<=0].+=N+1

    old=av1[1:N]
    new=av1[N+1:end]

    hline([L/2], ribbon=(L, 0), color=:turquoise, fillalpha=0.2, ylims=(-1.55,1.55))
    hline!([-L/2], color=:turquoise, xticks=xticks_Mb(N))
    hline!([0], color=:grey, linealpha=0.3)
    plot!(repl_inds,old[repl], ylabel="Long axis position [μm]",xlabel="Genomic position [Mb]", color=1)
    plot!(unrepl2.-(ter),old[unrepl2], color=:black, label="MaxEnt simulations")
    plot!(unrepl1.-(ter-N+1),old[unrepl1], color=:black)
    plot!(repl_inds,new[repl], color=2, size=[700,400])
end

"""
Given average long axis positions from simulations with and without smcs,
plot them on the same figure, with an indication of the cell size and fork positions
"""
function plot_av_zs_all(avs, forks, times, L=3.1)
    N=Int(size(avs)[1]/2)
    ter=floor(Int,N/2)
    nt=size(avs)[2]

    hline([0], color=:grey, linealpha=0.3, ylims=(-L/2,L/2), xticks=xticks_Mb(N))

    for i in 1:nt
        α=1-(i-1)/(nt)
        f1=max(1,forks[1,i]) #first fork; always between 0, N/2
        f2=min(N, forks[2,i]) #second fork: always between N/2, N-1

        unrepl1=f1:ter
        unrepl2=ter+1:f2

        repl=vcat(f2+1:N, 1:f1-1)
        repl_inds=repl.-(ter-1)
        repl_inds[repl_inds.<=0].+=N
        
        old=avs[1:N,i]
        new=avs[N+1:end,i]

        if i==1
            unrepl2=vcat(ter:N,1:ter)
            plot!(old[unrepl2], color=:black, label="$(times[i]) min")
        else 
            plot!(repl_inds,old[repl], ylabel="Long axis position [μm]",xlabel="Genomic position [Mb]", color=1, alpha=α)
            plot!(repl_inds,new[repl], color=2, size=[550,400], alpha=α)
            plot!(unrepl2.-(ter),old[unrepl2], color=:black, alpha=α, label="$(times[i]) min")
            plot!(unrepl1.-(ter-N),old[unrepl1], color=:black, alpha=α)
        end

    end
    return current()
end

################
# Local extensions

#Just one chromosome
function plot_abs_local_ext(data::Array,scale::Number, fork=(NaN,NaN))
    N=length(data)
    ter=N/2
    plot(shifted_vector(data), size=[500,300], ylabel="Local extension [μm]",
        ylims=[0.,0.4], color=:black, xticks=xticks_Mb(N), xlabel="Genomic position [Mb]")
    if !any(isnan,fork)
        vline!(fork+[ter,-ter], linestyle=:dash, color=4,linewidth=1,alpha=0.5)
    end
    return current()
end

#Partially replicated
function plot_abs_local_ext(old::Array,new::Array,scale::Number,fork)
    N=length(old)

    ter=ceil(Int,N/2)

    unrepl1=fork[1]-scale:ter-1
    unrepl2=ter:min(N,fork[2]+scale)

    repl=vcat(min(N,fork[2]+1+scale):N, 1:max(1,fork[1]-1-scale))
    repl_inds= repl.-(ter)
    repl_inds[repl_inds.<=0].+=N+1

    plot(unrepl1.-(ter-N+1), old[unrepl1], size=[5/3*270,270], ylabel="Local extension [μm]",
            ylims=[0.,0.4],color=:black, xticks=xticks_Mb(N), xlabel="Genomic position [Mb]")
    plot!(unrepl2.-ter,old[unrepl2], color=:black)
    plot!(repl_inds,old[repl],color=1)
    plot!(repl_inds,new[repl], color=2)
    vline!(fork+[ter,-ter], linestyle=:dash, color=4, linewidth=1, alpha=0.5)
end

#Change/original; partially replicated
function plot_rel_local_ext(old::Array,new::Array,data0::Array,scale::Number,fork;maxy=0.2, miny=-0.2)
    N=length(old)
    exts=(old.-data0)./data0
    exts_new=(new.-data0)./data0

    ter=ceil(Int,N/2)

    unrepl1=max(fork[1]-scale,1):ter-1
    unrepl2=ter:min(N,fork[2]+scale)

    repl=vcat(min(N,fork[2]+1+scale):N, 1:max(1,fork[1]-1-scale))
    repl_inds= repl.-(ter)
    repl_inds[repl_inds.<=0].+=N+1

    plot(unrepl1.-(ter-N+1), exts[unrepl1], size=[5/3*270,270], ylabel="Relative extension",
            ylims=(miny,maxy),color=:black, xticks=xticks_Mb(N), xlabel="Genomic position [Mb]")
    plot!(unrepl2.-ter,exts[unrepl2], color=:black)
    plot!(repl_inds,exts[repl],color=1)
    plot!(repl_inds,exts_new[repl], color=2)

    hline!([0], color=8, alpha=0.6)
    vline!(fork.+[ter,-ter], linestyle=:dash, color=4,linewidth=1,alpha=0.5)
end

#Change/original in a kymograph. Only plotted for ter-distal strand and unreplicated strand
function plot_local_ext_kymograph(kymograph_array,times, forks)
    N=size(kymograph_array)[1]
    nt=length(times)

    heatmap(circshift(kymograph_array, (floor(Int,N/2),0)), xticks=(1:nt-1, times[2:end]), xlabel="Time [min]", 
        ylabel="Genomic location [Mb]", color=cgrad(:bwr), clims=(-0.2,0.2),colorbartitle="Relative extension", 
        yticks=yticks_Mb(N))
    plot!(1:nt-1, forks[1,2:end].+N/2, color=:black, linestyle=:dot)
    plot!(1:nt-1, forks[2,2:end].-N/2, color=:black, linestyle=:dot)
end

#Calculate the increase in density over time, based on replication stage and volume, plus the mean local extension
function plot_ext_vs_density(forks, times, ext_n)
    N=Int(size(ext_n)[1]/2)
    densities=map(i->rel_monomer_density(forks[:,i],N,times[i]), 1:length(times))

    rel_compactions=copy(ext_n[:,:,2])
    for i in 2:length(times)
        rel_compactions[:,i]./=rel_compactions[:,1]
    end
    rel_compactions[:,1].=1

    mean_compaction=map(i->mean(rel_compactions[1:N,i]), 1:length(times))
    scatter(densities.^(-1/3), mean_compaction, xlabel=L"(\rho/\rho_0)^{-1/3}", 
            ylabel=L"\left\langle d_n(i,t)/d_n(i, t=0) \right\rangle", 
            size=[500,400], xticks=:auto, ylims=(0.9, 1.025), color=:black)
end

function ret_ext_vs_density(forks, times, ext_n)
    N=Int(size(ext_n)[1]/2)
    densities=map(i->rel_monomer_density(forks[:,i],N,times[i]), 1:length(times))

    rel_compactions=copy(ext_n[:,:,2])
    for i in 2:length(times)
        rel_compactions[:,i]./=rel_compactions[:,1]
    end
    rel_compactions[:,1].=1

    mean_compaction=map(i->mean(rel_compactions[1:N,i]), 1:length(times))
    return densities, mean_compaction
end

#Loop-extruder density plot
function plot_chipseq(chipseq, fork=NaN; kwargs...)
    N=length(chipseq)
    pl=plot(shifted_vector(chipseq), xlabel="Genomic position [Mb]", ylabel="Mean number loop-extruders", 
            color=:black, ylims=(0,1.1), size=[500,350], xticks=xticks_Mb(N); kwargs...)
    if ! any(isnan, fork)
            vline!(pl,fork.+[N/2,-N/2], linestyle=:dash, color=4, linewidth=1, alpha=0.5)
    end
    return pl
end

#Loop-extruder density compared to t=0
function plot_relative_chipseq(chipseq, chipseq_0, fork=NaN)
    rel_chipseq=copy(chipseq)
    N=length(chipseq)
    replicated=vcat(1:fork[1], fork[2]:N)
    rel_chipseq[replicated]./=2
    rel_chipseq./=chipseq_0
    pl=plot(shifted_vector(rel_chipseq), xlabel="Genomic position [Mb]", ylabel="Mean occupancy relative to t=0", 
            color=:black, size=[500,350], xticks=xticks_Mb(N), ylims=(0,5))
    hline!(pl,[1], color=8, alpha=0.6)
    if ! any(isnan, fork)
            vline!(pl,fork.+[N/2,-N/2], linestyle=:dash, color=4, linewidth=1, alpha=0.5)
    end
    return pl
end


"""
    plot_relative_counts(rel_col_sums, smoothed, found_forks, t, plot_folder)

Plots the relative Hi-C count sum, the smoothed vector, and the inferred forks, saves a png file.

# Arguments
- `rel_col_sums`: The relative Hi-C count sum.
- `smoothed`: The smoothed vector.
- `found_forks`: The inferred forks.
- `t`: The time point.
- `plot_folder`: The folder where the png file is saved.

# Returns
- nothing
"""
function plot_relative_counts(rel_col_sums, smoothed, found_forks, t, plot_folder)
    plot(shifted_vector(rel_col_sums), title="t=$t'", ylabel="Relative Hi-C count sum",
        ylims=(0.6,1.4), color=:black, xlabel="Genomic position [Mb]")
    plot!(shifted_vector(smoothed), color=:lightgray, size=(600,300), linewidth=2.3)
    vline!(found_forks.+[N/2,-N/2], linestyle=:dash, color=:red)
    savefig(plot_folder*"smoothed_$t.pdf")
end

"""
    plot_ders(ders, found_forks, t, plot_folder)

Plots the derivative, and the inferred forks, saves a png file.

# Arguments
- `ders`: The derivative.
- `found_forks`: The inferred forks.
- `t`: The time point.
- `plot_folder`: The folder where the png file is saved.

# Returns
- nothing
"""
function plot_ders(ders, found_forks, t, plot_folder)
    plot(shifted_vector(ders), title="t=$t'", ylims=(-0.006,0.006), ylabel="Derivative", xlabel="Genomic position [Mb]", color=:black)
    vline!(found_forks.+[N/2,-N/2], linestyle=:dash, color=:red, size=(600,300))
    savefig(plot_folder*"ders_$t.pdf")
end

"""
    plot_inferred_forks(times, inferred_forks_all, fit_right, fit_left, plot_folder)

Plots the inferred forks and the linear fits, saves a png file.

# Arguments
- `times`: The time points used for the fit.
- `inferred_forks_all`: The inferred forks.
- `fit_right`: The linear fit for the right arm.
- `fit_left`: The linear fit for the left arm.
- `plot_folder`: The folder where the png file is saved.

# Returns
- nothing
"""
function plot_inferred_forks(times, inferred_forks_all, fit_right, fit_left, plot_folder)
    scatter(times, inferred_forks_all[1,:], label="Right arm", ylims=(0,2),xticks=:auto, xlabel="Cell cycle progression [min]", ylabel="Fork progression [Mb]", color=1)
    plot!(times, fit_right[1].+fit_right[2].*times, label="", color=1)
    scatter!(times,N./100 .-inferred_forks_all[2,:], label="Left arm", color=2)
    plot!(times, fit_left[1].+fit_left[2].*times, label="", color=2, size=(600,300))
    savefig(plot_folder*"forks_plus_fit.pdf")
end

"""
        validation_plot(ts, oclock, dir_exp, dir_le, t_repl, cond_av=true)

Plots the average long axis position over time, for a given oclock, and compares to experimental data

# Arguments
- `ts::Array{Float64}`: Time points for the simulation
- `oclock::Int`: The oclock value for the locus
- `dir_exp::String`: The directory containing the experimental data
- `dir_le::String`: The directory containing the simulation data
- `t_repl::Float64`: The time point at which the replication fork passes locus
- `cond_av::Bool`: Whether to plot the conditional average (true) or the raw data (false)

# Returns
- `Plots.Plot`: The plot object
"""
function validation_plot(ts, oclock, dir_exp, dir_le, t_repl, cond_av=true)
    exp_m=readdlm(dir_exp*"$(oclock)_means.txt")[1:6,:]
    exp_std=readdlm(dir_exp*"$(oclock)_stds.txt")[1:6,:]

    if cond_av
        le_m=readdlm(dir_le*"/means_$oclock")
        le_std=readdlm(dir_le*"/stds_$oclock")
    else
        le_m=readdlm(dir_le*"/raw_means_$oclock")
        le_std=readdlm(dir_le*"/raw_stds_$oclock")
    end

    unrepl=exp_m[:,1].==exp_m[:,2]
    repl=exp_m[:,1].!=exp_m[:,2]

    plot(ts[unrepl], exp_m[unrepl,1],ribbon=exp_std[unrepl,1],color=:black,xticks=:auto, ylims=(0,1), size=[500,400])
    plot!(ts[repl], exp_m[repl,:],ribbon=exp_std[repl,:],color=[3 5],xlabel="Time [min]", ylabel="Scaled long axis position")
    scatter!(ts[unrepl],le_m[unrepl,1],yerror=le_std[unrepl,1],color=:black, markersize=10, markershape=:x,xticks=:auto)
    scatter!(ts[repl],le_m[repl,:],yerror=le_std[repl,:],color=[3 5], markersize=10, markershape=:x,
        markerstrokecolor=[3 5], yticks=0.2:0.2:1.0)
    
    if t_repl<ts[end]
        plot!([ts[unrepl][end],t_repl],[exp_m[unrepl,1][end],exp_m[unrepl,1][end]], color=:black,linestyle=:dash)
        plot!([t_repl, ts[repl][1]],[exp_m[unrepl,1][end],exp_m[repl,2][1]], color=5,linestyle=:dash)
        plot!([t_repl,ts[repl][1]],[exp_m[unrepl,1][end],exp_m[repl,1][1]], color=3,linestyle=:dash)
    end
    return current()
end

###################
# Main functions for plotting

function make_fast_plots(times, data_dir,plot_dir; skip_done=true, fig_type=".pdf")
    if fig_type[1]!='.'
        fig_type="."*fig_type
    end

    #make directories for plots that we don't want to fill up the simulation specific directories
    mkpath(plot_dir*"Local_ext")
    mkpath(plot_dir*"Rel_local_ext")

    nt=length(times)
    files=filter(contains("fast.h5"), readdir(data_dir, join=true))
    sample_counts=[]
    for filename in files
        if contains(filename, "_b_")
            monomer_size=parse_after(filename, "_b_")/1000
        else
            monomer_size=88/1000
        end

        #directories for single simulations
        out_dir=replace(filename, data_dir=>plot_dir)
        out_dir=replace(out_dir, ".h5"=>"/")

        if skip_done && isdir(out_dir)
            println("Found directory $out_dir. Skipping.")
        else
            mkpath(out_dir)
            println("Making plots for $filename...")

            #file handle for naming plots
            out_ext=replace(filename, data_dir=>"")
            out_ext=replace(out_ext, "/"=>"_")
            out_ext=replace(out_ext, ".h5"=>"")

            h5open(filename, "r") do file
                #All mean z positions over time in one plot
                zs=read(file, "mean_zs").*monomer_size
                N=Int(size(zs)[1]/2)
                forks=floor.(Int,read(file, "mean_forks"))
                push!(sample_counts, read(file, "num_samps"))
                forks[1,forks[1,:].>N/2].=floor(Int,N/2)
                forks[2,forks[2,:].<N/2].=floor(Int,N/2)

                pl=plot_av_zs_all(zs, forks, times)
                savefig(pl, out_dir*"mean_z_positions"*fig_type)

                ext_n=read(file,"mean_local_extensions_different_n")*monomer_size
                ns=read(file, "ns")
                
                z_histograms=read(file, "z_histograms")
                z_bins=read(file, "hist_bins")[1:end-1].*monomer_size

                #if the simulations included loop-extruders, we look at the loop sizes
                if haskey(file, "mean_chipseq")
                    chipseq=read(file, "mean_chipseq")
                    loop_size_bins=read(file, "loop_size_bins")
                    loop_size_weights=read(file, "loop_size_distr")

                    #Whisker plots of loop sizes
                    plot_loop_size_distributions(loop_size_weights, times)
                    savefig(out_dir*"loop_size_distributions"*fig_type)

                    #Check how loop size, extension, and density are related
                    plot_dens_small_loops(loop_size_weights, loop_size_bins, times)
                    savefig(out_dir*"loops_below_n_density_over_time"*fig_type)

                    plot_loop_size_histograms(loop_size_bins, loop_size_weights, times)
                    savefig(out_dir*"loop_size_hist"*fig_type)
                end

                #Local extension kymographs for different values of n
                for (n_index,n) in enumerate(ns)
                    ext=ext_n[:,:,n_index]
                    kymograph_array=extension_kymograph(ext)[:,2:end]
                    plot_local_ext_kymograph(kymograph_array,times, forks)
                    savefig(out_dir*"extension_kymograph_n_$(n)"*fig_type)
                end

                #extension vs density
                plot_ext_vs_density(forks, times, ext_n)
                savefig(out_dir*"ext_vs_density"*fig_type)

                #For each time point, plot local extension (changes), z-histogram, and chipseq
                for (index,t) in enumerate(times)
                    file_prefix=out_dir*"t_$(t)_"
                    file_prefix_shared_dir=out_ext*"_t_$(t)_"

                    #Local extension; line plots rather than kymographs
                   for (n_index,n) in enumerate(ns)
                       ext=ext_n[:,:,n_index]
                       if index==1
                           pl=plot_abs_local_ext(ext[1:N,1],2)
                       else
                           pl=plot_abs_local_ext(ext[1:N,index], ext[N+1:end,index], 2, forks[:,index])
                       end
                       if n==2
                           savefig(pl, file_prefix*"local_ext_n_$(n)"*fig_type)
                       end
                       savefig(pl, plot_dir*"Local_ext/"*file_prefix_shared_dir*"_n_$(n)"*fig_type)
                       if index>1
                           max_y= index==2 ? 1 : 0.5
                           #Relative local extension
                           pl=plot_rel_local_ext(ext[1:N,index], ext[N+1:end,index], ext[1:N,1], 2, forks[:,index], maxy=max_y)
                           savefig(pl, plot_dir*"Rel_local_ext/"*file_prefix_shared_dir*"_n_$(n)"*fig_type)
                       end

                       if n==2
                           max_y= index==2 ? 0.6 : 0.2
                           #Relative local extension
                           pl=plot_rel_local_ext(ext[1:N,index], ext[N+1:end,index], ext[1:N,1], 2, forks[:,index], maxy=max_y)
                           savefig(pl, file_prefix*"relative_local_ext_n_$(n)"*fig_type)
                       end
                   end

                    #Z histograms
                    if index>1
                        pl=plot_z_histograms(z_histograms[:,:,index], z_bins, time_to_height(t)+monomer_size, size=[380,320])
                        savefig(pl, file_prefix*"z_histograms"*fig_type)
                    end

                    #ChIP-seq
                    if haskey(file, "mean_chipseq")
                        pl=plot_chipseq(chipseq[:,index], forks[:,index])
                        savefig(pl,file_prefix*"chipseq"*fig_type)

                        if index>1
                            pl=plot_relative_chipseq(chipseq[:,index], chipseq[:,1], forks[:,index])
                            savefig(pl,file_prefix*"relative_chipseq"*fig_type)
                        end
                    end
                end
            end
        end
    end
    writedlm(plot_dir*"check_num_samps.txt", hcat(files,sample_counts)) 
end

function make_hic_plots(times, data_dir,plot_dir; skip_done=true, fig_type=".png")
    if fig_type[1]!='.'
        fig_type="."*fig_type
    end

    for subf in ["Hi_C/all", "Hi_C/old", "Hi_C/new", "Hi_C/inter"]
        if !isdir(plot_dir*subf)
            mkpath(plot_dir*subf)
        end
    end

    subdirs=readdir(data_dir, join=true)
    subdirs=subdirs[isdir.(subdirs)]

    Ls=time_to_height.(times)

    files=filter(contains("hic_distances.h5"), readdir(data_dir, join=true))
    for filename in files
        out_ext=replace(filename, data_dir=>"")
        out_ext=replace(out_ext, "/"=>"_")
        out_ext=replace(out_ext, ".h5"=>"")
        file_prefix_shared_dir=out_ext*"_t_$(times[end])_"

        if skip_done && isfile(plot_dir*"Hi_C/all/"*file_prefix_shared_dir*fig_type)
            println("Found plots for $filename. Skipping.")
        else
            println("Making plots for $filename...")

            h5open(filename, "r") do file
                forks=floor.(Int,read(file, "mean_forks"))
                #individual hic, local extension and chipseq maps for each time point
                contacts=read(file, "mean_hic")

                for (index,t) in enumerate(times)
                    file_prefix_shared_dir=out_ext*"_t_$(t)_"
                    all_hic=contacts[:,:,index]
                    N_hic=Int(size(all_hic)[1]/2)

                    #Hi-C maps
                    hic_old=all_hic[1:N_hic,1:N_hic]
                    hic_new=all_hic[N_hic+1:2*N_hic, N_hic+1:2*N_hic]
                    hic_inter=all_hic[N_hic+1:2*N_hic, 1:N_hic]
                    hic_all=hic_old.+hic_new.+hic_inter.+transpose(hic_inter)

                    pl=plot_hic(hic_all, fork=forks[:,index])
                    savefig(pl, plot_dir*"Hi_C/all/"*file_prefix_shared_dir*fig_type)

                    pl=plot_hic(hic_old, fork=forks[:,index])
                    savefig(pl, plot_dir*"Hi_C/old/"*file_prefix_shared_dir*fig_type)

                    pl=plot_hic(hic_new, fork=forks[:,index])
                    savefig(pl, plot_dir*"Hi_C/new/"*file_prefix_shared_dir*fig_type)

                    if t>0
                        pl=plot_inter_hic(hic_inter, fork=forks[:,index])
                        savefig(pl, plot_dir*"Hi_C/inter/"*file_prefix_shared_dir*fig_type)
                    end
                end
            end
        end
    end
end

function make_spring_length_histograms(stat_parent_dir="Stats_spring_lengths/", plot_dir="Plots/Spring_lengths/"; fig_type=".pdf")
    stat_files=readdir(stat_parent_dir, join=true)
    stat_files=filter(contains("10_kb_lengths.h5"), stat_files)
    for f in stat_files
        if contains(f, "_b_")
            monomer_size=parse_after(f, "_b_")
        else
            monomer_size=88
        end

        ds=h5read(f, "d_bins")
        bin_vals=h5read(f, "spring_lengths")
        bin_width=ds[2]-ds[1]
        mean_length=mean(ds, weights(bin_vals))
        std_l=std(ds, weights(bin_vals))

        max_x=500
        min_x=0

        bar(ds, bin_vals, width=bin_width, xticks=:auto, xlabel="10 kb separation [nm]", ylabel="Frequency", size=[500,400], xlims=(min_x, max_x),ylims=(0,0.6),linealpha=0, color=8, legend=:topright)

        out_file_name=replace(f, stat_parent_dir=>plot_dir*"/10_kb_distances/")
        out_file_name=replace(out_file_name, ".h5"=>fig_type)
        mkpath(dirname(out_file_name))
        savefig(out_file_name)
    end
end


"""
    make_validation_plots(times,oclocks, locus_replication_times, validation_data_dir, experimental_dir, plot_parent_dir)

Make plots comparing the experimental data to the MaxEnt predictions

# Arguments
- `times`: the times to calculate the mean separations
- `oclocks`: the locus positions
- `locus_replication_times`: the times of replication in minutes
- `validation_data_dir`: the directory containing the validation data
- `experimental_dir`: the directory containing the experimental data
- `plot_parent_dir`: the directory to output the plots

# Output
- nothing, but writes the plots to files
"""
function make_validation_plots(times,oclocks, locus_replication_times, validation_data_dir, experimental_dir, plot_parent_dir, fig_type=".pdf")
    if fig_type[1]!='.'
        fig_type="."*fig_type
    end

    for dir in subdirs(validation_data_dir)
        if contains(dir, "Experiment")
            continue
        end
        plotdir=replace(dir, validation_data_dir=>plot_parent_dir)
        mkpath(plotdir)

        for (ind,oclock) in enumerate(oclocks)
            pl=validation_plot(times,oclock,experimental_dir,dir,locus_replication_times[ind])
            savefig(plotdir*"/$(oclock)_oclock_cond_av"*fig_type)

            pl=validation_plot(times,oclock,experimental_dir,dir,locus_replication_times[ind], false)
            savefig(plotdir*"/$(oclock)_oclock_raw_av"*fig_type)
        end
    end           
end

"""
    plot_energies(energies, mid=1, clim=2)

Plot the effective energies of a MaxEnt model.

# Arguments
- `energies`: The effective energies.
- `mid`: The middle of the energy map.
- `clim`: The color limit.

# Returns
- `Plots.Plot`: The plot object
"""
function plot_energies(energies, mid=1, clim=2)
    N=size(energies)[1]
    ter=floor(Int,N/2)
    heatmap(shifted_map(energies,mid+ter), color=cgrad(:bwr), clims=(-clim,clim), 
            aspect_ratio=1, size=(450,400), xlabel="Genomic position [Mb]", ylabel="Genomic position [Mb]",
            colorbartitle="Effective interaction energy", xticks=xticks_Mb(N), yticks=yticks_Mb(N))
end

"""
    visualize_energies(folder::String, out_dir::String, best_try::Int)

Plot the energy map for the simulation with the highest Pearson correlation coefficient and the energies for the constrained sites over simulation time.
Also plot the effective energy for ori separation

# Arguments
- `folder::String`: The path to the simulation output directory.
- `out_dir::String`: The output directory to save the generated plot.

# Notes
The function saves the plots as PNGs in the provided output directory.
"""
function visualize_energies(in_dir, out_dir)
    ts=[]
    ori_es=[]
    for folder in subdirs(in_dir) #loop over simulation output files
        time_stamp = parse(Int, folder[end-1:end])
        push!(ts, time_stamp)
        #First plot the energy map for the best try
        en_file=filter(x->occursin("/energies_", x), readdir(folder, join=true))[1]
        energies =readdlm(en_file)

        plot_energies(energies,306)
        savefig(out_dir*"energies_t_$time_stamp.png")

        ori_en_file=filter(x->occursin("/sep_energ_", x), readdir(folder, join=true))
        if length(ori_en_file)==0
            push!(ori_es, 0)
        else
            push!(ori_es, readdlm(ori_en_file[1])[1,1])
        end

    end
    plot(ts, ori_es, xlabel="Time [min]", ylabel="Origin separation energy", color=:black, size=(480,400))
    savefig(out_dir*"ori_energies.pdf")
end

"""
    plot_strain_forks(sampled_ts, extension_front, fit_fork1, fit_fork2)

Plot the strain in front of the forks, as well as exponential fits.

# Arguments
- `sampled_ts`: The time points.
- `extension_front`: The extension in front of the forks.
- `fit_fork1`: The fit for the first fork.
- `fit_fork2`: The fit for the second fork.
"""
function plot_strain_forks(sampled_ts, extension_front, fit_fork1, fit_fork2)
    plot(t->fit_fork1.k+fit_fork1.p[1]*exp(fit_fork1.λ[1]*t), label="Fit 1", color=:black,  ylabel="Rel. extension ahead of fork", xlabel="Time [min]", size=(700,400), xlims=(8,78), ylims=(-0.08, 0.05), yscale=:identity)
    scatter!(sampled_ts[2:end],extension_front[2:end,1], label="Fork 1",color=:black)
    plot!(t->fit_fork2.k+fit_fork2.p[1]*exp(fit_fork2.λ[1]*t), label="Fit 2", color=:black, linestyle=:dash)
    scatter!(sampled_ts[2:end],extension_front[2:end,2], label="Fork 2", color=:black, marker=:diamond,alpha=0.5)
    #inset with linlog scale
    plot!(t->fit_fork1.p[1]*exp(fit_fork1.λ[1]*t), color=:black, ylabel="", xlabel="",yscale=:log10, ylims=(1e-4,0.5), xlims=(1,76), xticks=15:15:75, inset = (1, bbox(0.55,0.05,0.4,0.4)), subplot=2)
    plot!(t->fit_fork2.p[1]*exp(fit_fork2.λ[1]*t), color=:black, linestyle=:dash, subplot=2)
    scatter!(sampled_ts[2:end],extension_front[2:end,1].-fit_fork1.k, ylabel="", xlabel="", color=:black, subplot = 2)
    scatter!(sampled_ts[2:end],extension_front[2:end,2].-fit_fork2.k, alpha=0.5,color=:black, subplot = 2, marker=:diamond)
end

"""
    plot_strain_forks(sampled_ts, extension_front)

Plot the strain in front of the forks. No exponential fits, but a horizontal line for mean extension.

# Arguments
- `sampled_ts`: The time points.
- `extension_front`: The extension in front of the forks.
"""
function plot_strain_forks(sampled_ts, extension_front)
    scatter(sampled_ts[2:end],extension_front[2:end,1], label="Fork 1", ylabel="Rel. extension ahead of fork", xlabel="Time [min]", color=:black, size=(750,400), ylims=(-0.078, 0.05))
    scatter!(sampled_ts[2:end],extension_front[2:end,2], label="Fork 2", color=:black, marker=:diamond,alpha=0.5)
    hline!([0], linewidth=0.5, color=:black)
    hline!([mean(extension_front[3:end,:])], linestyle=:dash, color=:grey, label="Mean t>10 min", ribbon=std(extension_front[3:end,:]))
end

"""
    plot_local_extensions_ahead_of_forks(forks, extensions, sampled_ts, plot_directory, ranges_ahead)

Plot the local extensions in front of the forks.

# Arguments
- `forks`: The fork positions.
- `extensions`: The extensions.
- `sampled_ts`: The time points.
- `plot_directory`: The directory to save the plots.
- `ranges_ahead`: The ranges ahead of the forks.
"""
function plot_local_extensions_ahead_of_forks(forks, extensions, sampled_ts, plot_directory, ranges_ahead)
    nt=length(sampled_ts)
    #Fetch the extensions in front of the forks
    for average_over in ranges_ahead #bins
        fork_ranges_1=map(i->forks[1,i]:forks[1,i]+average_over, 1:nt) #ranges ahead of forks
        fork_ranges_2=map(i->forks[2,i]-average_over:forks[2,i], 1:nt)
        extension_front_1=map(i->mean(extensions[fork_ranges_1[i],i]), 1:nt) #mean extension in front of forks
        extension_front_2=map(i->mean(extensions[fork_ranges_2[i],i]), 1:nt)
        extension_front=hcat(extension_front_1, extension_front_2)

        #Fit exponential to the extension in front of the forks; see if it builds up
        fit_fork1=expsum_fit(sampled_ts[2:end], extension_front_1[2:end], 1, withconst=true)
        fit_fork2=expsum_fit(sampled_ts[2:end], extension_front_2[2:end], 1, withconst=true)

        plot_strain_forks(sampled_ts, extension_front, fit_fork1, fit_fork2)
        savefig(plot_directory*"extension_front_with_fit_av_over_$(average_over).pdf")

        plot_strain_forks(sampled_ts, extension_front)
        savefig(plot_directory*"extension_front_const_av_over_$(average_over).pdf")
    end
end

"""
    plot_strain_over_time(sampled_ts, mean_extensions, fit, replicated_at, monomer_ind)

Plot the strain over time for different genomic regions.

# Arguments
- `sampled_ts`: The time points.
- `mean_extensions`: The mean extensions.
- `fit`: The fit.
- `replicated_at`: The time of replication.
- `monomer_ind`: The monomer index.
"""
function plot_strain_over_time(sampled_ts, mean_extensions, fit, replicated_at, monomer_ind)
    #First make the plot on a linear scale
    plot(t->fit.k[1]+fit.p[1]*exp(fit.λ[1]*t), label="Exponential fit", color=:gray, ylims=(-0.04,0.2), xlims=(0,85), size=(600,400))
    scatter!(sampled_ts, mean_extensions, title="Center $(monomer_ind*10) kb", ylabel="Relative extension" ,xlabel="Time [min]", color=:black)
    vline!([replicated_at], label="Replicated at", color=:red, linestyle=:dash)
    hline!([0], linecolor=:black, linewidth=0.5)

    #Then include an inset in the top right corner with a linlog scale
    scatter!(sampled_ts, mean_extensions.-fit.k[1], label="",ylabel="" ,xlabel="", color=:black, size=(700,400), yscale=:log10,xlims=(1,76), xticks=15:15:75,ylims=(1e-5,1),
        inset = (1, bbox(0.55,0.05,0.4,0.4)), subplot = 2)
    plot!(t->fit.p[1]*exp(fit.λ[1]*t), color=:black, subplot=2)
    vline!([replicated_at], color=:red, linestyle=:dash, subplot=2)
end

"""
    plot_extension_decay_after_replication(sampled_ts, extensions, stat_dir, plot_directory, region_widths)

Plot the extension decay after replication.

# Arguments
- `sampled_ts`: The time points.
- `extensions`: The extensions.
- `stat_dir`: The directory to save the statistics.
- `plot_directory`: The directory to save the plots.
- `region_widths`: The region widths.
"""
function plot_extension_decay_after_replication(sampled_ts, extensions, stat_dir, plot_directory, region_widths)
    N=Int(size(extensions)[1]/2)
    nt=length(sampled_ts)

    for region_width in region_widths #bins
        mid_points=range(region_width,N-region_width, step=region_width)
        exponents=zeros(length(mid_points))
        coefs_CIDs=zeros(length(mid_points))
        constants_CIDs=zeros(length(mid_points))
        mean_ext_t0=zeros(length(mid_points))
        strain_when_replicated=zeros(length(mid_points))
        time_replicated=zeros(length(mid_points))

        #plot extension of regions as a function of time
        for (index,mid_point) in enumerate(mid_points)
            if mid_point+region_width<=N/2 #right arm
                replicated_at=(mid_point+region_width)/replication_speed #in minutes
            else
                replicated_at=(N-mid_point-region_width)/(replication_speed)+5 #in minutes
            end
            time_replicated[index]=replicated_at
            replicated_from=findfirst(sampled_ts.>replicated_at)

            region=max(1, mid_point-region_width):min(N, mid_point+region_width)

            mean_strain=vec(mean(extensions[region,:], dims=1))
            mean_ext_t0[index]=mean_strain[1]

            if !isnothing(replicated_from) && replicated_from<nt-2
                fit=expsum_fit(sampled_ts[replicated_from:end], mean_strain[replicated_from:end], 1, withconst=true)
                exponents[index]=fit.λ[1]
                coefs_CIDs[index]=fit.p[1]
                constants_CIDs[index]=fit.k[1]
                
                strain_when_replicated[index]=fit.k+fit.p[1]*exp(fit.λ[1]*replicated_at)
                plot_strain_over_time(sampled_ts, mean_strain, fit, replicated_at, mid_point)
                savefig(plot_directory*"extension_midpoint_$(mid_point)_width_$(region_width).pdf")
            else
                exponents[index]=NaN
                coefs_CIDs[index]=NaN
                strain_when_replicated[index]=NaN
                constants_CIDs[index]=NaN
            end
        end

        #save the time scale, ie 1/exponent, for each region
        writedlm(stat_dir*"strain_exponents_region_width_$(region_width).txt", hcat(mid_points, exponents))
        writedlm(stat_dir*"strain_timescales_region_width_$(region_width).txt", hcat(mid_points, -1 ./exponents))
    end
end
