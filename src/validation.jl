"""
This file contains scripts that take experimental data in csv files to calculate the average z-axis positions of given foci and their replicated copies

It then calculates the MaxEnt predictions for these same means and stds

Finally, we construct plots that compare the prediction and experiment
"""

using StatsBase, LinearAlgebra, HypothesisTests, DelimitedFiles
include("helpers.jl")

"""
    calc_mean_std_experiment(oclock, repl_time, exp_dir)

Calculate the mean and standard deviation of the z-axis positions of the given loci and their replicated copies

# Arguments
- `oclock`: the locus position
- `repl_time`: the time of replication in minutes
- `exp_dir`: the directory containing the experimental data

# Output
- `means`: a 2D array of mean z-axis positions, with the first column being the near locus and the second column being the far locus
- `stds`: a 2D array of standard deviations, with the first column being the near locus and the second column being the far locus
"""
function calc_mean_std_experiment(oclock, repl_time, exp_dir)
    data=readdlm(exp_dir*"$(oclock)oclock_ori.csv",',',skipstart=1)
    n_samples=size(data)[1]
    times=data[:,1]
    ls=data[:,2]

    times_good=[]
    near=[]
    far=[]

    for i in 1:n_samples
        pos=data[i,3:6]./ls[i]
        replicated=times[i]>=repl_time

        if replicated && !any(iszero,pos)
            low,low_i=findmin(pos)
            push!(near,low)
            push!(far,pos[(low_i+1)%4+1]) #11->21 etc
            push!(times_good,times[i])
        elseif !replicated
            p=minimum(pos[1:2])
            push!(near,p)
            push!(far,p)
            push!(times_good,times[i])
        end
    end

    ts=sort(unique(times_good))
    means=zeros(length(ts),2)
    stds=zeros(length(ts),2)
    counts=zeros(length(ts),2)

    for (index,t) in enumerate(ts)
        inds=findall(==(t),times_good)
        count=length(inds)
        means[index,:]=[mean(near[inds]), mean(far[inds])]
        stds[index,:]=[std(near[inds]),std(far[inds])]
        counts[index,:]=[count,count]
    end
    return means,stds,counts
end

"""
    calc_separations_experiment(oclock, repl_time, exp_dir)

Calculate the z-axis separations of the given loci and their replicated copies

# Arguments
- `oclock`: the locus position
- `repl_time`: the time of replication in minutes
- `exp_dir`: the directory containing the experimental data

# Output
- `ds`: a 1D array of z-axis separations
"""
function calc_separations_experiment(oclock, repl_time, exp_dir)
    if !(oclock in [0,12])
        data=readdlm(exp_dir*"$(oclock)oclock_ori.csv",',',skipstart=1)
        n_samples=size(data)[1]
        times=data[:,1]
        ls=data[:,2]

        times_good=[]
        ds=[]

        for i in 1:n_samples
            pos1=data[i,3]./ls[i]
            pos2=data[i,4]./ls[i]
            replicated=times[i]>=repl_time

            if replicated && !iszero(pos2) && !iszero(pos1)
                push!(ds,abs(pos1-pos2))
                push!(times_good,times[i])
            elseif !replicated
                push!(ds,0)
                push!(times_good,times[i])
            end
        end

        ts=sort(unique(times_good))
        means=zeros(length(ts))

        for (index,t) in enumerate(ts)
            inds=findall(==(t),times_good)
            means[index]=mean(ds[inds])
        end
        return means
    else
        data=readdlm(exp_dir*"$(oclock)oclock_mean_distances.txt",skipstart=1)
        times=data[:,1]
        ds=data[:,2]
        return ds
    end
end

"""
    calc_validation_stats_experiment(oclocks, replication_times,exp_dir, out_dir)

Calculate the mean and standard deviation of the z-axis positions of the given loci and their replicated copies

# Arguments
- `oclocks`: the locus positions
- `replication_times`: the times of replication in minutes
- `exp_dir`: the directory containing the experimental data
- `out_dir`: the directory to output the means and standard deviations

# Output
- nothing, but writes the means and standard deviations to files
"""
function calc_validation_stats_experiment(oclocks, replication_times,exp_dir, out_dir)
    out_exp=out_dir*exp_dir
    if !isdir(out_exp)
        mkpath(out_exp)
    end

    for (index,oclock) in enumerate(oclocks)
        means, stds, counts= calc_mean_std_experiment(oclock,replication_times[index],exp_dir)
        writedlm(out_exp*"$(oclock)_means.txt", means)
        writedlm(out_exp*"$(oclock)_stds.txt",stds)
        writedlm(out_exp*"$(oclock)_counts.txt",counts)
    end
end

"""
    calc_all_separations_experiment(oclocks, positions, total_length, replication_times,exp_dir, out_dir)

Calculate the z-axis separations of the given loci and their replicated copies

# Arguments
- `oclocks`: the locus positions
- `positions`: the locus positions in Mb
- `total_length`: the total length of the chromosome in Mb
- `replication_times`: the times of replication in minutes
- `exp_dir`: the directory containing the experimental data
- `out_dir`: the directory to output the means and standard deviations

# Output
- nothing, but writes the means to files
"""
function calc_all_separations_experiment(oclocks, positions, total_length, replication_times,exp_dir, out_dir)
    out_exp=out_dir*exp_dir
    if !isdir(out_exp)
        mkpath(out_exp)
    end

    for (index,oclock) in enumerate(oclocks)
        ds= calc_separations_experiment(oclock,replication_times[index],exp_dir)
        monomer_index=periodic_ind(round(Int,positions[index]*(405*4/total_length) +1220), 405*4)
        writedlm(out_exp*"$(oclock)_oclock_$(monomer_index)_monomer_index_separations.txt", ds)
    end
end

"""
    calc_mean_std_ME(ts, loci_me, replication_times, Ls, subdir, threshold_for_bias)

Calculate the mean and standard deviation of the z-axis positions of the given loci and their replicated copies

# Arguments
- `ts`: the times to calculate the mean and standard deviation
- `loci_me`: the loci to calculate the mean and standard deviation
- `replication_times`: the times of replication in minutes
- `Ls`: the simulation lengths
- `subdir`: the directory containing the simulation data
- `threshold_for_bias`: if the separation is below this, don't count it

# Output
- `me_means`: a 3D array of mean z-axis positions, with the first column being the near locus and the second column being the far locus
- `me_stds`: a 3D array of standard deviations, with the first column being the near locus and the second column being the far locus
- `me_means_cond`: a 3D array of mean z-axis positions, with the first column being the near locus and the second column being the far locus, conditioned on the separation being above `threshold_for_bias`
- `me_stds_cond`: a 3D array of standard deviations, with the first column being the near locus and the second column being the far locus, conditioned on the separation being above `threshold_for_bias`
"""
function calc_mean_std_ME(ts, loci_me, replication_times, Ls, subdir, threshold_for_bias)
    me_means=zeros(length(ts),2,length(loci_me))
    me_stds=zeros(length(ts),2,length(loci_me))
    me_means_cond=zeros(length(ts),2,length(loci_me))
    me_stds_cond=zeros(length(ts),2,length(loci_me))

    for (ind,t) in enumerate(ts)
        L=Ls[ind]

        lin_file=h5open(subdir*"/lin_t$t.h5", "r")
        ring_file=h5open(subdir*"/ring_t$t.h5", "r")

        num_samples=length(keys(lin_file))
        pos1=zeros(3,length(loci_me),num_samples)

        for i in 1:num_samples
            chr=read(ring_file,"$i")
            pos1[:,:,i].=chr[:,loci_me]
        end
        close(ring_file)
    
        if t>0
            pos2=zeros(3,length(loci_me),num_samples)
            for i in 1:num_samples
                chr=read(lin_file,"$i")
                for (li,l) in enumerate(loci_me)
                    if t>=replication_times[li]
                        pos2[:,li,i].=chr[:,l]
                    else
                        pos2[:,li,i].=pos1[:,li,i]
                    end
                end
            end
        else
            pos2=copy(pos1)
        end
        close(lin_file)

        #all positive
        pos1[3,:,:].+=L/2
        pos2[3,:,:].+=L/2
    
        for (j,locus) in enumerate(loci_me)
            near=[]
            far=[]
            near_cond=[]
            far_cond=[]
            for i in 1:num_samples
                if t<replication_times[j]
                    d1=min(pos1[3,j,i],L-pos1[3,j,i])/L
                    push!(near,d1)
                    push!(far,d1)
                    push!(near_cond,d1)
                    push!(far_cond,d1)
                else
                    d1=min(pos1[3,j,i],L-pos1[3,j,i])/L
                    d2=min(pos2[3,j,i],L-pos2[3,j,i])/L
                    push!(near,min(d1,d2))
                    push!(far,min(d1,d2)+abs(pos1[3,j,i]-pos2[3,j,i])/L)

                    if norm(pos1[2:3,j,i].-pos2[2:3,j,i])>threshold_for_bias #wouldn't be resolvable in experiment
                        push!(near_cond,min(d1,d2))
                        push!(far_cond,min(d1,d2)+abs(pos1[3,j,i]-pos2[3,j,i])/L)
                    end

                end
            end
            me_means[ind,1,j]=mean(near)
            me_means[ind,2,j]=mean(far)

            me_stds[ind,1,j]=std(near)
            me_stds[ind,2,j]=std(far)
        
            me_means_cond[ind,1,j]= !isempty(near_cond) ? mean(near_cond) : NaN
            me_means_cond[ind,2,j]= !isempty(far_cond) ? mean(far_cond) : NaN

            me_stds_cond[ind,1,j]= !isempty(near_cond) ? std(near_cond) : NaN
            me_stds_cond[ind,2,j]= !isempty(far_cond) ? std(far_cond) : NaN
        end

    end

    return me_means, me_stds, me_means_cond,me_stds_cond

end

"""
    calc_separations_ME(ts, Ls, loci_me, replication_times, subdir, threshold_for_bias)

Calculate the z-axis separations of the given loci and their replicated copies

# Arguments
- `ts`: the times to calculate the mean and standard deviation
- `Ls`: the simulation cell lengths
- `loci_me`: the loci to calculate the mean and standard deviation
- `replication_times`: the times of replication in minutes
- `subdir`: the directory containing the simulation data
- `threshold_for_bias`: if the separation is below this, don't count it

# Output
- `me_means`: a 2D array of mean z-axis separations, rows are times, columns are loci
- `me_means_cond`: a 2D array of mean z-axis separations, rows are times, columns are loci, conditioned on the separation being above `threshold_for_bias`
"""
function calc_separations_ME(ts, Ls, loci_me, replication_times, subdir, threshold_for_bias)
    me_means=zeros(length(ts),length(loci_me))
    me_means_cond=zeros(length(ts),length(loci_me))
    me_counts_cond=zeros(length(ts),length(loci_me))

    for (ind,t) in enumerate(ts)
        lin_file=h5open(subdir*"/lin_t$t.h5", "r")
        ring_file=h5open(subdir*"/ring_t$t.h5", "r")
        L=Ls[ind]
        num_samples=length(keys(lin_file))
        N=Int(size(read(lin_file,"1"),2))
    
        if t>0
            for i in 1:num_samples
                #load the z coordinate for each copy
                chr=read(ring_file,"$i")
                z1s=chr[3,loci_me]
                chr=read(lin_file,"$i")
                z2s=chr[3,loci_me]
                for (li,l) in enumerate(loci_me)
                    if t>=replication_times[li]
                        d=abs(z1s[li]-z2s[li])/L
                        me_means[ind,li]+=d
                        if d>threshold_for_bias/L
                            me_means_cond[ind,li]+=d
                            me_counts_cond[ind,li]+=1
                        end
                    end
                end
            end
        end
        close(ring_file)
        close(lin_file)
        me_means[ind,:]./=num_samples
    end
    me_means_cond./=me_counts_cond

    return me_means, me_means_cond
end

"""
    calc_validation_stats_MaxEnt(times, Ls, oclocks, locus_positions, replication_times, genomic_length, N, threshold_for_bias, data_dir, out_dir)

Calculate the mean and standard deviation of the z-axis positions of the given loci and their replicated copies

# Arguments
- `times`: the times to calculate the mean and standard deviation
- `Ls`: the simulation cell lengths
- `oclocks`: the locus positions
- `locus_positions`: the locus positions in Mb
- `replication_times`: the times of replication in minutes
- `genomic_length`: the total length of the chromosome in Mb
- `N`: the number of monomers
- `threshold_for_bias`: if the separation is below this, don't count it
- `data_dir`: the directory containing the simulation data
- `out_dir`: the directory to output the means and standard deviations

# Output
- nothing, but writes the means and standard deviations to files
"""
function calc_validation_stats_MaxEnt(times, Ls, oclocks, locus_positions, replication_times, genomic_length, N, threshold_for_bias, data_dir, out_dir)
    loci_me=round.(Int,locus_positions.*N/genomic_length)

    for subdir in subdirs(data_dir)
        println("Analyzing $subdir")
        out_sub_dir=replace(subdir, data_dir=>out_dir)
        if !isdir(out_sub_dir)
            mkpath(out_sub_dir)
        end

        me_means, me_stds, me_means_cond,me_stds_cond=calc_mean_std_ME(times, loci_me, replication_times, Ls, subdir, threshold_for_bias)

        for (index,oclock) in enumerate(oclocks)
            writedlm(out_sub_dir*"/raw_means_$oclock", me_means[:,:,index])
            writedlm(out_sub_dir*"/means_$oclock", me_means_cond[:,:,index])
            writedlm(out_sub_dir*"/raw_stds_$oclock", me_stds[:,:,index])
            writedlm(out_sub_dir*"/stds_$oclock", me_stds_cond[:,:,index])
        end
    end
end

"""
    calc_all_separations_MaxEnt(times, Ls, locus_positions, replication_times,genomic_length,N, threshold_for_bias, data_dir, out_dir)

Calculate the z-axis separations of the given loci and their replicated copies

# Arguments
- `times`: the times to calculate the mean separations
- `Ls`: the simulation cell lengths
- `locus_positions`: the locus positions in Mb
- `replication_times`: the times of replication in minutes
- `genomic_length`: the total length of the chromosome in Mb
- `N`: the number of monomers
- `threshold_for_bias`: if the separation is below this, don't count it
- `data_dir`: the directory containing the simulation data
- `out_dir`: the directory to output the means and standard deviations

# Output
- nothing, but writes the means to files
"""
function calc_all_separations_ME(times, Ls, locus_positions, replication_times,genomic_length,N, threshold_for_bias, data_dir, out_dir)
    loci_me=round.(Int,locus_positions.*N/genomic_length).+1

    for subdir in subdirs(data_dir) #loop over different simulations
        println("Analyzing mean distances for $subdir")
        out_sub_dir=replace(subdir, data_dir=>out_dir)
        if !isdir(out_sub_dir)
            mkpath(out_sub_dir)
        end
        #calculate 2D array (time, locus), with mean distances between replicates
        ds, ds_cond= calc_separations_ME(times, Ls, loci_me, replication_times, subdir, threshold_for_bias)
        difs=ds_cond.-ds
        difs=difs[.!isnan.(difs)]
        difs=difs[.!iszero.(difs)]
        mean_dif=mean(difs)
        writedlm(out_sub_dir*"/all_separations.txt", ds)
        writedlm(out_sub_dir*"/separations_cond.txt", ds_cond)
        writedlm(out_sub_dir*"/separations_difference_all_conditioned_mean_$(mean_dif).txt", ds_cond.-ds)
    end
end


"""
    calc_mean_std_LE(ts, loci_le, replication_times, Ls, subdir, threshold_for_bias)

Calculate the mean and standard deviation of the z-axis positions of the given loci and their replicated copies

# Arguments
- `ts`: the times to calculate the mean and standard deviation
- `loci_le`: the loci to calculate the mean and standard deviation
- `replication_times`: the times of replication in minutes
- `Ls`: the simulation lengths
- `subdir`: the directory containing the simulation data
- `threshold_for_bias`: if the separation is below this, don't count it

# Output
- `le_means`: a 3D array of mean z-axis positions, with the first column being the near locus and the second column being the far locus
- `le_stds`: a 3D array of standard deviations, with the first column being the near locus and the second column being the far locus
- `le_means_cond`: a 3D array of mean z-axis positions, with the first column being the near locus and the second column being the far locus, conditioned on the separation being above `threshold_for_bias`
- `le_stds_cond`: a 3D array of standard deviations, with the first column being the near locus and the second column being the far locus, conditioned on the separation being above `threshold_for_bias`
"""
function calc_mean_std_LE(ts, loci_le, monomer_size, subdir, threshold_for_bias, N=4050)
    le_means=zeros(length(ts),2,length(loci_le))
    le_stds=zeros(length(ts),2,length(loci_le))
    le_means_cond=zeros(length(ts),2,length(loci_le))
    le_stds_cond=zeros(length(ts),2,length(loci_le))
    num_samples=zeros(Int, length(ts), 2, length(loci_le))
    num_samples_resolved=zeros(Int, length(ts),2,length(loci_le))
    for filename in filter(contains(".h5"), readdir(subdir, join=true))
        h5open(filename, "r") do file
            list_inds=get_indices_at_times(file,times)
            for (ind,t) in enumerate(ts)
                for i in list_inds[ind]
                    positions, fork=fetch_pos_fork_at_ind(file,i, ME_orient=true, fork_orient=false)
                    pos1=positions[:,1:N]
                    pos2=positions[:,N+1:end]
                    L=time_to_height(t)/monomer_size
                    #all positive
                    pos1[3,:].+=L/2
                    pos2[3,:].+=L/2

                    for (j,locus) in enumerate(loci_le)
                        if !is_replicated(locus,fork) #only one copy
                            d1=min(pos1[3,locus],L-pos1[3,locus])/L

                            le_means[ind,1,j]+=d1
                            le_stds[ind,1,j]+=d1^2
                            num_samples[ind,1,j]+=1

                            le_means_cond[ind,1,j]+=d1
                            le_stds_cond[ind,1,j]+=d1^2
                            num_samples_resolved[ind,1,j]+=1
                        else #two copies, might or might not be resolved
                            d1=min(pos1[3,locus],L-pos1[3,locus])/L #locus 1 distances to poles
                            d2=min(pos2[3,locus],L-pos2[3,locus])/L #locus 2 distances to poles
                            d1=min(d1,d2)
                            d2=d1+abs(pos1[3,locus]-pos2[3,locus])/L
                            
                            le_means[ind,:,j].+=[d1,d2]
                            le_stds[ind,:,j].+=[d1,d2].^2
                            num_samples[ind,:,j].+=1

                            if norm(pos1[2:3,locus].-pos2[2:3,locus])>threshold_for_bias #resolved
                                le_means_cond[ind,:,j].+=[d1,d2]
                                le_stds_cond[ind,:,j].+=[d1,d2].^2
                                num_samples_resolved[ind,:,j].+=1
                            end
                        end
                    end
                end
            end
        end
    end
    #divide by number of samples to get means
    le_means_cond./=num_samples_resolved
    le_stds_cond./=num_samples_resolved #<x^2>
    le_stds_cond.-=le_means_cond.^2 #<x^2>-<x>^2

    le_means./=num_samples
    le_stds./=num_samples
    le_stds.-=le_means.^2

    return le_means, sqrt.(le_stds), le_means_cond, sqrt.(le_stds_cond)
end

"""
    calc_validation_stats_LE(times, Ls, oclocks, locus_positions, replication_times, genomic_length, N, threshold_for_bias, data_dir, out_dir)

Calculate the mean and standard deviation of the z-axis positions of the given loci and their replicated copies

# Arguments
- `times`: the times to calculate the mean and standard deviation
- `Ls`: the simulation cell lengths
- `oclocks`: the locus positions
- `locus_positions`: the locus positions in Mb
- `replication_times`: the times of replication in minutes
- `genomic_length`: the total length of the chromosome in Mb
- `N`: the number of monomers
- `threshold_for_bias`: if the separation is below this, don't count it
- `data_dir`: the directory containing the simulation data
- `out_dir`: the directory to output the means and standard deviations

# Output
- nothing, but writes the means and standard deviations to files
"""
function calc_validation_stats_LE(times, oclocks, locus_positions, genomic_length, threshold_for_bias_microns, data_dir, out_dir)

    for subdir in subdirs(data_dir)
        N=parse_after(subdir, "N_")
        monomer_size_microns=parse_after(subdir, "mon_size_")/1000

        loci_le=round.(Int,locus_positions.*N/genomic_length)
        threshold_for_bias=threshold_for_bias_microns/monomer_size_microns

        println("Analyzing $subdir")
        out_sub_dir=replace(subdir, data_dir=>out_dir)
        if !isdir(out_sub_dir)
            mkpath(out_sub_dir)
        end

        le_means, le_stds, le_means_cond,le_stds_cond=calc_mean_std_LE(times, loci_le, monomer_size_microns, subdir, threshold_for_bias, N)

        for (index,oclock) in enumerate(oclocks)
            writedlm(out_sub_dir*"/raw_means_$oclock", le_means[:,:,index])
            writedlm(out_sub_dir*"/means_$oclock", le_means_cond[:,:,index])
            writedlm(out_sub_dir*"/raw_stds_$oclock", le_stds[:,:,index])
            writedlm(out_sub_dir*"/stds_$oclock", le_stds_cond[:,:,index])
        end
    end
end

"""
    calc_p_value_mean(means, stds, exp_means, exp_stds)

Given a list of means and standard deviations, calculate the p-value of the experimental mean and standard deviation

# Arguments
- `means`: a 2D array of mean z-axis positions, with the first column being the near locus and the second column being the far locus
- `counts` : a 2D array of counts, with the first column being the near locus and the second column being the far locus
- `stds`: a 2D array of standard deviations, with the first column being the near locus and the second column being the far locus
- `exp_means`: aa 2D array of mean z-axis positions, with the first column being the near locus and the second column being the far locus

# Output
- `p_values`: an array of p-value of the means
"""
function calc_p_value_mean(means, counts, exp_means, exp_stds)
    #use OneSampleTTest to calculate p-values of the mean_zs
    p_values=pvalue.(OneSampleTTest.(exp_means, exp_stds, counts, means))
    return p_values
end

function collect_p_values(ts, oclocks, dir_exp, dir_LE, cond_av=true)
    all_p_values=zeros(length(ts), 2,length(oclocks))
    for (ind,oclock) in enumerate(oclocks)
        exp_m=readdlm(dir_exp*"$(oclock)_means.txt")[1:6,:]
        exp_std=readdlm(dir_exp*"$(oclock)_stds.txt")[1:6,:]
        exp_counts=Int.(readdlm(dir_exp*"$(oclock)_counts.txt")[1:6,:])

        if cond_av
            le_m=readdlm(dir_LE*"/means_$oclock")
        else
            le_m=readdlm(dir_LE*"/raw_means_$oclock")
        end

        unrepl=le_m[:,1].==le_m[:,2]

        all_p_values[:,:,ind].=calc_p_value_mean(le_m, exp_counts, exp_m, exp_std)
        all_p_values[unrepl,2,ind].=NaN
    end
    mean_p_val=mean(filter(!isnan,all_p_values))
    return all_p_values, mean_p_val
end
