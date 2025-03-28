using LinearAlgebra, StatsBase, Distances
include("helpers.jl") #helper functions for fetching time point data from time trajectory simulations


"""
    fetch_zs(pos, fork)

Fetch the z-positions of the monomers, including replicated regions based on fork position

# Arguments
- `pos`: 3x2N array of positions
- `fork`: 2-element array of fork positions

# Returns
- `zs`: 2N array of z-positions
"""
function fetch_zs(pos, fork)
    N=Int(size(pos)[2]/2)
    zs=pos[3,:]
    #unreplicated; just use coordinate of first
    zs[N+fork[1]+1:N+fork[2]-1].=pos[3,fork[1]+1:fork[2]-1]
    return zs
end

"""
    fetch_zs_resolution_limited(pos, fork, resolution_limit=3.6363)

Fetch the z-positions of the monomers, including replicated regions based on fork position, with a resolution limit

# Arguments
- `pos`: 3x2N array of positions
- `fork`: 2-element array of fork positions
- `resolution_limit`: resolution limit

# Returns
- `zs`: 2N array of z-positions
- `resolved`: 2N array of whether the monomer is resolved
"""
function fetch_zs_resolution_limited(pos, fork, resolution_limit=3.6363)
    N=Int(size(pos)[2]/2)
    zs=pos[3,:]
    replicated=vcat(1:fork[1], fork[2]:N)
    resolved=vcat(trues(N),falses(N))
    for focus in replicated
        replicate=focus+N
        if norm(pos[:,focus].-pos[:,replicate])>resolution_limit
            resolved[replicate]=true
        else
            resolved[focus]=false
        end
    end
    return zs, resolved
end

function fetch_sister_distances(pos,fork, N=405)
    replicated=union(1:max(1,fork[1]), min(N,fork[2]):N)
    unreplicated=max(1,fork[1]):min(N, fork[2])
    distances=zeros(N)
    distances[unreplicated].=NaN

    for i in replicated
        distances[i]=norm(pos[:,i].-pos[:,i+N])
    end
    
    return distances
end

"""
    fetch_contacts(d_map, fork, contact_r=2)


# Arguments
- `d_map`: distance matrix
- `fork`: replication forks
- `contact_r`: contact radius

# Returns
- `contacts`: 2Nx2N array of contacts
"""
function fetch_contacts(d_map, fork, contact_r=2)
    N=Int(size(d_map)[1]/2)
    contacts=zeros(size(d_map))
    contacts[d_map.<=contact_r].=1
    unrepl=unreplicated(fork)
    contacts[:,unrepl.+N].=0
    contacts[unrepl.+N,:].=0
    return contacts 
end

"""
    fetch_z_histograms(pos, fork, bins)

Fetch the z-histograms of the monomers, including replicated regions based on fork position

# Arguments
- `pos`: 3x2N array of positions
- `fork`: 2-element array of fork positions
- `bins`: edges for positional bins along the cell length

# Returns
- `hist1`: histogram of monomer positions along cell length, replicated strand 1
- `hist2`: histogram of monomer positions along cell length, replicated strand 2 
- `hist3`: histogram of monomer positions along cell length, unreplicated region
"""
function fetch_z_histograms(pos, fork, bins)    
    N=Int(size(pos)[2]/2)
    repl=replicated(fork,N)
    unrepl=fork[1]+1:fork[2]-1
    total_mon=N+length(repl)
    zs_1=pos[3,repl]
    zs_2=pos[3,repl.+N]
    zs_3=pos[3, unrepl]
    hist1=StatsBase.fit(Histogram, zs_1, bins).weights ./total_mon
    hist2=StatsBase.fit(Histogram, zs_2, bins).weights ./total_mon
    hist3=StatsBase.fit(Histogram, zs_3, bins).weights ./total_mon
    return hist1, hist2, hist3
end

"""
    fetch_chipseq(smcs, N=4050)

Fetch the ChIP-seq signal for loop-extruders; an array where the ith element is number of loop-extruder legs at the ith monomer

# Arguments
- `smcs`: array of SMC indices
- `N`: number of monomers

# Returns
- `chipseq`: N array of ChIP-seq signal
"""
function fetch_chipseq(smcs, N=4050)
    inds=reshape(smcs, :,1)
    inds[inds.>N].-=N
    chipseq=zeros(N)
    chipseq[inds].+=1
    return chipseq
end
function fetch_loop_size_distribution(smcs, loop_size_bins, N=4050, total_kb=4050)
    proj_smcs=copy(smcs)
    @assert all(proj_smcs.>0) "got negative smc positions"
    proj_smcs[smcs.>N].-=N #project onto 1st chromosome
    @assert all(proj_smcs.<=N)

    loop_sizes=abs.(proj_smcs[1,:].-proj_smcs[2,:])
    loop_sizes[loop_sizes.>N/2].= N .- loop_sizes[loop_sizes.>N/2]

    #in kb instead of bins
    loop_sizes=loop_sizes.*total_kb/N

    loop_size_hist=StatsBase.fit(Histogram, loop_sizes, loop_size_bins).weights

    #return normalized by number of loops
    return loop_size_hist./sum(loop_size_hist)
end

"""
    fetch_local_extensions(pos, fork, d=2)

Fetch the local extensions of the monomers, including replicated regions based on fork position

# Arguments
- `pos`: 3x2N array of positions
- `fork`: 2-element array of fork positions
- `d`: genomic distance for which we measure local extension

# Returns
- `local_extensions`: 2N array of local extensions
"""
function fetch_local_extensions(pos, fork, d=2)
    N=Int(size(pos)[2]/2)
    local_extensions=zeros(2*N)
    for i in 1:N
        local_extensions[i]=norm(pos[:,periodic_ind(i-d,N)].-pos[:,periodic_ind(i+d,N)])
    end
    local_extensions[N+1:2*N].=local_extensions[1:N]
    for i in N+1:N+fork[1]-d
        local_extensions[i]=norm(pos[:,periodic_ind(i-d,N)+N].-pos[:,periodic_ind(i+d,N)+N])
    end
    for i in N+fork[2]+d+1:2*N
        local_extensions[i]=norm(pos[:,periodic_ind(i-d,N)+N].-pos[:,periodic_ind(i+d,N)+N])
    end
    return local_extensions
end

function relative_extensions(extensions, index=1)
    rel_extensions=zeros(size(extensions))
    nt=size(extensions,2)
    for i in 1:nt
        rel_extensions[:,i].=(extensions[:,i].-extensions[:,index])./extensions[:,index]
    end
    return rel_extensions
end

"""
Given an array of local extensions over time, calculate the relative extensions,
smoothen them using a moving average, and return the smoothened extensions as a kymograph
"""
function extension_kymograph(local_extensions, smoothening_width=5)
    nt=size(local_extensions)[2]
    N=Int(size(local_extensions)[1]/2)

    #First, calculate the relative extension
    rel_extensions=relative_extensions(local_extensions)

    #Then, smoothen the curves
    smoothened_extensions=zeros((N,nt))
    for t in 1:nt
        smoothened_extensions[:,t]=map(i->mean(rel_extensions[(i-smoothening_width+N)%N+1:(i+smoothening_width)%N+1,t]), 1:N)
    end

    return smoothened_extensions 
end




"""
    fetch_distance_matrix(pos, fork)

Fetch the distance matrix of the monomers, including replicated regions based on fork position

# Arguments
- `pos`: 3x2N array of positions
- `fork`: 2-element array of fork positions

# Returns
- `d_matrix`: 2N x 2N array of distances
"""
function fetch_distance_matrix(pos, fork)
    N=Int(size(pos)[2]/2)
    not_replicated=N+fork[1]+1:N+fork[2]-1

    d_matrix =pairwise(Euclidean(), pos, pos, dims=2)
    d_matrix[not_replicated,1:N].=d_matrix[not_replicated.-N,1:N]
    d_matrix[1:N,not_replicated].=d_matrix[1:N,not_replicated.-N]
    d_matrix[not_replicated,N+1:2*N].=d_matrix[not_replicated.-N,1:N]
    d_matrix[N+1:2*N,not_replicated].=d_matrix[1:N,not_replicated.-N]
    return d_matrix
end

"""
    fetch_all_spring_lengths(pos, bin_ds, incr=10)

For comparing the MD simulation 10 kb distance to MaxEnt model's constraints. Basically, 10 kb genomic distance should have a 129 nm mean 3D distance.

# Arguments
- `pos`: 3x2N array of positions
- `bin_ds`: bin edges for spring lengths
- `incr`: monomer number increment for spring lengths; depends on coarse-graining scale of MD simulation

# Returns
- `all_lengths`: 2N array of spring lengths
"""
function fetch_all_spring_lengths(pos, bin_ds, incr=10)
    N=size(pos,2)
    all_lengths= map(i->norm(pos[:,periodic_ind(i+incr,N)].-pos[:,i]), 1:N)
    return map(i->searchsortedfirst(bin_ds, i), all_lengths) #bin number to increment in
end


"""
    fetch_radial_positions(pos, fork, N=405)

Fetch the radial positions of the monomers, including replicated regions based on fork position

# Arguments
- `pos`: 3x2N array of positions
- `fork`: 2-element array of fork positions

# Returns
- `radial_positions`: 2N array of radial positions
"""
function fetch_radial_positions(pos, fork)
    N=Int(size(pos)[2]/2)
    radial_positions=zeros(2*N)
    replicated=union(N+1:max(1,fork[1]-1+N), min(2*N,fork[2]+1+N):2*N)
    for i in 1:N
        radial_positions[i]=norm(pos[1:2,i])
    end
    radial_positions[N+1:2*N].=radial_positions[1:N]
    for i in replicated
        radial_positions[i]=norm(pos[1:2,i])
    end
    return radial_positions
end

"""
    fetch_seg_fraction(pos, fork)

Measure used in Harju et al. 2024 to quantify the segregation of chromosomes within the cell. Gives a value between 0 and 1.

# Arguments
- `pos`: 3x2N array of positions
- `fork`: 2-element array of fork positions

# Returns
- `seg_fraction`: fraction of monomers in the "right" cell half, minus fraction in the "wrong" cell half
"""
function fetch_seg_fraction(pos, fork)
    N=Int(size(pos)[2]/2)
    repl1=vcat(1:fork[1]-1, fork[2]+1:N)
    N_repl=length(repl1)

    if N_repl<=0
        return 0
    end

    repl2=repl1.+N
    
    below1=length(findall(pos[3,repl1].<0))
    below2=length(findall(pos[3,repl2].<0))

    above2=length(findall(pos[3,repl2].>0))
    above1=length(findall(pos[3,repl1].>0))

    return max((2*below1+2*above2-2*N_repl)/(2*N_repl), (2*below2+2*above1-2*N_repl)/(2*N_repl))
end
