#Name some colors used for plotting
using HDF5

#############################################
#Parsing and  other general helpers
subdirs(dir) = filter(x -> isdir(x), readdir(dir, join=true)) 

function find_first_or_end(s::String, chars)
    # Find the first occurrence of any of the characters
    index = findfirst(chars, s)
    # If no characters are found, return the length of the string + 1
    return index !== nothing ? index[1] : length(s) + 1
end

function find_first_or_end(s::String, chars)
    # Find the first occurrence of any of the characters
    index = findfirst(chars, s)
    
    # If no characters are found, return the length of the string + 1
    return index !== nothing ? index[1] : length(s) + 1
end

function parse_after(dir_name, what="/R_")
    from=findfirst(what, dir_name)[end]+1
    last_chars=r"(/|\.h5|_|\.txt)"

    to=from+find_first_or_end(dir_name[from:end], last_chars)-2
    try
        return parse(Int,dir_name[from:to])
    catch
        try
            return parse(Float64, dir_name[from:to])
        catch
            @error "Could not parse $what from $(dir_name)"
        end
    end
end

#############################################
#Navigating monomer indices

periodic_ind(i,N)= i>1 ? Int((i-1)%N+1) : periodic_ind(i+N,N)

replicated(fork,N) = vcat(max(fork...)+1:N, 1:min(fork...)-1)

function unreplicated(fork)
    return min(fork...)+1:max(fork...)-1
end

is_replicated(locus, fork)= (locus<min(fork...) || locus>max(fork...))

#############################################
# Shifting vectors and matrices so that ori is centered

#Shift map so that the middle of the map is at the middle of the matrix
function shifted_map(M,mid=1)
    N=size(M)[1]
    return circshift(M,(floor(Int,N/2)-mid, floor(Int,N/2)-mid))
end

function shifted_vector(v,mid=1)
    N=size(v)[1]
    return circshift(v,floor(Int,N/2)-mid) 
end

function shifted_vector_replicated(v,mid=1)
    N=Int(size(v)[1]/2)
    return vcat(shifted_vector(v[1:N],mid),shifted_vector(v[N+1:end],mid))
end

################################
# Converting between replicated length, time, cell length

function fork_to_R(fork,N)
    return sort(fork)[1]+N-sort(fork)[2]
end

function rel_monomer_density(fork,N,t)
    replicated_monomers=N+fork[1]-1+(N-fork[2])
    height=time_to_height(t)
    new_density=replicated_monomers/height
    original_density=N/time_to_height(0)
    return new_density/original_density
end

time_to_height(t,r_0=25.6*88/1000, growth_exponent=0.0055)=r_0*exp(growth_exponent*t)


################################
# Fetching data from MD simulation output

time_to_block(time, blocks_per_minute=1)=round(Int,time/blocks_per_minute)
block_to_time(block, blocks_per_minute=1)=block*blocks_per_minute

function fetch_time(n, file)
    read_attribute(file["$n"], "block")
    return block_to_time(block)
end

"""
Turn chromosome so that the ori further from the ter is lower, and labeled as the first replicated strand
This ensures that statistics are comparable to MaxEnt model
"""
function turn_chromosome!(array,fork,N=4050)
    @assert size(array)[2]==2*N
    ori_1=1
    ori_2=N+1
    ter=floor(Int,N/2)
    #check if need to relabel
    if norm(array[3,ori_1]-array[3,ter])<norm(array[3,ori_2]-array[3,ter]) #ori 1 is closer to ter than ori 2
        swapped=copy(array)
        s1=vcat(1:fork[1]-1, fork[2]+1:N)
        array[:,s1].=swapped[:,s1.+N]
        array[:,s1.+N].=swapped[:,s1]
    end
    #turn so that ori 1 is always in first cell half
    if array[3,ori_1]>0
        array[3,:].*=-1
    end
end

"""
This function calculates the time for each index based on the block number,
returns the indices for the time-points closest to those in array
'list_times'
"""
function get_indices_at_times(file, list_times)
    ks=keys(file)
    blocks=map(k-> read_attribute(file["$k"], "block"), ks)
    required_blocks=time_to_block.(list_times)
    list_inds=[]
    for required_block in required_blocks
        i=findall(blocks .== required_block)
        if !isempty(i)
            push!(list_inds, ks[i])
        else
            push!(list_inds, [])
        end
    end
    return list_inds
end

"""
Get the array for monomer positions
    - file : .h5 file to get data from
    - inds : list of indices (usually from get_indices_at_times) to get positions for
    - ME_orient : turn positions to match MaxEnt data?
    - N : Number of files
"""
function fetch_pos_fork_at_inds(file, inds; ME_orient=false, fork_orient=false)
    data=[]
    forks=[]
    for i in inds
        if isnothing(i)
            push!(data, nothing)
            push!(forks, nothing)
        else
            pos=read(file["$i/pos"])
            N=Int(size(pos)[2]/2)
            fork=fetch_fork_at_ind(file,i)
            unrepl=unreplicated(fork)
            pos[:,unrepl.+N].=pos[:,unrepl] #assign replicated positions to unreplicated monomers.
            if ME_orient
                turn_chromosome!(pos,fork,N)
            end
            if fork_orient
                fork_down_chromosome!(pos,fork,N)
            end
            push!(data,pos)
            push!(forks,fork)
        end
    end
    return data, forks
end

"""
Fetch position array for a single index
"""
function fetch_pos_fork_at_ind(file, ind; ME_orient=false, fork_orient=false)
    if !haskey(file, "$ind") || isnothing(ind)
        return nothing
    else
        pos=read(file["$ind/pos"])
        N=Int(size(pos)[2]/2)
        fork=fetch_fork_at_ind(file,ind)
        unrepl=unreplicated(fork)
        pos[:,unrepl.+N].=pos[:,unrepl]
        if ME_orient
            turn_chromosome!(pos,fork,N)
        end
        if fork_orient
            fork_down_chromosome!(pos,fork,N)
        end
        return pos, fork
    end
end

"""
For a set of indices, fetch fork positions
from .h5 file

These can then be used as indices for replicated and unreplicated
sections of the chromosome
"""
function fetch_fork_at_inds(file, inds)
    data=[]
    for i in inds
        try
            f=read(file["$i/Fork"])
        catch
            f=read(file["$i/fork"])
        end
        sort!(f)
        push!(data,f.+1)
    end
    return data
end

function fetch_fork_at_ind(file, ind)
    try
        f=read(file["$ind/Fork"])
        sort!(f)
        return f.+1
    catch
        f=read(file["$ind/fork"])
        sort!(f)
        return f.+1
    end
end

"""
Given .h5 file, fetch the positions of SMC links along chromosome
"""
function fetch_smcs_at_inds(file, inds)
    data=[]
    for i in inds
        if isnothing(i)
            push!(data, nothing)
        else
            push!(data,read(file["$i/SMCs"]).+1)
        end
    end
    return data
end

function fetch_smcs_at_ind(file, ind)
    return read(file["$ind/SMCs"]).+1
end

