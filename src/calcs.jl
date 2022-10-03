using FFTW
using StaticArrays
using Statistics

"""
    msd_direct(traj::Trajectory, targetatom::String, resolution::Integer, timerange::Real)

Calculate the mean squared displacement (MSD) averaged over all atoms of type `targetatom` in trajectory `traj`.

Correlation times will be sampled in `resolution` equidistant points up until `timerange` times trajectory length. For a trajectory of length 100, a timerange of 0.3 and resolution of 10 would produce correlation times in 0:3:27.
"""
function msd_direct(traj::Trajectory, targetatom, resolution, timerange)
    if !traj.unwrapped
        println("Trajectory is not unwrapped -> MSD might be faulty.")
    end

    # Correlation time in ps, MSD in A^2
    τ, msd = msd_direct(traj.coords, traj.timestep_in_fs, traj.atomlabels, targetatom, resolution, timerange)
    return τ, msd
end

function msd_direct(coords, timestep_fs, atom, targetatom, resolution, timerange)
    timesteps = size(coords)[2]
    corr_time_list = range(0, length = resolution, step = round(Int, timesteps * timerange / resolution))
    coords_target = @view coords[atom .== targetatom, :]
    msd_result = Array{Float64, 1}(undef, resolution)
    Threads.@threads for i in 1:resolution
        msd_result[i] = msd_for_corr_time(coords_target, corr_time_list[i])
    end
    # correlation time in ps, MSD in A^2
    return corr_time_list * timestep_fs / 1000, msd_result
end

function msd_for_corr_time(coords, corr_time, unwrapped = false)
    a = 0
    atom_no, frame_no = size(coords)
    for atom in 1:atom_no
        for frame in 1:frame_no-corr_time
            @inbounds a += sum( x -> x^2, coords[atom, frame] .- coords[atom, frame + corr_time] )
        end
    end
    return a /= (frame_no - corr_time) * atom_no
end

"""
    msd_xyz_direct(traj::Trajectory, targetatom::String, resolution::Integer, timerange::Real)

Calculate the mean squared displacement (MSD) in x-, y- and z-direction averaged over all atoms of type `targetatom` in trajectory `traj`.

Correlation times will be sampled in `resolution` equidistant points up until `timerange` times trajectory length. For a trajectory of length 100, a timerange of 0.3 and resolution of 10 would produce correlation times in 0:3:27.
"""
function msd_xyz_direct(traj::Trajectory, targetatom, resolution, timerange)
    if !traj.unwrapped
        println("Trajectory is not unwrapped -> MSD might be faulty.")
    end

    # Correlation time in ps, MSD in A^2
    τ, msd = msd_xyz_direct(traj.coords, traj.timestep_in_fs, traj.atomlabels, targetatom, resolution, timerange)
    return τ, msd
end

function msd_xyz_direct(coords, timestep_fs, atom, targetatom, resolution, timerange)
    timesteps = size(coords)[2]
    corr_time_list = range(0, length = resolution, step = round(Int, timesteps * timerange / resolution))
    coords_target = @view coords[atom .== targetatom, :]
    msd_result = Array{Float64, 2}(undef, resolution, 3)
    Threads.@threads for i in 1:resolution
        msd_result[i, :] = msd_xyz_for_corr_time(coords_target, corr_time_list[i])
    end
    # correlation time in ps, MSD in A^2
    return corr_time_list * timestep_fs / 1000, msd_result
end

function msd_xyz_for_corr_time(coords, corr_time, unwrapped = false)
    a = @SVector [0, 0, 0]
    atom_no, frame_no = size(coords)
    for atom in 1:atom_no
        for frame in 1:frame_no-corr_time
            @inbounds a += (coords[atom, frame] .- coords[atom, frame + corr_time]).^2
        end
    end
    return a /= (frame_no - corr_time) * atom_no
end


function msd_fft(traj, timestep_fs, atom, targetatom, timerange)
    traj_msd = @view traj[:, atom .== targetatom, :]
    traj_msd = permutedims(traj_msd, [3, 1, 2])

    # compute S2
    fft_result = rfft(cat(dims = 1, traj_msd, zeros(Float64, size(traj_msd))), 1)
    fft_result = abs2.(fft_result) * size(traj_msd)[1]
    fft_result = irfft(fft_result, size(traj_msd)[1] * 2, 1)
end

"""
    rdf(traj::Trajectory, atomtype1::String, atomtype2::String, d_min::Real, d_max::Real, bins::Integer, central_atoms::Nothing)
    rdf(traj::Trajectory, atomtype1::String, atomtype2::String, d_min::Real, d_max::Real, bins::Integer, central_atoms::Vector{String})

Calculate the radial distribution function (RDF) between the two specified atom types in the trajectory.

If central_atoms is specified, then atoms will be assigned to molecules based on their closest central atom and only intermolecular distances are considered. Otherwise normal RDF is calculated.
"""
function rdf(traj::Trajectory, atomtype1, atomtype2, d_min, d_max, bins, central_atoms = nothing)
    coord, atom, pbc = traj.coords, traj.atomlabels, traj.mdbox.pbc_matrix
    # Each bin i contains distances from distance_limits[i] to distance_limits[i+1]; distnace_mean[i] will be the average of both values
    distance_limits = range(d_min, d_max, length = bins + 1)
    distance_mean = (distance_limits[1:end-1] + distance_limits[2:end]) / 2
    # Shell volume of each bin i and overall ideal density used to normalize RDF
    shell_volumes = 4/3 * pi * (distance_limits[2:end].^3 - distance_limits[1:end-1].^3)
    cell_volume = det(pbc) #pbc[1,:] ⋅ (pbc[2,:] × pbc[3,:])
    ideal_density = count(atom .== atomtype1) * count(atom .== atomtype2) / cell_volume
    # Create distance histogram; delegate calculation of all distances and binning to create_dist_histogram_multi
    dist_histogram = create_dist_histogram_multi(traj, atomtype1, atomtype2, distance_limits, central_atoms)
    rdf = dist_histogram ./ (shell_volumes .* ideal_density .* size(coord)[2])
    return distance_mean, rdf
end

function create_dist_histogram_multi(traj::Trajectory, atomtype1, atomtype2, distance_limits, central_atoms::Nothing)
    coord, atom, mdbox = traj.coords, traj.atomlabels, traj.mdbox
    dist_histogram = zeros(Int, distance_limits.len - 1)
    # Copy trajectories of specified atom types for later use; using @view here causes a lot more memory usage
    coord1 = coord[atom .== atomtype1, :]
    coord2 = coord[atom .== atomtype2, :]
    d_min, d_max, d_step = distance_limits[1], distance_limits[end], convert(Float64, distance_limits.step)
    # Split length of trajecgory $size(coord)[3] into roughly equal chunks, one for each available thread
    thread_limits = convert.(Int, floor.(range(1, size(coord)[2] + 1, length = Threads.nthreads() + 1)))
    # Lock to prevent two threads from editing dist_histogram
    Rlock = ReentrantLock()
    # Each thread will process one iteration of this loop; each will produce its own dist_histogram_thread and add it to the overall dist_histogram at the end
    Threads.@threads for i_thread in 1:Threads.nthreads()
        dist_histogram_thread = zeros(Int, distance_limits.len - 1)
        mdbox_thread = deepcopy(mdbox)
        i_limits = thread_limits[Threads.threadid()]:(thread_limits[Threads.threadid() + 1] - 1)
        # Iterate through chunk of timesteps set aside for this thread
        @views for i in i_limits
            # Iterate through atoms of both atom types
            for x in coord1[:, i]
                for y in coord2[:, i]
                    distance = pbc_dist(x, y, mdbox_thread)
                    bin!(dist_histogram_thread, distance, d_min, d_max, d_step)
                end
            end
        end
        # Add thread-local dist_histogram_thread to overall dist_histogram; lock it in the meantime in order to prevent threads from interfering with each other
        lock(Rlock)
        dist_histogram .+= dist_histogram_thread
        unlock(Rlock)
    end
    return dist_histogram
end

function create_dist_histogram_multi(traj::Trajectory, atomtype1, atomtype2, distance_limits, central_atoms)
   coord, atom, mdbox = traj.coords, traj.atomlabels, traj.mdbox
   dist_histogram = zeros(Int, distance_limits.len - 1)
   # Copy trajectories of specified atom types for later use; using @view here causes a lot more memory usage
   coord1 = coord[atom .== atomtype1, :]
   coord2 = coord[atom .== atomtype2, :]
   d_min, d_max, d_step = distance_limits[1], distance_limits[end], convert(Float64, distance_limits.step)
   # Split length of trajecgory $size(coord)[3] into roughly equal chunks, one for each available thread
   thread_limits = convert.(Int, floor.(range(1, size(coord)[2] + 1, length = Threads.nthreads() + 1)))
   # Lock to prevent two threads from editing dist_histogram
   Rlock = ReentrantLock()
   # Each thread will process one iteration of this loop; each will produce its own dist_histogram_thread and add it to the overall dist_histogram at the end
   Threads.@threads for i_thread in 1:Threads.nthreads()
       dist_histogram_thread = zeros(Int, distance_limits.len - 1)
       mdbox_thread = deepcopy(mdbox)
       molecules_1 = zeros(Int, size(coord1, 1))
       molecules_2 = zeros(Int, size(coord2, 1))
       i_limits = thread_limits[Threads.threadid()]:(thread_limits[Threads.threadid() + 1] - 1)
       # Iterate through chunk of timesteps set aside for this thread
       @views for i in i_limits
           frame1 = coord1[:, i]
           frame2 = coord2[:, i]
           frame_central = coord[ in(central_atoms).(atom), i]
           for (yi, y) in enumerate(frame2)
               molecules_2[yi] = next_neighbor(y, frame_central, mdbox)[1]
           end
           for (xi, x) in enumerate(frame1)
               molecules_1[xi] = next_neighbor(x, frame_central, mdbox)[1]
           end
           # Iterate through atoms of both atom types
           for (x, m1) in zip(frame1, molecules_1)
               for (y, m2) in zip(frame2, molecules_2)
                   if m1 != m2
                       distance = pbc_dist(x, y, mdbox_thread)
                       bin!(dist_histogram_thread, distance, d_min, d_max, d_step)
                   end
               end 
           end 
       end
       # Add thread-local dist_histogram_thread to overall dist_histogram; lock it in the meantime in order to prevent threads from interfering with each other
       lock(Rlock)
       dist_histogram .+= dist_histogram_thread
       unlock(Rlock)
   end 
   return dist_histogram
end

"""
    bin!(dist_histogram::Vector{R}, distance::Real, d_range::D) where {R <: Real, D <: StepRangeLen}
    bin!(dist_histogram::Vector{R}, distance::Real, d_min::Real, d_max::Real, d_step::Real)

Add a value `distance` to the appropriate bin in `dist_histogram` according to the limits in `d_range` or `d_min`, `d_max` and `d_step`.

# Examples
```jldoctest
julia> d_limits = range(1; stop = 5, length = 5)
1.0:1.0:10.0
julia> hist = zeros(Int, 10)
5-element Vector{Int64}:
 0
 0
 0
 0
 0
julia> TrajTools.bin!(hist, 2.5, d_limits)
1
julia> hist
10-element Vector{Int64}:
 0
 1
 0
 0
 0
```
"""
function bin!(dist_histogram, distance, d_range::D) where {D <: StepRangeLen}
    bin!(dist_histogram, distance, d_range[1], d_range[end], convert(Float64, d_range.step))
end

function bin!(dist_histogram, distance, d_min, d_max, d_step)
    # each bin i ∈ [1, bins] covers distances from (d_min + d_step * (i-1)) to (d_min + d_step * i)
    if d_min < distance < d_max
        index = floor(Int, (distance - d_min) / d_step + 1)
        dist_histogram[index] += 1
    end
end

"""
    autocorrFFT(inputarray::Abstractarray)

Autocorrelates the content of 3D array `inputarray` along the last axis.
"""
function autocorrFFT(inputarray)
    # Possible improvements: remove allocation of paddedarray
    #                        make it agnostic to dimensionality of inputarray
    arraylength = size(inputarray)[end]
    paddedlength = arraylength * 2 + 1
    paddedarray = zeros(eltype(inputarray), size(inputarray)[1:end-1]..., paddedlength)
    @views paddedarray[:, :, 1:arraylength] = inputarray

    FTresult = rfft(paddedarray, 3)
    FTresult .= FTresult .* conj(FTresult)
    autocorr = irfft(FTresult, paddedlength, 3)[:, :, 1:arraylength]
    for i in 1:arraylength
        @views autocorr[:, :, i] ./= (arraylength - i + 1)
    end
    return autocorr
end

"""
    oacf(traj::Trajectory, atomtype1::String, atomtype2::String)

Calculate the autocorrelation of the vectors connecting each atom of `atomtype1` to its nearest neighbor of `atomtype2` (orientational autocorrelation function, OACF).
"""
function oacf(traj::Trajectory, atomtype1, atomtype2)
    atomno1 = count(traj.atomlabels .== atomtype1)
    atomno2 = count(traj.atomlabels .== atomtype2)
    frame_no = size(traj.coords)[2]

    neighborsfirstframe = zeros(Int64, atomno1)
    frame1type2 = traj.coords[traj.atomlabels .== atomtype2, 1]
    @views for (i, pointtype1) in enumerate(traj.coords[traj.atomlabels .== atomtype1, 1])
        neighborsfirstframe[i] = TrajTools.next_neighbor(pointtype1, frame1type2, traj.mdbox)[1]
    end
    coordstype1 = @view traj.coords[traj.atomlabels .== atomtype1, :]
    coordstype2 = @view traj.coords[traj.atomlabels .== atomtype2, :]
    bondvectors = zeros(Float64, 3, atomno1, frame_no)
    @views for frame in 1:frame_no
        for atom1 in 1:count(traj.atomlabels .== atomtype1)
            # Find nearest neighbor of atomtype2 for each atomtype1
            neighborindex = neighborsfirstframe[atom1]
            coord1 = coordstype1[atom1, frame]
            coord2 = coordstype2[neighborindex, frame]
            # Store normalized vector connecting atomtype1 with its atomtype2 neighbor
            bondvectors[:, atom1, frame] .= TrajTools.minimum_image_vector(coord1, coord2, traj.mdbox)
            normalize!(bondvectors[:, atom1, frame])
        end
    end
    # Autocorrelate vectors over time (autocorrelations are independent of each other along all vectors and all spatial directions)
    auto1 = autocorrFFT(bondvectors)
    # orientational auto correlation -> sum over all spatial dimensions and average over all atoms
    oacf = vec(mean(sum(auto1, dims=1), dims = 2))
    # time in ps
    time_oacf = (1:frame_no) * traj.timestep_in_fs / 1E3
    return time_oacf, oacf
end
