using FFTW
using StaticArrays

function msd_for_corr_time(traj_msd, corr_time, unwrapped = false)
    a = 0
    for k in 1:size(traj_msd)[3]
        for j in 1:size(traj_msd)[2]
            for i in 1:size(traj_msd)[1]-corr_time
                a += (traj_msd[i+corr_time, j, k] - traj_msd[i, j, k])^2
            end
        end
    end
    return a /= (size(traj_msd)[1] - corr_time) * size(traj_msd)[3]
end

function msd_direct(traj::Trajectory, targetatom, resolution, timerange)
    if !traj.unwrapped
        println("Trajectory is not unwrapped -> MSD might be faulty.")
    end

    return msd_direct(traj.coords, traj.timestep_in_fs, traj.atomlabels, targetatom, resolution, timerange)
end

function msd_direct(coords, timestep_fs, atom, targetatom, resolution, timerange)
    timesteps = size(coords)[3]
    corr_time_list = range(0, length = resolution, step = round(Int, timesteps * timerange / resolution))
    # Permuting the trajectory to shape (timesteps, dimensions, atoms) is faster, because it allows the innermost loop over the timesteps to be along the first dimension
    traj_msd = permutedims(coords[:, atom .== targetatom, :], [3, 1, 2])
    msd_result = Array{Float64, 1}(undef, resolution)
    Threads.@threads for i in 1:resolution
        msd_result[i] = msd_for_corr_time(traj_msd, corr_time_list[i])
    end
    # correlation time in ps, MSD in A^2
    return corr_time_list * timestep_fs / 1000, msd_result
end

function msd_fft(traj, timestep_fs, atom, targetatom, timerange)
    traj_msd = @view traj[:, atom .== targetatom, :]
    traj_msd = permutedims(traj_msd, [3, 1, 2])

    # compute S2
    fft_result = rfft(cat(dims = 1, traj_msd, zeros(Float64, size(traj_msd))), 1)
    fft_result = abs2.(fft_result) * size(traj_msd)[1]
    fft_result = irfft(fft_result, size(traj_msd)[1] * 2, 1)
end

function rdf(traj::Trajectory, atomtype1, atomtype2, d_min, d_max, bins)
    coord, atom, pbc = traj.coords, traj.atomlabels, traj.mdbox.pbc_matrix
    # Each bin i contains distances from distance_limits[i] to distance_limits[i+1]; distnace_mean[i] will be the average of both values
    distance_limits = range(d_min, d_max, length = bins + 1)
    distance_mean = (distance_limits[1:end-1] + distance_limits[2:end]) / 2
    # Shell volume of each bin i and overall ideal density used to normalize RDF
    shell_volumes = 4/3 * pi * (distance_limits[2:end].^3 - distance_limits[1:end-1].^3)
    cell_volume = det(pbc) #pbc[1,:] ⋅ (pbc[2,:] × pbc[3,:])
    ideal_density = count(atom .== atomtype1) * count(atom .== atomtype2) / cell_volume
    # Create distance histogram; delegate calculation of all distances and binning to create_dist_histogram_multi
    dist_histogram = create_dist_histogram_multi(traj, atomtype1, atomtype2, distance_limits)
    rdf = dist_histogram ./ (shell_volumes .* ideal_density .* size(coord)[3])
    return distance_mean, rdf
end

function create_dist_histogram_multi(traj::Trajectory, atomtype1, atomtype2, distance_limits)
    coord, atom, mdbox = traj.coords, traj.atomlabels, traj.mdbox
    dist_histogram = zeros(Int, distance_limits.len - 1)
    # Copy trajectories of specified atom types for later use; using @view here causes a lot more memory usage
    coord1 = coord[:, atom .== atomtype1, :]
    coord2 = coord[:, atom .== atomtype2, :]
    d_min, d_max, d_step = distance_limits[1], distance_limits[end], convert(Float64, distance_limits.step)
    # Split length of trajecgory $size(coord)[3] into roughly equal chunks, one for each available thread
    thread_limits = convert.(Int, floor.(range(1, size(coord)[3] + 1, length = Threads.nthreads() + 1)))
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
            for x in eachcol(coord1[:, :, i])
                for y in eachcol(coord2[:, :, i])
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

function bin!(dist_histogram, distance, d_min, d_max, d_step)
    # each bin i ∈ [1, bins] covers distances from (d_min + d_step * (i-1)) to (d_min + d_step * i)
    if d_min < distance < d_max
        index = floor(Int, (distance - d_min) / d_step + 1)
        dist_histogram[index] += 1
    end
end

function oacf(traj::Trajectory, atomtype1, atomtype2)
   neighborsfirstframe = zeros(Int64, count(traj.atomlabels .== atomtype1))
   frame1type2 = traj.coords[:, traj.atomlabels .== atomtype2, 1]
   @views for (i, pointtype1) in enumerate(eachcol(traj.coords[:, traj.atomlabels .== atomtype1, 1]))
       neighborsfirstframe[i] = TrajTools.next_neighbor(pointtype1, frame1type2, traj.mdbox)[1]
   end
   coordstype1 = @view traj.coords[:, traj.atomlabels .== atomtype1, :]
   coordstype2 = @view traj.coords[:, traj.atomlabels .== atomtype2, :]
   bondvectors = zeros(Float64, 3, count(traj.atomlabels .== atomtype1), size(traj.coords)[3])
   @views for frame in 1:size(traj.coords)[3]
       for atom1 in 1:count(traj.atomlabels .== atomtype1)
           neighborindex = neighborsfirstframe[atom1]
           coord1 = coordstype1[:, atom1, frame]
           coord2 = coordstype2[:, neighborindex, frame]
           bondvectors[:, atom1, frame] .= coord1 .- coord2
       end
   end
end
