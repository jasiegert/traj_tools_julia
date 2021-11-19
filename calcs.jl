using FFTW
using Statistics

function msd_for_corr_time(traj, corr_time)
    start_coord = @view traj[(1 + corr_time):end, :, :]
    end_coord = @view traj[1:(end - corr_time), :, :]
    #mean(sum((end_coord - start_coord).^2, dims = 2))
    mean((end_coord - start_coord).^2) * 3
end

function msd_direct(traj, timestep_fs, atom, targetatom, resolution, timerange)
    timesteps = size(traj)[3]
    corr_time_list = range(0, length = resolution, step = round(Int, timesteps * timerange / resolution))
    traj_msd = traj[:, atom .== targetatom, :]
    traj_msd = permutedims(traj_msd, [3, 1, 2])
    corr_time_list * timestep_fs / 1000, [msd_for_corr_time(traj_msd, corr_time) for corr_time in corr_time_list]
end

function msd_fft(traj, timestep_fs, atom, targetatom, timerange)
    traj_msd = @view traj[:, atom .== targetatom, :]
    traj_msd = permutedims(traj_msd, [3, 1, 2])

    # compute S2
    fft_result = rfft(cat(dims = 1, traj_msd, zeros(Float64, size(traj_msd))), 1)
    fft_result = abs2.(fft_result) * size(traj_msd)[1]
    fft_result = irfft(fft_result, size(traj_msd)[1] * 2, 1)

end
