using JLD
using LinearAlgebra
using StaticArrays
using DelimitedFiles


"""
    Trajectory(coords, atomlabels, pbc, unwrapped = false, timestep_in_fs = 0.5, inv_pbc = inv(pbc))

Combines trajectory information read from xyz-file.

# Arguments
- `coords::Array{Float64, 3}`: coordinates
- `atomlabels::Vector{String}`: atom labels read from the first frame
- `pbc::StaticArrays.SMatrix{3, 3, Float64, 9}`: periodic boundary conditions as a 3x3 StaticArray
- `unwrapped::Bool`: whether or not the trajectory has been unwrapped
- `timestep_in_fs::Real`: length of each time step in fs
- `inv_pbc::StaticArrays.SMatrix{3, 3, Float64, 9}`: inverse of the pbc (stored seperately for performance reasons)
"""
struct Trajectory
    coords::Array{Float64, 3}
    atomlabels::Vector{String}
    pbc::StaticArrays.SMatrix{3, 3, Float64, 9}
    unwrapped::Bool
    timestep_in_fs::Real
    inv_pbc::StaticArrays.SMatrix{3, 3, Float64, 9}
end

Trajectory(coords, atomlabels, pbc, unwrapped, timestep_in_fs) = Trajectory(coords, atomlabels, pbc, unwrapped, timestep_in_fs, inv(pbc))
Trajectory(coords, atomlabels, pbc, unwrapped) = Trajectory(coords, atomlabels, pbc, unwrapped, 0.5)
Trajectory(coords, atomlabels, pbc) = Trajectory(coords, atomlabels, pbc, true)

"""
    read_pbc(pbc_path)

Reads the periodic boundary conditions from a properly formatted text file.

Periodic boundary conditions have to be represented as 3 lines corresponding to the 3 cell vectors, each containing 3 float values.
pbc_path is a string containing the path to the pbc-file.
"""
function read_pbc(pbc_path)
    pbc = transpose(readdlm(pbc_path))
    return SMatrix{3,3}(pbc)
end

function read_trajectory(xyz_path::String, pbc_path::String, com = true)
    pbc = read_pbc(pbc_path)
    return read_trajectory(xyz_path, pbc, com)
end

function read_trajectory(xyz_path::String, pbc::AbstractArray, com = true)
    coord, atom = xyz_or_jld_to_arrays(xyz_path, com)
    traj = Trajectory(coord, atom, pbc)
    return traj
end

function xyz_to_coord(path, noa)
    frame_no = Int( countlines(path) / (noa + 2) )
    return open(path) do traj
           coord = Array{Float64, 3}(undef, 3, noa, frame_no)
           for i in 1:frame_no
               readline(traj)
               readline(traj)
               for j in 1:noa
                   coord[:, j, i] = parse.(Float64, split(readline(traj))[2:end])
               end
           end
           return coord
     end
end

function xyz_to_atomlabels(path, noa)
    atom = open(path) do traj
        # skip first two rows (number of atoms and comment-line)
        readline(traj)
        readline(traj)
        return [string(split(readline(traj))[1]) for i in 1:noa]
    end
    return atom
end

function xyz_to_ar(path)
    # Read number of atoms from first line 
    firstline = readline(path)
    noa = parse(Int, strip(firstline))

    # Read atom labels from first frame

    atomlabels = xyz_to_atomlabels(path, noa)
    coord = xyz_to_coord(path, noa)
    
    return coord, atomlabels
end

function remove_com!(trajectory, atom)
    atomic_mass_dict = Dict{String, Float64}(
        "H" => 1.0079,
        "Li" => 6.941,
        "C" => 12.01,
        "N" => 14.0067,
        "O" => 15.9994,
        "F" => 18.9984,
        "Si" => 28.0855,
        "P" => 30.9738,
        "S" => 32.065,
        "Cs" => 132.9055,
        "Sn" => 118.71,
    )

    atom_masses = [atomic_mass_dict[entry] for entry in atom]
    atom_masses /= sum(atom_masses)

    com = zeros(MVector{3})
    for i in 1:size(trajectory)[3]
        mul!(com, trajectory[:, :, i], atom_masses)
        trajectory[:, :, i] .-= com
    end
    return trajectory
end

function xyz_or_jld_to_arrays(path, com = true)
    auxiliarypath = path * ".jld"
    if isfile(auxiliarypath)
        println("Auxiliary file $auxiliarypath exists. Reading...")
        auxiliary = JLD.load(auxiliarypath)
        coord, atom = auxiliary["coord"], auxiliary["atom"]
        println("Coord: $(size(coord)), atom: $(size(atom))")
    elseif isfile(path)
        println("Auxiliary file $auxiliarypath doesn't exist. Reading xyz-file...")
        coord, atom = xyz_to_ar(path)
        println("Coord: $(size(coord)), atom: $(size(atom))")
        JLD.save(auxiliarypath, "coord", coord, "atom", atom)
    else
        println("Trajectory $path was not found.")
        return nothing, nothing
    end

    if com == true
        remove_com!(coord, atom)
    end

    return coord, atom
end

