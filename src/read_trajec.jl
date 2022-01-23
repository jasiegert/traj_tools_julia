using JLD
using LinearAlgebra
using StaticArrays
using DelimitedFiles

abstract type MDBox end


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
struct Trajectory{BoxType<:MDBox}
    coords::Array{Float64, 3}
    atomlabels::Vector{String}
    mdbox::BoxType
    unwrapped::Bool
    timestep_in_fs::Float64
end

Trajectory(coords, atomlabels, box, unwrapped) = Trajectory(coords, atomlabels, box, unwrapped, 0.5)
Trajectory(coords, atomlabels, box) = Trajectory(coords, atomlabels, box, false)

struct OrthorhombicBox <: MDBox
    pbc_matrix::SMatrix{3, 3, Float64, 9}
    cell_parameters::SVector{3, Float64}
    dist_tmp::MVector{3, Float64}
end

OrthorhombicBox(pbc::AbstractArray) = OrthorhombicBox(SMatrix{3, 3, Float64, 9}(diagm(pbc)), SVector{3, Float64}(pbc), zeros(MVector{3, Float64}))

struct TriclinicBox <: MDBox
    pbc_matrix::SMatrix{3, 3, Float64, 9}
    inv_pbc_matrix::SMatrix{3, 3, Float64, 9}
    realspace_tmp::MVector{3, Float64}
    inversespace_tmp::MVector{3, Float64}
end

TriclinicBox(pbc::AbstractArray) = TriclinicBox(pbc, inv(pbc), zeros(MVector{3, Float64}), zeros(MVector{3, Float64}))


"""
    read_pbc(pbc_path)

Reads the periodic boundary conditions from a properly formatted text file.

The file pointed to by pbc_path has to be of either of two formats:
 - arbitrary cell: 3 lines corresponding to the 3 cell vectors, each containing 3 float values
 - orthorhombic cell: 3 lines containing one float value each or 1 line containing three float values
"""
function read_pbc(pbc_path)
    pbc = readdlm(pbc_path)
    # if 3 float values are found, turn them into a vector and return an orthorhombic box
    if size(pbc) == (1,3) || size(pbc) == (3,1)
        pbc = vec(pbc)
        return OrthorhombicBox(pbc)
    # if 3x3 float values of correct shape, return a triclinic box
    elseif size(pbc) == (3,3)
        # pbc has to be determined by 6 float values -> should be lower triangular:
        # x1  0  0
        # x2 y2  0
        # x3 y3 z3
        # note: pbc in mdbox will be transposed, thus upper triangular
        pbc = transpose(pbc)
        if pbc == UpperTriangular(pbc)
            return TriclinicBox(pbc)
        else
            error("pbc-file not formatted correctly -> 3x3 matrix has to be lower triangular")
        end
    else
        error("pbc-file not formatted correctly")
        return nothing
    end
end

"""
    read_trajectory(xyz_path, pbc_path, removecom = true)

Reads a trajectory from an xyz-file and its periodic boundary conditions from a text file.

For the format of the pbc-file, see [`read_pbc(pbc_path)`](@ref).
"""
function read_trajectory(xyz_path::String, pbc_path::String, com = true, unwrap = true)
    mdbox = read_pbc(pbc_path)
    coords, atomlabels = xyz_or_jld_to_arrays(xyz_path)
    
    if com
        remove_com!(coords, atomlabels)
    end

    if unwrap
        unwrap_trajectory!(coords, mdbox)
    end

    return Trajectory(coords, atomlabels, mdbox, unwrap)
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

"""
    remove_com!(coords, atomlabels)

Removes the center of mass (com) movement from the trajectory coords by setting the com to (0, 0, 0) in each frame.
"""
function remove_com!(coords, atom)
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
    for i in 1:size(coords)[3]
        mul!(com, coords[:, :, i], atom_masses)
        coords[:, :, i] .-= com
    end
    return coords
end

function unwrap_trajectory(coords, mdbox::MDBox)
    coords_unwrapped = deepcopy(coords)
    unwrap_trajectory!(coords_unwrapped, mdbox)
    return coords_unwrapped
end

function unwrap_trajectory!(coords, mdbox::MDBox)
    pbc = mdbox.pbc
    if size(pbc) == (3,)
        pbc = diagm(pbc)
    end
    inv_pbc = inv(pbc)

    to_be_wrapped = mdbox.inversespace_tmp
    for frame in 1:size(coords)[3]-1
        for atom in 1:size(coords)[2]
            @views mdbox.realspace_tmp .= coords[:, atom, frame + 1] .- coords[:, atom, frame]
            mul!(mdbox.inversespace_tmp, inv_pbc, mdbox.realspace_tmp)
            to_be_wrapped .= floor.(mdbox.inversespace_tmp .+ 0.5)
            mul!(mdbox.realspace_tmp, pbc, to_be_wrapped)
            coords[:, atom, frame + 1] .-= mdbox.realspace_tmp
        end
    end
end


function xyz_or_jld_to_arrays(path)
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

    return coord, atom
end
