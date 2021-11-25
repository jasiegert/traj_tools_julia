#using Pandas
using JLD
using LinearAlgebra
using StaticArrays
using DelimitedFiles

function read_pbc(pbc_path)
    pbc = transpose(readdlm(pbc_path))
    return SMatrix{3,3}(pbc)
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


function xyz_to_ar(path)
    # Read number of atoms from first line 
    firstline = readline(path)
    noa = parse(Int, strip(firstline))

    # Read atom labels from first frame
    atom = open(path) do traj
        # skip first two rows (number of atoms and comment-line)
        readline(traj)
        readline(traj)
        return [string(split(readline(traj))[1]) for i in 1:noa]
    end
    println(atom)

    # Read coordiantes from entire file; currently just a wrapper for Pandas read_csv function
    #coord = Array(read_csv(path, header = nothing, usecols=[1,2,3], delimiter = "\\s+", skiprows = (x) -> x % (noa + 2) in (0, 1)))
    #coord = reshape(coord, Int(length(coord)/noa/3), noa, 3)
    coord = xyz_to_coord(path, noa)
    
    return coord, atom
end

function remove_com(trajectory, atom)
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

     return trajectory .- sum(trajectory .* repeat(atom_masses', 3, 1), dims = 2)
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

function read_trajectory(path, com = true)
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
        coord = remove_com!(coord, atom)
    end

    return coord, atom
end

function max_distance_for_pbc_dist(pbc)
    volume = det(pbc)
    a, b, c = eachcol(pbc)
    dist1 = volume / norm(cross(a, b))
    dist2 = volume / norm(cross(b, c))
    dist3 = volume / norm(cross(c, a))
    return min(dist1, dist2, dist3) / 2
end

function pbc_dist(point1, point2, pbc, inv_pbc = inv(pbc), dist_tmp = zeros(MVector{3}), matmul_tmp = zeros(MVector{3}))
    dist_tmp .= point1 .- point2
    mul!(matmul_tmp, inv_pbc,dist_tmp)
    matmul_tmp .-= floor.(matmul_tmp .+ 0.5)
    mul!(dist_tmp, pbc, matmul_tmp)
    return norm(dist_tmp)
end

function pbc_dist_triclinic(point1, point2, pbc, inv_pbc = inv(pbc), dist_tmp = zeros(MVector{3}), matmul_tmp = zeros(MVector{3}))
    dist_tmp .= abs.(point1 .- point2)
    mul!(matmul_tmp, inv_pbc,dist_tmp)
    matmul_tmp .-= floor.(matmul_tmp .+ 0.5) .+ 1
    distance = pbc[1, 1] + pbc[2, 2] + pbc[3, 3]
    for i in 1:3
        for j in 1:3
            for k in 1:3
                mul!(dist_tmp, pbc, matmul_tmp)
                new_distance = norm(dist_tmp)
                if new_distance < distance
                    distance = new_distance
                end
                matmul_tmp[3] += 1
            end
            matmul_tmp[3] -= 3
            matmul_tmp[2] += 1
        end
        matmul_tmp[2] -= 3
        matmul_tmp[1] += 1
    end
    return distance
end

#function next_neighbor(point1::Array{Float64, 1}, group2::Array{Float64, 2}, pbc, inv_pbc = inv(pbc), dist_tmp = zeros(MVector{3}), matmul_tmp = zeros(MVector{3}))
function next_neighbor(point1, group2, pbc, inv_pbc = inv(pbc), dist_tmp = zeros(MVector{3}), matmul_tmp = zeros(MVector{3}))
    point2 = @view group2[:, 1]
    index, distance = 1, pbc_dist(point1, point2, pbc, inv_pbc, dist_tmp, matmul_tmp)
    #if size(group2)[2] > 1
        for i in 2:size(group2)[2]
            point2 = @view group2[:, i]
            new_distance = pbc_dist(point1, point2, pbc, inv_pbc, dist_tmp, matmul_tmp)
            if new_distance < distance
                distance = new_distance
                index = i
            end
        end
    #end
    return index, distance
end

function next_neighbor(group1::Array{Float64, 2}, group2, pbc, inv_pbc = inv(pbc), dist_tmp = zeros(MVector{3}), matmul_tmp = zeros(MVector{3}))
    atom_no_1 = size(group1)[2]
    indices = zeros(Int, atom_no_1)
    distances = zeros(Float64, atom_no_1)
    for i in 1:atom_no_1
        indices[i], distances[i] = next_neighbor(group1[:, i], group2, pbc, inv_pbc, dist_tmp, matmul_tmp)
    end
    return indices, distances
end

function next_neighbor(point1::Array{Float64, 1}, group2::Array{Float64, 1}, pbc, inv_pbc = inv(pbc), dist_tmp = zeros(MVector{3}), matmul_tmp = zeros(MVector{3}))
    return 1, pbc_dist(point1, point2, pbc, inv_pbc, dist_tmp, matmul_tmp)
end

function next_neighbor!(indices_out, distances_out, group1::Array{Float64, 2}, group2, pbc, inv_pbc = inv(pbc), dist_tmp = zeros(MVector{3}), matmul_tmp = zeros(MVector{3}))
    atom_no_1 = size(group1)[2]
    indices = zeros(Int, atom_no_1)
    distances = zeros(Float64, atom_no_1)
    for i in 1:atom_no_1
        indices_out[i], distances_out[i] = next_neighbor(group1[:, i], group2, pbc, inv_pbc, dist_tmp, matmul_tmp)
    end
    return nothing #indices, distances
end
