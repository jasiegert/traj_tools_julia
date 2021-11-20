#using Pandas
using JLD
using LinearAlgebra

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
        )

     atom_masses = [atomic_mass_dict[entry] for entry in atom]
     atom_masses /= sum(atom_masses)

     return trajectory .- sum(trajectory .* repeat(atom_masses', 3, 1), dims = 2)
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
        coord = remove_com(coord, atom)
    end

    return coord, atom
end

function pbc_dist(point1, point2, pbc)
    inv_pbc = inv(pbc)
    dist = abs.(point1 .- point2)
    rel_dist = transpose(inv_pbc) * dist
    rel_dist -= floor.(rel_dist .+ 0.5)
    new_dist = transpose(pbc) * rel_dist
    return norm(new_dist)
end

function pbc_dist_inv!(dist_tmp, point1, point2, pbc, inv_pbc)
    dist_tmp .= abs.(point1 .- point2)
    dist_tmp = inv_pbc * dist_tmp
    dist_tmp .-= floor.(dist_tmp .+ 0.5)
    dist_tmp = pbc * dist_tmp
    return norm(dist_tmp)
end

function pbc_dist_inv_matmul!(dist_tmp, matmul_tmp, point1, point2, pbc, inv_pbc)
    dist_tmp .= abs.(point1 .- point2)
    mul!(matmul_tmp, inv_pbc,dist_tmp)
    matmul_tmp .-= floor.(matmul_tmp .+ 0.5)
    mul!(dist_tmp, pbc, matmul_tmp)
    return norm(dist_tmp)
end
