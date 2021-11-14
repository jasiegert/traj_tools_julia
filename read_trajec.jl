using Pandas
using JLD

#path = "csh2po4_444_standard-pos-1.xyz"
#path = "small.xyz"

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
    # Read number of lines from first line 
    firstline = readline(path)
    noa = parse(Int, strip(firstline))

    # Read atom labels from first frame
    atom = open(path) do traj
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

function read_trajectory(path)
    auxiliarypath = path * ".jld"
    if isfile(auxiliarypath)
        println("Auxiliary file $auxiliarypath exists. Reading...")
        auxiliary = JLD.load(auxiliarypath)
        coord, atom = auxiliary["coord"], auxiliary["atom"]
        println("Coord: $(size(coord)), atom: $(size(atom))")
        return coord, atom
    elseif isfile(path)
        println("Auxiliary file $auxiliarypath doesn't exist. Reading xyz-file...")
        coord, atom = xyz_to_ar(path)
        println("Coord: $(size(coord)), atom: $(size(atom))")
        JLD.save(auxiliarypath, "coord", coord, "atom", atom)
        return coord, atom
    else
        println("Trajectory $path was not found.")
        return nothing, nothing
    end
end


path = ARGS[1]
coord, atom = read_trajectory(path)
if coord == nothing || atom == nothing
    println("Nothing is fine.")
else
    println("Everything is fine.")
end
