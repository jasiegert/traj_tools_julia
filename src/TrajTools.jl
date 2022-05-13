module TrajTools

export Trajectory
export dist_tmp, matmul_tmp
export pbc_dist, next_neighbor, minimum_image_vector

include("read_trajec.jl")
include("distances.jl")
include("calcs.jl")

const dist_tmp = zeros(MVector{3})
const matmul_tmp = zeros(MVector{3})

end
