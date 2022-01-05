module TrajTools

export Trajectory
export dist_tmp, matmul_tmp

include("read_trajec.jl")
include("distances.jl")
include("calcs.jl")

const dist_tmp = zeros(MVector{3})
const matmul_tmp = zeros(MVector{3})

end
