using LinearAlgebra
using StaticArrays

function max_distance_for_pbc_dist(traj::Trajectory)
    max_distance_for_pbc_dist(traj.mdbox)
end

function max_distance_for_pbc_dist(mdbox::MDBox)
    max_distance_for_pbc_dist(mdbox.pbc)
end

function max_distance_for_pbc_dist(pbc::AbstractArray)
    volume = det(pbc)
    a, b, c = eachcol(pbc)
    dist1 = volume / norm(cross(a, b))
    dist2 = volume / norm(cross(b, c))
    dist3 = volume / norm(cross(c, a))
    return min(dist1, dist2, dist3) / 2
end

function pbc_dist(point1, point2, mdbox::TriclinicBox)
    pbc_dist(point1, point2, mdbox.pbc, mdbox.inv_pbc, mdbox.realspace_tmp, mdbox.inversespace_tmp)
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
