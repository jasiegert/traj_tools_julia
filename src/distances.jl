using LinearAlgebra
using StaticArrays

function max_distance_for_pbc_dist(traj::Trajectory)
    max_distance_for_pbc_dist(traj.mdbox)
end

function max_distance_for_pbc_dist(mdbox::MDBox)
    max_distance_for_pbc_dist(mdbox.pbc_matrix)
end

function max_distance_for_pbc_dist(pbc::AbstractArray)
    volume = det(pbc)
    a, b, c = eachcol(pbc)
    dist1 = volume / norm(cross(a, b))
    dist2 = volume / norm(cross(b, c))
    dist3 = volume / norm(cross(c, a))
    return min(dist1, dist2, dist3) / 2
end

function pbc_dist(point1, point2, mdbox::OrthorhombicBox)
    for i in 1:3
        mdbox.dist_tmp[i] = (point1[i] - point2[i]) / mdbox.cell_parameters[i]
        mdbox.dist_tmp[i] -= floor(mdbox.dist_tmp[i] + 0.5)
        mdbox.dist_tmp[i] *= mdbox.cell_parameters[i]
    end
    return norm(mdbox.dist_tmp)
end

function pbc_dist(point1, point2, mdbox::TriclinicBox)
    pbc_dist(point1, point2, mdbox.pbc_matrix, mdbox.inv_pbc_matrix, mdbox.realspace_tmp, mdbox.inversespace_tmp)
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

function next_neighbor(point1, group2, mdbox::MDBox)
    point2 = @view group2[:, 1]
    index, distance = 1, pbc_dist(point1, pointi2, mdbox)
    for i in 2:size(group2)[2]
        point2 = @view group2[:, i]
        new_distance = pbc_dist(point1, point2, mdbox)
        if new_distance < distance
            distance = new_distance
            index = i
        end
    end
    return index, distance
end
