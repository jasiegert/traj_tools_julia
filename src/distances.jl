using LinearAlgebra
using StaticArrays

"""
    max_distance_for_pbc_dist(traj::Trajectory)
    max_distance_for_pbc_dist(mdbox::MDBox)
    max_distance_for_pbc_dist(pbc::AbstractArray)

Calculates half of the shortest box diameter. For distances above this value, [`pbc_dist(point1, point2, mdbox)`](@ref) might give results that are too large.
"""
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

"""
    pbc_dist(point1, point2, mdbox)

Calculates the distance between two points in a simulation box according to the minimum image convention.

For non-orthorhombic boxes, this is only accurate up to half of the shortest box diameter. This limit can be calculated with [`max_distance_for_pbc_dist(traj::Trajectory)`](@ref).

# Arguments
    - point1::AbstractVector: coordinates in a vector of length 3
    - point2::AbstractVector: coordinates in a vector of length 3
    - mdbox::MDBox: simulation box
"""
function pbc_dist(point1, point2, mdbox::OrthorhombicBox)
    for i in 1:3
        mdbox.dist_tmp[i] = (point1[i] - point2[i]) / mdbox.cell_parameters[i]
        mdbox.dist_tmp[i] -= floor(mdbox.dist_tmp[i] + 0.5)
        mdbox.dist_tmp[i] *= mdbox.cell_parameters[i]
    end
    return norm(mdbox.dist_tmp)
end

#function pbc_dist(point1, point2, mdbox::TriclinicBox)
#    pbc_dist(point1, point2, mdbox.pbc_matrix, mdbox.inv_pbc_matrix, mdbox.realspace_tmp, mdbox.inversespace_tmp)
#end

#function pbc_dist(point1, point2, pbc, inv_pbc = inv(pbc), dist_tmp = zeros(MVector{3}), matmul_tmp = zeros(MVector{3}))
#    dist_tmp .= point1 .- point2
#    mul!(matmul_tmp, inv_pbc,dist_tmp)
#    matmul_tmp .-= round.(matmul_tmp)
#    mul!(dist_tmp, pbc, matmul_tmp)
#    return norm(dist_tmp)
#end

function pbc_dist(p1, p2, mdbox::TriclinicBox)
    pbc_dist(p1, p2, mdbox.pbc_matrix, mdbox.inv_pbc_matrix)
end

# ugliest implementation yet, but also fast
# who needs arrays anyways
function pbc_dist(p1, p2, pbc, inv_pbc)
    # pbc (and inv_pbc by extension) have to be upper-triagonal now!
    # direct distance between the two points
        # d = p1 .- p300
    d1, d2, d3 = p1[1] - p2[1], p1[2] - p2[2], p1[3] - p2[3]
    # convert distance to relative coordinates and save rounded values
        # m = round.(inv_pbc * d)
    m1 = round(inv_pbc[1, 1] * d1 + inv_pbc[1, 2] * d2 + inv_pbc[1, 3] * d3)
    m2 = round(inv_pbc[2, 2] * d2 + inv_pbc[2, 3] * d3)
    m3 = round(inv_pbc[3, 3] * d3)
    # convert rounded relative distance into real distance and subtract it from the real distance; return norm
        # d -= pbc * m
        # n = norm(d)
    n = sqrt(
        (d1 - pbc[1, 1] * m1 - pbc[1, 2] * m2 - pbc[1, 3] * m3)^2 + 
        (d2 - pbc[2, 2] * m2 - pbc[2, 3] * m3)^2 + 
        (d3 - pbc[3, 3] * m3)^2
        )
    return n
end

"""
    minimum_image_vector(point1, point2, mdbox::MDBox)

Returns the vector connecting two points according to the minimum image convenction.

For distances above half of the shortest box diameter (see [`max_distance_for_pbc_dist(traj::Trajectory)`](@ref)), this might give erroneous results. 
"""
function minimum_image_vector(point1, point2, mdbox::TriclinicBox)
    minimum_image_vector(point1, point2, mdbox.pbc_matrix, mdbox.inv_pbc_matrix, mdbox.realspace_tmp, mdbox.inversespace_tmp)
end

function minimum_image_vector(point1, point2, pbc, inv_pbc = inv(pbc), dist_tmp = zeros(MVector{3}), matmul_tmp = zeros(MVector{3}))
    dist_tmp .= point1 .- point2
    mul!(matmul_tmp, inv_pbc,dist_tmp)
    matmul_tmp .-= floor.(matmul_tmp .+ 0.5)
    mul!(dist_tmp, pbc, matmul_tmp)
    return dist_tmp
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
    index, distance = 1, pbc_dist(point1, point2, mdbox)
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
