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
           old_dist = point1 .- point2
           offsets = round.( old_dist .* mdbox.inv_cell_parameters )
           dist = old_dist .- offsets .* mdbox.cell_parameters
           return norm(dist)
       end


function pbc_dist(p1, p2, mdbox::TriclinicBox)
    pbc_dist(p1, p2, mdbox.pbc_matrix, mdbox.inv_pbc_matrix)
end

# ugliest implementation yet, but also fast
# who needs arrays anyways
function pbc_dist(p1, p2, pbc, inv_pbc)
    # pbc (and inv_pbc by extension) have to be upper-triagonal now!
    # direct distance between the two points
        # d = p1 .- p2
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
    minimum_image_vector(point1, point2, mdbox.pbc_matrix, mdbox.inv_pbc_matrix)
end

function minimum_image_vector(point1, point2, pbc, inv_pbc = inv(pbc))
    old_dist = point1 .- point2
    rel_dist = inv_pbc * old_dist
    rel_dist -= round.(rel_dist)
    new_dist = pbc * rel_dist
    return new_dist
end

function pbc_dist_triclinic(point1, point2, pbc, inv_pbc = inv(pbc))
    rel_dist = inv_pbc * (point2 - point1)
    rel_dist -= round.(rel_dist)
    direct_dist = pbc * rel_dist
    distance = Inf
    a, b, c = pbc[:, 1], pbc[:, 2], pbc[:, 3]
    for i in -1:1
        for j in -1:1
             for k in -1:1
                  test_dist = norm(direct_dist + i * a + j*b + k*c)
                  if test_dist < distance
                       distance = test_dist
                  end
             end
        end
    end
    return distance
end

function next_neighbor(point1, group2, mdbox::MDBox)
    point2 = group2[1]
    index, distance = 1, pbc_dist(point1, point2, mdbox)
    for i in 2:size(group2)[1]
        point2 = group2[i]
        new_distance = pbc_dist(point1, point2, mdbox)
        if new_distance < distance
            distance = new_distance
            index = i
        end
    end
    return index, distance
end

function next_neighbor(point1::SVector{3, Float64}, group2::Vector{SVector{3, Float64}}, mdbox::MDBox, atomtypes::Vector{String}, targettype::Vector{String})
    index, distance = -1, Inf
    for i in 1:length(group2)
        if !(atomtypes[i] in targettypes)
            continue
        end
        point2 = group2[i]
        new_distance = pbc_dist(point1, point2, mdbox)
        if new_distance < distance
            distance = new_distance
            index = i
        end
    end
    return index, distance
end
