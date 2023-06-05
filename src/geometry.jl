"""
    MBGeometry{T}

Definition of the geometry for a moving boundary problem problem.

# Fields
- `mesh_points::T`: The mesh points. Must be sorted and satisfy `mesh_points[begin] = 0` and `mesh_points[end] = 1` (if they do not satisfy this, they are scaled first).
- `spacings::T`: The spacings between the mesh points. 
- `volumes::T`: The volumes of the cells defined by the mesh points.

# Constructors

To construct the geometry, you can directly call the default constructor, 

    MBGeometry(mesh_points, spacings, volumes)

or you can call the convenience constructor,

    MBGeometry(mesh_points)

which will compute the spacings and volumes for you.

See also [`MBProblem`](@ref).
"""
struct MBGeometry{T}
    mesh_points::T
    spacings::T
    volumes::T
    function MBGeometry(mesh_points, spacings, volumes)
        @assert issorted(mesh_points) "mesh_points is not sorted."
        @assert mesh_points[begin] == 0 && mesh_points[end] == 1 "mesh_points must satisfy mesh_points[begin] = 0 and mesh_points[end] = 1."
        @assert length(mesh_points) == length(volumes) == length(spacings) + 1 "mesh_points and volumes must have the same length, and spacings must have one less element."
        mesh_points, spacings, volumes = promote(collect(mesh_points), spacings, volumes)
        T = typeof(mesh_points)
        return new{T}(mesh_points, spacings, volumes)
    end
end
function MBGeometry(mesh_points)
    if !(mesh_points[begin] == 0 && mesh_points[end] == 1)
        a = mesh_points[begin]
        b = mesh_points[end]
        mesh_points = (mesh_points .- a) ./ (b - a)
    end
    spacings = compute_spacings(mesh_points)
    volumes = compute_volumes(mesh_points, spacings)
    return MBGeometry(mesh_points, spacings, volumes)
end
function compute_spacings(mesh_points)
    return diff(mesh_points)
end
function compute_volumes(mesh_points, spacings)
    V = similar(mesh_points)
    V[begin] = 0.5spacings[begin]
    V[end] = 0.5spacings[end]
    for i in (firstindex(mesh_points)+1):(lastindex(mesh_points)-1)
        V[i] = 0.5 * (spacings[i-1] + spacings[i])
    end
    return V
end