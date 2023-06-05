using ..MovingBoundaryProblems1D
MB = MovingBoundaryProblems1D

mesh_points = rand(100)
@test_throws AssertionError MBGeometry(mesh_points)
mesh_points = sort(mesh_points)
scaled_mesh_points = (mesh_points .- mesh_points[1]) ./ (mesh_points[end] - mesh_points[1])
geo = MBGeometry(mesh_points)
@test geo.mesh_points == scaled_mesh_points 
@test geo.spacings == diff(scaled_mesh_points) == MB.compute_spacings(scaled_mesh_points)
@test geo.volumes == MB.compute_volumes(scaled_mesh_points, geo.spacings)
@test geo.volumes == 0.5 * [geo.spacings[1]; geo.spacings[1:end-1] .+ geo.spacings[2:end]; geo.spacings[end]]
@test_throws AssertionError MBGeometry(mesh_points, geo.spacings, geo.volumes)
