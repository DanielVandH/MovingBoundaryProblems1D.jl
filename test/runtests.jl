using MovingBoundaryProblems1D
using Test
using SafeTestsets

@safetestset "MBGeometry" begin
    include("geometry.jl")
end
@safetestset "BoundaryConditions" begin
    include("boundary_conditions.jl")
end
@safetestset "MBProblem" begin
    include("problem.jl")
end
@safetestset "ODEProblem" begin
    include("ode_problem.jl")
end
