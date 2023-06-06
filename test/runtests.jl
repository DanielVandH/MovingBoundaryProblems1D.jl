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

@testset "Examples" begin
    @safetestset "Heat equation" begin
        include("heat_equation.jl")
    end
    @safetestset "Stefan problem" begin
        include("heat_stefan.jl")
    end
    @safetestset "Fisher-Stefan" begin
        include("fisher_stefan.jl")
    end
    @safetestset "Parabolic" begin
        include("parabolic.jl")
    end
end