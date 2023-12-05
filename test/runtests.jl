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
    println("Starting the Heat equation example.")
    @safetestset "Heat equation" begin
        include("heat_equation.jl")
    end
    println("Starting the Stefan problem example.")
    @safetestset "Stefan problem" begin
        include("heat_stefan.jl")
    end
    println("Starting the Fisher-Stefan problem example.")
    @safetestset "Fisher-Stefan" begin
        include("fisher_stefan.jl")
    end
    println("Starting the Epithelial problem example.")
    @safetestset "Epithelial" begin
        include("epithelial.jl")
    end
end