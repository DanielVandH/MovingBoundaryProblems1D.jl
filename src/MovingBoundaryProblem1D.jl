module MovingBoundaryProblems1D

using SparseArrays
using SciMLBase
using CommonSolve

export FVMGeometry, BoundaryConditions, FVMProblem
export Dirichlet, Neumann, Robin
export solve

include("geometry.jl")
include("boundary_conditions.jl")
include("problem.jl")
include("equations.jl")
include("solve.jl")

@static if VERSION < v"1.7"
    _stable_typeof(x) = typeof(x)
    _stable_typeof(::Type{T}) where {T} = @isdefined(T) ? Type{T} : DataType
    struct Returns{V} <: Function
        value::V
        Returns{V}(value) where {V} = new{V}(value)
        Returns(value) = new{_stable_typeof(value)}(value)
    end
    (obj::Returns)(@nospecialize(args...); @nospecialize(kw...)) = obj.value
end

end # module 