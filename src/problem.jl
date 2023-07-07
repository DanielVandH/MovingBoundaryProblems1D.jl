"""
    MBProblem{T,DF,DP,RF,RP,L,R,M,IC,IE,FT}

Definition of an `MBProblem`.

# Fields 
- `geometry::MBGeometry{T}`: The geometry of the problem.
- `boundary_conditions::BoundaryConditions{L, R, M}`: The boundary conditions.
- `diffusion_function::DF`: The diffusion function, of the form `(u, x, t, p) -> Number`, where `p = diffusion_parameters`.
- `diffusion_parameters::DP`: The parameters for the diffusion function.
- `reaction_function::RF`: The reaction function, of the form `(u, x, t, p) -> Number`, where `p = reaction_parameters`.
- `reaction_parameters::RP`: The parameters for the reaction function.
- `initial_condition::IC`: The initial condition, with `initial_condition[i]` corresponding to the value at `geometry.mesh_points[i]` and `t = initial_time`.
- `initial_endpoint::IE`: The initial endpoint. This is the value of the moving boundary at `t = initial_time`.
- `initial_time::FT`: The initial time.
- `final_time::FT`: The final time.

See also [`SteadyMBProblem`](@ref).

# Constructors

You can use the default constructor, but we also provide the constructor 

    MBProblem(;
        geometry, 
        boundary_conditions,
        diffusion_function,
        diffusion_parameters = nothing,
        reaction_function = Returns(0.0),
        reaction_parameters = nothing,
        initial_condition,
        initial_endpoint,
        initial_time = 0.0,
        final_time)

which provides some default values. Moreover, instead of providing `geometry` and `boundary_conditions`, you can use 

    MBProblem(mesh_points, lhs, rhs, moving_boundary; kwargs...)

which will construct `geometry = MBGeometry(mesh_points)` and `boundary_conditions = BoundaryConditions(lhs, rhs, moving_boundary)`. 
The `kwargs...` are as above, except without `geometry` and `boundary_conditions` of course.

To solve the `MBProblem`, just use `solve` as you would in DifferentialEquations.jl. For example, 

    sol = solve(prob, TRBDF2(linsolve=KLUFactorization()), saveat=0.1)

solves the problem using the `TRBDF2()` algorithm, with the associated linear systems solved with the sparse solver `KLUFactorization()`,
and the solution will be saved every `0.1` units of time from `initial_time` up to, and including, `final_time`. The solution in this case 
will be such that `sol.u[i]` has the values of `u` in `sol.u[i][begin:(end-1)]`, and the position of the moving boundary at `sol.t[i]` in 
`sol.u[i][end]`. See also [`scaled_mesh_points`](@ref) for the corresponding grid points.
"""
Base.@kwdef struct MBProblem{T,DF,DP,RF,RP,L,R,M,IC,IE,FT}
    geometry::MBGeometry{T}
    boundary_conditions::BoundaryConditions{L,R,M}
    diffusion_function::DF
    diffusion_parameters::DP = nothing
    reaction_function::RF = Returns(0.0)
    reaction_parameters::RP = nothing
    initial_condition::IC
    initial_endpoint::IE
    initial_time::FT = 0.0
    final_time::FT
    function MBProblem(
        gemoetry::MBGeometry{T}, boundary_conditions::BoundaryConditions{L,R,M},
        diffusion_function::DF, diffusion_parameters::DP,
        reaction_function::RF, reaction_parameters::RP,
        initial_condition::IC, initial_endpoint::IE,
        initial_time::FT, final_time::FT) where {T,DF,DP,RF,RP,L,R,M,IC,IE,FT}
        @assert initial_endpoint > 0.0 "initial_endpoint must be positive."
        return new{T,DF,DP,RF,RP,L,R,M,IC,IE,FT}(
            gemoetry, boundary_conditions,
            diffusion_function, diffusion_parameters,
            reaction_function, reaction_parameters,
            initial_condition, initial_endpoint,
            initial_time, final_time)
    end
end
MBProblem(mesh_points, lhs, rhs, moving_boundary; kwargs...) = MBProblem(;
    geometry=MBGeometry(mesh_points),
    boundary_conditions=BoundaryConditions(lhs, rhs, moving_boundary),
    kwargs...)

"""
    SteadyMBProblem{M}

Defines a steady-state problem for a moving boundary problem. Only has a 
single field, `prob::MBProblem`.

You can `solve` this problem as you would a `NonlinearProblem`, e.g. with 

    solve(prob, alg)

where `alg` is e.g. `DynamicSS(TRBDF2()))` from SteadyStateDiffEq.jl.
"""
struct SteadyMBProblem{M}
    prob::M
end