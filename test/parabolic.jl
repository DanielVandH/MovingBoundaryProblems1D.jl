using ..MovingBoundaryProblems1D
using CairoMakie
using OrdinaryDiffEq
using LinearSolve
using ReferenceTests
using SteadyStateDiffEq
MB = MovingBoundaryProblems1D

@static if VERSION < v"1.9"
    stack(x) = reduce(hcat, x)
end

mesh_points = LinRange(0, 1, 1_000)
diffusion_function = (u, x, t, p) -> one(u)
reaction_function = (u, x, t, p) -> u * (one(u) - u)
lhs = Dirichlet(0.0)
rhs = Dirichlet(0.0)
moving_boundary = Robin(0.0, -1/2)
ic = x -> 2x * (1 - x)
initial_condition = ic.(mesh_points)
initial_endpoint = 1.0
final_time = 0.4
prob = MBProblem(mesh_points, lhs, rhs, moving_boundary;
    diffusion_function,
    reaction_function,
    initial_condition,
    initial_endpoint,
    final_time)
sol = solve(prob, TRBDF2(linsolve=KLUFactorization()))
@test sol[end, :][end] ≈ solve(SteadyMBProblem(prob), DynamicSS(TRBDF2())).u[end] rtol=1e-2

let x = (stack ∘ scaled_mesh_points)(prob, sol), t = repeat(sol.t, inner=length(mesh_points)), u = sol[begin:(end-1), :]
    fig = Figure(fontsize=33)
    ax = Axis(fig[1, 1], width=600, height=300, xlabel=L"x", ylabel=L"t")
    tricontourf!(ax, vec(x), t, vec(u), levels=0.0:0.05:0.5, extendlow=:auto)
    resize_to_layout!(fig)
    fig_path = normpath(@__DIR__, "..", "docs", "src", "figures")
    @test_reference joinpath(fig_path, "parabolic.png") fig
    fig
end