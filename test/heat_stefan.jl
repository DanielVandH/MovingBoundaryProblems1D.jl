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

# See e.g. https://doi.org/10.1098/rspa.2019.0378 with λ = 0 

mesh_points = LinRange(0, 1, 500)

D = 0.1
diffusion_function = (u, x, t, p) -> p
diffusion_parameters = D
lhs = Neumann(0.0)
rhs = Dirichlet(0.0)
κ = 0.1
moving_boundary = Robin(0.0, -κ)

α = 0.5
β = 1.0
ic = x -> x < β ? α : 0.0
initial_condition = ic.(mesh_points)
initial_endpoint = β
final_time = 20.0

prob = MBProblem(mesh_points, lhs, rhs, moving_boundary;
    diffusion_function,
    diffusion_parameters,
    initial_condition,
    initial_endpoint=initial_endpoint,
    final_time)
sol = solve(prob, TRBDF2(linsolve=KLUFactorization()), saveat=final_time / 250)

let x = (stack ∘ scaled_mesh_points)(prob, sol), t = repeat(sol.t, inner=length(mesh_points)), u = sol[begin:(end-1), :]
    fig = Figure(fontsize=33)
    ax = Axis(fig[1, 1], width=600, height=300, xlabel=L"x", ylabel=L"t")
    tricontourf!(ax, vec(x), t, vec(u))
    resize_to_layout!(fig)
    fig_path = normpath(@__DIR__, "..", "docs", "src", "figures")
    @test_reference joinpath(fig_path, "heat_stefan.png") fig
    fig
end

sprob = SteadyMBProblem(prob)
ssol = solve(sprob, DynamicSS(TRBDF2()))
@test ssol.u[begin:(end-1)] ≈ zeros(length(mesh_points)) atol=1e-6
@test ssol.u[end] ≈ 3/2 rtol=1e-3