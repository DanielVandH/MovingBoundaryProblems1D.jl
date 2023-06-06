using ..MovingBoundaryProblems1D
using CairoMakie
using OrdinaryDiffEq
using LinearSolve
using ReferenceTests
MB = MovingBoundaryProblems1D

@static if VERSION < v"1.9"
    stack(x) = reduce(hcat, x)
end

# Heat equation with a Stefan BC: https://doi.org/10.1016/S0377-0427(97)00034-4
#   ∂u/∂t = ∂²u/∂x²         0 < x < L(t) 
#   ∂u/∂t = -exp(t)                x = 0
#       u = 0                   x = L(t) 
#   dL/dt = -∂u/∂x              x = L(t)
#    L(0) = 0                 
# u(x, 0) = 0               0 ≤ x ≤ L(0) = 0
#
# The exact exact solution is 
#   u(x, t) = exp(t - x) - 1, 0 ≤ x ≤ L(t), 0 < t < 1
#      L(t) = t

# Define the exact solutions for gettig the initial data
exact_u = (x, t) -> exp(t - x) - 1
exact_L = t -> t

# Define the boundary conditions
lhs = Neumann((u, t, p) -> -exp(t))
rhs = Dirichlet(0.0)
moving_boundary = Robin(0.0, -1.0)

# Setup the initial data and the PDE
initial_time = 0.1
final_time = 0.5
initial_endpoint = exact_L(initial_time)
mesh_points = LinRange(0, initial_endpoint, 500)
initial_condition = exact_u.(mesh_points, initial_time)
diffusion_function = (u, x, t, p) -> one(u)

# Define the problem
prob = MBProblem(mesh_points, lhs, rhs, moving_boundary;
    diffusion_function,
    initial_time,
    final_time,
    initial_endpoint,
    initial_condition)

sol = solve(prob, TRBDF2(linsolve=KLUFactorization()), saveat=0.1)

scaled_mesh = scaled_mesh_points(prob, sol)
@test reduce(hcat, scaled_mesh) ≈ prob.geometry.mesh_points .* [sol.u[i][end] for i in eachindex(sol)]'

exact_solu = [vcat(exact_u.(x, t), exact_L(t)) for (x, t) in zip(scaled_mesh, sol.t)]
@test sol.u ≈ exact_solu rtol = 1e-4
@test sol[end-1, :] ≈ zeros(length(sol)) atol = 1e-4

fig = Figure(fontsize=35)
ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"u(x, t)", width=600, height=300)
colors = [:red, :blue, :black, :magenta, :darkgreen]
[lines!(ax, x, u[begin:(end-1)], color=clr, linewidth=1.5) for (x, u, clr) in zip(scaled_mesh, sol.u, colors)]
[lines!(ax, x, exact_u.(x, t), color=clr, linestyle=:dash, linewidth=3) for (x, t, clr) in zip(scaled_mesh, sol.t, colors)]
resize_to_layout!(fig)
fig_path = normpath(@__DIR__, "..", "docs", "src", "figures")
@test_reference joinpath(fig_path, "heat_equation.png") fig
