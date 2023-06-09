using ..MovingBoundaryProblems1D
using CairoMakie
using OrdinaryDiffEq
using LinearSolve
using ReferenceTests
using DataInterpolations
using SpecialFunctions
MB = MovingBoundaryProblems1D

@static if VERSION < v"1.9"
    stack(x) = reduce(hcat, x)
end

## Define the parameters 
k, s, η, β = 10.0, 1.0, 1.0, 0.00577
F = (q, p) -> p.k * (p.s - q)
D = (q, p) -> p.k / (p.η * q^2)
G = (q, p) -> p.β
L₀ = 10.0
N₀ = 40.0
σ = sqrt(3)

## Define the initial condition 
mesh_points = LinRange(0, L₀, 1000)
q₀ = x -> N₀ * exp(-(2x - L₀)^2 / (8σ^2)) / (erf(L₀ * sqrt(2) / (4σ)) * sqrt(2π * σ^2))
initial_condition = q₀.(mesh_points)

## Define the PDE 
diffusion_function = (q, x, t, p) -> p.D(q, p)
diffusion_parameters = (D=D, k=k, η=η)
reaction_function = (q, x, t, p) -> q * p.G(inv(q), p)
reaction_parameters = (G=G, β=β)

## Define the boundary conditions 
lhs = Neumann(0.0)
rhs_f = (q, t, p) -> -2q * p.F(inv(q), p) / (p.η * p.D(q, p))
rhs_p = (F=F, η=η, D=D, s=s, k=k)
rhs = Neumann(rhs_f, rhs_p)
moving_boundary_f = (q, t, p) -> (zero(q), -p.D(q, p) / q)
moving_boundary_p = (D=D, k=k, η=η)
moving_boundary = Robin(moving_boundary_f, moving_boundary_p)

## Define the problem 
prob = MBProblem(mesh_points, lhs, rhs, moving_boundary;
    diffusion_function,
    diffusion_parameters,
    reaction_function,
    reaction_parameters,
    initial_condition,
    initial_endpoint=L₀,
    final_time=400.0)

## Solve the problem 
sol = solve(prob, TRBDF2(linsolve=KLUFactorization()))

## Compute the cell numbers 
function integrate_solution(prob, sol)
    N = zeros(length(sol))
    x = prob.geometry.mesh_points
    for i in eachindex(sol)
        q = @views sol.u[i][begin:(end-1)]
        L = sol.u[i][end]
        interp = LinearInterpolation(q, x)
        N[i] = DataInterpolations.integral(interp, 0.0, 1.0) * L # ∫₀ᴸ q(x) dx = L∫₀¹ q(ξ) dξ 
    end
    return N
end
Nt = integrate_solution(prob, sol)

## Plot 
fig = Figure(fontsize=33)
colors = (:red, :blue, :black, :magenta, :darkgreen)

ax1 = Axis(fig[1, 1], xlabel=L"x", ylabel=L"q(x, t)",
    title=L"(a): $q(x, t)$", titlealign=:left,
    width=600, height=300)
ax2 = Axis(fig[1, 2], xlabel=L"t", ylabel=L"N(t)",
    title=L"(b): $N(t)$", titlealign=:left,
    width=600, height=300)
ax3 = Axis(fig[1, 3], xlabel=L"t", ylabel=L"L(t)",
    title=L"(c): $L(t)$", titlealign=:left,
    width=600, height=300)

t = [0.0, 100.0, 200.0, 300.0, 400.0]
qL = sol.(t)
q = [@views qL[begin:(end-1)] for qL in qL]
L = [qL[end] for qL in qL]
ξ_grid = prob.geometry.mesh_points
[lines!(ax1, ξ_grid .* L[i], q[i], color=colors[i]) for i in eachindex(q)]
lines!(ax2, sol.t, Nt, color=:black, linewidth=3)
lines!(ax3, sol.t, @views(sol[end, :]), color=:black, linewidth=3)
resize_to_layout!(fig)
fig_path = normpath(@__DIR__, "..", "docs", "src", "figures")
@test_reference joinpath(fig_path, "epithelial.png") fig