using ..MovingBoundaryProblems1D
using CairoMakie
using OrdinaryDiffEq
using ReferenceTests
using LinearSolve
using StatsBase
using SteadyStateDiffEq
using DataInterpolations
MB = MovingBoundaryProblems1D

@static if VERSION < v"1.9"
    stack(x) = reduce(hcat, x)
end

# Fisher-Stefan model: See e.g. https://doi.org/10.1137/090771089
#    ∂u/∂t = ∂²u/∂x² + u(1-u)      0 < x < L(t)
#    ∂u/∂x = 0                            x = 0
#        u = 0                         x = L(t)
#    dL/dt = -κ∂u(L(t), t)/∂x          x = L(t) 
#  u(x, 0) = x < β ? α : 0.0         0 ≤ x ≤ L₀
#     L(0) = β

function construct_problem(κ, α, β, T, n=500)
    mesh_points = LinRange(0, 1, n)
    diffusion_function = (u, x, t, p) -> one(u)
    reaction_function = (u, x, t, p) -> u * (one(u) - u)
    lhs = Neumann(0.0)
    rhs = Dirichlet(0.0)
    moving_boundary = Robin(0.0, -κ)
    ic = x -> x < β ? α : 0.0
    initial_condition = ic.(mesh_points)
    prob = MBProblem(
        mesh_points, lhs, rhs, moving_boundary;
        diffusion_function,
        reaction_function,
        initial_condition,
        initial_endpoint=β,
        final_time=T
    )
    return prob, mesh_points
end
function solve_problem(κ, α, β, T, n=500)
    prob, mesh_points = construct_problem(κ, α, β, T, n)
    sol = solve(prob, TRBDF2(linsolve=KLUFactorization()), saveat=T / 250)
    scaled_mesh = scaled_mesh_points(prob, sol)
    u = sol[begin:(end-1), :]
    L = sol[end, :]
    t = sol.t

    sprob = SteadyMBProblem(prob)
    ssol = solve(sprob, DynamicSS(TRBDF2()))
    su = ssol.u[begin:(end-1)]
    sL = ssol.u[end]

    return prob, stack(scaled_mesh), u, t, L, su, sL, mesh_points .* sL, maximum(L), sol
end

prob1, x1, u1, t1, L1, su1, sL1, sx1, mL1, msol1 = solve_problem(20.0, 0.5, 1.0, 20.0);
prob2, x2, u2, t2, L2, su2, sL2, sx2, mL2, msol2 = solve_problem(0.45, 0.5, 1.0, 20.0);
@test u1[end, :] ≈ zero(u1[end, :]) atol = 1e-4
@test u2[end, :] ≈ zero(u2[end, :]) atol = 1e-4

fig = Figure(fontsize=33)
idx_rng = (1, 50, 100, 150, 200, 250)
colors = (:red, :blue, :black, :magenta, :darkgreen, :orange)
ax1 = Axis(fig[1, 1], width=600, height=300, xlabel=L"x", ylabel=L"u(x, t)", title=L"(a): $\kappa = 20$", titlealign=:left)
ax2 = Axis(fig[1, 2], width=600, height=300, xlabel=L"x", ylabel=L"u(x, t)", title=L"(b): $\kappa = 0.45$", titlealign=:left)
ax3 = Axis(fig[2, 1], width=600, height=300, xlabel=L"x", ylabel=L"t", title=L"(c): $\kappa = 20$", titlealign=:left)
ax4 = Axis(fig[2, 2], width=600, height=300, xlabel=L"x", ylabel=L"t", title=L"(d): $\kappa = 20$", titlealign=:left)
[lines!(ax1, x1[:, idx], u1[:, idx], color=clr, linewidth=1.5) for (idx, clr) in zip(idx_rng, colors)]
[lines!(ax2, x2[:, idx], u2[:, idx], color=clr, linewidth=1.5) for (idx, clr) in zip(idx_rng, colors)]
tricontourf!(ax3, vec(x1), repeat(t1, inner=size(x1, 1)), vec(u1), levels=100)
tricontourf!(ax4, vec(x2), repeat(t2, inner=size(x2, 1)), vec(u2), levels=20)
resize_to_layout!(fig)
fig
fig_path = normpath(@__DIR__, "..", "docs", "src", "figures")
@test_reference joinpath(fig_path, "fisher_stefan_1.png") fig

using StatsBase
function estimate_wave_speed(sol)
    Ls = @views sol[end, (end÷2):end] # only take late time 
    ts = @views sol.t[(end÷2):end]
    cs = diff(Ls) ./ diff(ts)
    return mean(cs)
end
c = estimate_wave_speed(msol1)
@test c ≈ 1.2 rtol = 1e-1

@test su2 ≈ zero(su2) atol = 1e-6
@test sL2 ≈ 1.45465 rtol = 1e-2

#=
function compute_wave_speed(κ, α, β, T)
    prob, x, u, t, L, su, sL, sx, mL, sol = solve_problem(κ, α, β, T, 2_000)
    if sL > 200.0 # the steady state solution makes no sense in this case, so we have a travelling wave
        c = estimate_wave_speed(sol)
        return c
    else
        return NaN
    end
end

κ = LinRange(0.4, 0.5, 20)
c = compute_wave_speed.(κ, 0.5, 1.0, 5.0)

fig = Figure(fontsize=33)
ax = Axis(fig[1, 1], xlabel=L"\kappa", ylabel=L"c", width=600, height=300)
lines!(ax, κ, c, linewidth=3, color=:red)
resize_to_layout!(fig)
xlims!(ax, 0, 2)
ylims!(ax, 0, 1/2)
@test_reference joinpath(fig_path, "fisher_stefan_2.png") fig

κc = κ[findfirst(!isnan, c)]
@test κc ≈ 0.48 rtol=1e-2
=#