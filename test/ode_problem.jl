using ..MovingBoundaryProblems1D
using SciMLBase
using SparseArrays
MB = MovingBoundaryProblems1D

mesh_points = sort(rand(100))
geo = MBGeometry(mesh_points)
diffusion_function = (u, t, p) -> u * p[1] + 2.0
diffusion_parameters = (1.0,)
reaction_function = (u, t, p) -> u * p[1] + 2.0 + p[2]^2
reaction_parameters = (1.0, 2.0)
initial_condition = rand(100)
initial_endpoint = 2.8
initial_time = 0.3
final_time = 2.8771
lhs = Dirichlet(0.5)
rhs = Neumann((u, t, p) -> 0.39u + p, 0.5)
moving_boundary = Robin(0.5, 0.3)
prob = MBProblem(;
    geometry=geo,
    boundary_conditions=BoundaryConditions(lhs, rhs, moving_boundary),
    diffusion_function=diffusion_function,
    diffusion_parameters=diffusion_parameters,
    reaction_function=reaction_function,
    reaction_parameters=reaction_parameters,
    initial_condition=initial_condition,
    initial_endpoint=initial_endpoint,
    final_time=final_time,
    initial_time=initial_time)
J = MB.jacobian_sparsity(prob)
@test J == MB.jacobian_sparsity(mesh_points)
_J = zeros(length(mesh_points) + 1, length(mesh_points) + 1)
n = length(mesh_points)
for i in 1:(n+1)
    if i == 1
        idx = [1, 2, n - 1, n, n + 1]
    elseif i == n
        idx = [n, n - 1, n + 1]
    elseif i == n + 1
        idx = [n, n - 1, n + 1]
    else
        idx = unique([i - 1, i, i + 1, n, n - 1, n + 1])
    end
    _J[i, idx] .= 1
end
@test J == _J
@test nnz(J) == sum(_J) == 6n - 4
ode_prob = ODEProblem(prob)
@test ode_prob.p == prob
@test ode_prob.f.f.f == MB.mb_odes!
@test ode_prob.u0 == [initial_condition; initial_endpoint]
@test ode_prob.tspan == (prob.initial_time, final_time)
@test ode_prob.f.jac_prototype == J
cb = ode_prob.kwargs.data.callback.discrete_callbacks
lhs_cb, rhs_cb = cb
@test lhs_cb.condition(0.0, 0.0, 0.0)
@test !rhs_cb.condition(0.0, 0.0, 0.0)
@test lhs_cb.affect! == MB.update_lhs!
@test rhs_cb.affect! == MB.update_rhs!