function mb_odes!(duLdt, uL, prob::P, t) where {P}
    u = @views uL[begin:(end-1)]
    L = uL[end]
    dudt = @views duLdt[begin:(end-1)]
    n = length(u)

    mesh = prob.geometry
    mesh_points = mesh.mesh_points
    V = mesh.volumes
    h = mesh.spacings

    boundary_conditions = prob.boundary_conditions
    lhs = boundary_conditions.lhs
    rhs = boundary_conditions.rhs
    moving_boundary = boundary_conditions.moving_boundary

    D = prob.diffusion_function
    Dp = prob.diffusion_parameters
    R = prob.reaction_function
    Rp = prob.reaction_parameters

    # Start by computing dLdt
    a₂, b₂ = moving_boundary(u[end], t)
    if !is_dirichlet(rhs)
        a₁ = rhs(u[end], t)
        dLdt = a₂ + a₁ * b₂
    else
        hₙ₋₁ = h[end-1]
        uₙ₋₁ = u[end-1]
        uₙ = u[end]
        dLdt = a₂ + (b₂ / L) * (uₙ - uₙ₋₁) / hₙ₋₁
    end
    duLdt[end] = dLdt

    # Now compute the dudt at the interior nodes
    for i in (firstindex(dudt)+1):(lastindex(dudt)-1)
        Vᵢ = volumes[i]
        ξᵢ₋₁ = mesh_points[i-1]
        ξᵢ = mesh_points[i]
        ξᵢ₊₁ = mesh_points[i+1]
        xᵢ₋₁ = ξᵢ₋₁ * L
        xᵢ = ξᵢ * L
        xᵢ₊₁ = ξᵢ₊₁ * L
        hᵢ₋₁ = h[i-1]
        hᵢ = h[i]
        wᵢ = (ξᵢ₋₁ + ξᵢ) / 2
        eᵢ = (ξᵢ + ξᵢ₊₁) / 2
        uᵢ₋₁ = u[i-1]
        uᵢ = u[i]
        uᵢ₊₁ = u[i+1]
        Dᵢ₋₁ = D(uᵢ₋₁, xᵢ₋₁, t, Dp)
        Dᵢ = D(uᵢ, xᵢ, t, Dp)
        Dᵢ₊₁ = D(uᵢ₊₁, xᵢ₊₁, t, Dp)
        Rᵢ = R(uᵢ, xᵢ, t, Rp)
        term1 = 1 / (Vᵢ * L) * dLdt * (eᵢ * (qᵢ + qᵢ₊₁) / 2 - wᵢ * (qᵢ₋₁ + qᵢ) / 2)
        term2 = -(1 / L) * dLdt * uᵢ
        term3 = Rᵢ
        term4 = 1 / (Vᵢ * L^2) * (((Dᵢ + Dᵢ₊₁) / 2) * ((qᵢ₊₁ - qᵢ) / hᵢ) - ((Dᵢ₋₁ + Dᵢ) / 2) * ((qᵢ - qᵢ₋₁) / hᵢ₋₁))
        dudt[i] = term1 + term2 + term3 + term4
    end

    # Get the LHS dudt 
    if is_dirichlet(lhs)
        dudt[begin] = lhs(u[begin], t)
    else
        V₁ = volumes[begin]
        ξ₁ = mesh_points[begin]
        ξ₂ = mesh_points[begin+1]
        x₁ = ξ₁ * L
        x₂ = ξ₂ * L
        h₁ = h[begin]
        u₁ = u[begin]
        u₂ = u[begin+1]
        e₁ = (ξ₁ + ξ₂) / 2
        a₀ = lhs(u₁, t)
        D₁ = D(u₁, x₁, t, Dp)
        D₂ = D(u₂, x₂, t, Dp)
        R₁ = R(u₁, x₁, t, Rp)
        term1 = 1 / (V₁ * L) * dLdt * e₁ * ((q₁ + q₂) / 2)
        term2 = -(1 / L) * dLdt * u₁
        term3 = R₁
        term4 = 1 / (V₁ * L^2) * (((D₁ + D₂) / 2) * ((q₂ - q₁) / h₁) - D₁ * L * a₀)
        dudt[begin] = term1 + term2 + term3 + term4
    end

    # Get the RHS dudt
    if is_dirichlet(rhs)
        dudt[end] = rhs(u[end], t)
    else
        Vₙ = volumes[end]
        ξₙ₋₁ = mesh_points[end-1]
        ξₙ = mesh_points[end]
        xₙ₋₁ = ξₙ₋₁ * L
        xₙ = ξₙ * L
        hₙ₋₁ = h[end-1]
        uₙ₋₁ = u[end-1]
        uₙ = u[end]
        wₙ = (ξₙ₋₁ + ξₙ) / 2
        a₁ = rhs(uₙ, t)
        Dₙ₋₁ = D(uₙ₋₁, xₙ₋₁, t, Dp)
        Dₙ = D(uₙ, xₙ, t, Dp)
        Rₙ = R(uₙ, xₙ, t, Rp)
        term1 = 1 / (Vₙ * L) * dLdt * (uₙ - wₙ * ((qₙ₋₁ + qₙ) / 2))
        term2 = -(1 / L) * dLdt * uₙ
        term3 = Rₙ
        term4 = 1 / (Vₙ * L^2) * (Dₙ * L * a₁ - ((Dₙ₋₁ + Dₙ) / 2) * ((qₙ - qₙ₋₁) / hₙ₋₁))
        dudt[end] = term1 + term2 + term3 + term4
    end
    return nothing
end

# similar to a FVMProblem, but we also add dependency on the last boundary point and on L since each equation has dLdt and L
jacobian_sparsity(prob::MBProblem) = jacobian_sparsity(prob.geometry.mesh_points)
function jacobian_sparsity(pts)
    n = length(pts)
    num_nnz = 5n - 2
    I = zeros(Int64, num_nnz)
    J = zeros(Int64, num_nnz)
    V = ones(num_nnz)
    I[1] = 1
    J[1] = 1
    I[2] = 1
    J[2] = 2
    I[3] = 1
    J[3] = n
    I[4] = 1
    J[4] = n + 1
    ctr = 5
    for i in 2:(n-2)
        I[ctr] = i
        J[ctr] = i - 1
        ctr += 1
        I[ctr] = i
        J[ctr] = i
        ctr += 1
        I[ctr] = i
        J[ctr] = i + 1
        ctr += 1
        I[ctr] = i
        J[ctr] = n
        ctr += 1
        I[ctr] = i
        J[ctr] = n + 1
        ctr += 1
    end
    I[ctr] = n - 1
    J[ctr] = n - 2
    ctr += 1
    I[ctr] = n - 1
    J[ctr] = n - 1
    ctr += 1
    I[ctr] = n - 1
    J[ctr] = n
    ctr += 1
    I[ctr] = n - 1
    J[ctr] = n + 1
    ctr += 1
    I[ctr] = n
    J[ctr] = n - 1
    ctr += 1
    I[ctr] = n
    J[ctr] = n
    ctr += 1
    I[ctr] = n
    J[ctr] = n + 1
    ctr += 1
    I[ctr] = n + 1
    J[ctr] = n
    ctr += 1
    I[ctr] = n + 1
    J[ctr] = n + 1
    return sparse(I, J, V)
end