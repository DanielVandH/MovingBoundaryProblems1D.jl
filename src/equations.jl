function mb_odes!(duLdt, uL, prob::P, t) where {P}
    u = @views uL[begin:(end-1)]
    L = uL[end]
    dudt = @views duLdt[begin:(end-1)]

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
        hₙ₋₁ = h[end]
        hₙ₋₂ = h[end-1]
        uₙ = u[end]
        uₙ₋₁ = u[end-1]
        uₙ₋₂ = u[end-2]
        ∂q = hₙ₋₁ / (hₙ₋₂ * (hₙ₋₂ + hₙ₋₁)) * uₙ₋₂ - (hₙ₋₂ + hₙ₋₁) / (hₙ₋₂ * hₙ₋₁) * uₙ₋₁ + (hₙ₋₂ + 2hₙ₋₁) / (hₙ₋₂ * (hₙ₋₁ + hₙ₋₂)) * uₙ
        dLdt = a₂ + (b₂ / L) * ∂q
    end
    duLdt[end] = dLdt

    # Now compute the dudt at the interior nodes
    for i in (firstindex(dudt)+1):(lastindex(dudt)-1)
        Vᵢ = V[i]
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
        term1 = 1 / (Vᵢ * L) * dLdt * (eᵢ * (uᵢ + uᵢ₊₁) / 2 - wᵢ * (uᵢ₋₁ + uᵢ) / 2)
        term2 = -(1 / L) * dLdt * uᵢ
        term3 = Rᵢ
        term4 = 1 / (Vᵢ * L^2) * (((Dᵢ + Dᵢ₊₁) / 2) * ((uᵢ₊₁ - uᵢ) / hᵢ) - ((Dᵢ₋₁ + Dᵢ) / 2) * ((uᵢ - uᵢ₋₁) / hᵢ₋₁))
        dudt[i] = term1 + term2 + term3 + term4
    end

    # Get the LHS dudt 
    if is_dirichlet(lhs)
        dudt[begin] = zero(dudt[begin])
    else
        V₁ = V[begin]
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
        term1 = 1 / (V₁ * L) * dLdt * e₁ * ((u₁ + u₂) / 2)
        term2 = -(1 / L) * dLdt * u₁
        term3 = R₁
        term4 = 1 / (V₁ * L^2) * (((D₁ + D₂) / 2) * ((u₂ - u₁) / h₁) - D₁ * L * a₀)
        dudt[begin] = term1 + term2 + term3 + term4
    end

    # Get the RHS dudt
    if is_dirichlet(rhs)
        dudt[end] = zero(dudt[end])
    else
        Vₙ = V[end]
        ξₙ₋₁ = mesh_points[end-1]
        ξₙ = mesh_points[end]
        xₙ₋₁ = ξₙ₋₁ * L
        xₙ = ξₙ * L
        hₙ₋₁ = h[end]
        uₙ₋₁ = u[end-1]
        uₙ = u[end]
        wₙ = (ξₙ₋₁ + ξₙ) / 2
        a₁ = rhs(uₙ, t)
        Dₙ₋₁ = D(uₙ₋₁, xₙ₋₁, t, Dp)
        Dₙ = D(uₙ, xₙ, t, Dp)
        Rₙ = R(uₙ, xₙ, t, Rp)
        term1 = 1 / (Vₙ * L) * dLdt * (uₙ - wₙ * ((uₙ₋₁ + uₙ) / 2))
        term2 = -(1 / L) * dLdt * uₙ
        term3 = Rₙ
        term4 = 1 / (Vₙ * L^2) * (Dₙ * L * a₁ - ((Dₙ₋₁ + Dₙ) / 2) * ((uₙ - uₙ₋₁) / hₙ₋₁))
        dudt[end] = term1 + term2 + term3 + term4
    end
    return duLdt
end

#=
Similar to a FVMProblem, but we also add dependency on the last and second last boundary point 
(if Dirichlet is used at RHS, uₙ and uₙ₋₁ are used)  and on L since each equation has dLdt and L.
=#
function get_nnz_itr(i, n)
    if i == 1
        return (1, 2, n, n - 1, n + 1)
    elseif 1 < i < n - 2
        return (i - 1, i, i + 1, n, n - 1, n + 1)
    elseif i == n - 2
        return (n - 3, n - 2, n - 1, n, n + 1)
    elseif i == n - 1
        return (n - 2, n - 1, n, n + 1)
    elseif i == n
        return (n, n - 1, n + 1)
    else # i == n + 1 
        return (n, n - 1, n + 1)
    end
end
jacobian_sparsity(prob::MBProblem) = jacobian_sparsity(prob.geometry.mesh_points)
function jacobian_sparsity(pts)
    n = length(pts)
    num_nnz = 6n - 4
    I = zeros(Int64, num_nnz)
    J = zeros(Int64, num_nnz)
    V = ones(num_nnz)
    ctr = 1
    for i in 1:(n+1)
        itr = get_nnz_itr(i, n)
        for j in itr
            I[ctr] = i
            J[ctr] = j
            ctr += 1
        end
    end
    return sparse(I, J, V)
end