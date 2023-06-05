abstract type AbstractBoundaryCondition{F,P} end
(bc::AbstractBoundaryCondition{F,P})(u, t) where {F,P} = bc.f(u, t, bc.p)

@doc raw"""
    Dirichlet{F,P} <: AbstractBoundaryCondition{F,P}

A Dirichlet boundary condition with fields `f` and `p` (default `p = nothing`),
with `f` being a function of the form `f(u, p)` and `p` being
the parameters for `f`. 

A Dirichlet boundary condition takes the form

```math
u(a, t) ↤ f(u(a, t), t, p),
```

where `a` is one of the endpoints. 

# Constructors 

    Dirichlet(f::Function, p = nothing) -> Dirichlet(f, p)
    Dirichlet(; f, p = nothing)         -> Dirichlet(f, p)
    Dirichlet(v::Number)                -> Dirichlet((u, t, p) -> oftype(u, v), nothing)
"""
Base.@kwdef struct Dirichlet{F,P} <: AbstractBoundaryCondition{F,P}
    f::F
    p::P = nothing
    Dirichlet(f::F, p::P=nothing) where {F,P} = new{F,P}(f, p)
end
Dirichlet(f::Function) = Dirichlet(f, nothing)
Dirichlet(v::Number) =
    let v = v
        Dirichlet((u, t, p) -> oftype(u, v))
    end

@doc raw"""
    Neumann{F,P} <: AbstractBoundaryCondition{F,P}

A Neumann boundary condition with fields `f` and `p` (default `p = nothing`),
with `f` being a function of the form `f(u, t, p)` and `p` being
the parameters for `f`.

A Neumann boundary condition takes the form

```math
\dfrac{\partial u}{\partial x}(a, t) = f(u(a, t), t, p),
```

where `a` is one of the endpoints. 

# Constructors 

    Neumann(f::Function, p = nothing) -> Neumann(f, p)
    Neumann(; f, p = nothing)         -> Neumann(f, p)
    Neumann(v::Number)                -> Neumann((u, t, p) -> oftype(u, v), nothing)
"""
Base.@kwdef struct Neumann{F,P} <: AbstractBoundaryCondition{F,P}
    f::F
    p::P = nothing
    Neumann(f::F, p::P=nothing) where {F,P} = new{F,P}(f, p)
end
Neumann(v::Number) =
    let v = v
        Neumann((u, t, p) -> oftype(u, v))
    end

@doc raw"""
    Robin{F,P} <: AbstractBoundaryCondition{F,P}

A Robin boundary condition with fields `f` and `p` (default `p = nothing`),
with `f` being a function of the form `f(u, t, p)` that returns a `Tuple` `(a₂, b₂)` 
and `p` being the parameters for `f`.

A Robin boundary condition takes the form

```math
dL/dt = a₂(u, t, p) + b₂(u, t, p)∂u/∂x,
```

where `f(u, t, p) = (a₂(u, t, p), b₂(u, t, p))`. 

# Constructors 

    Robin(f::Function, p = nothing) -> Robin(f, p)
    Robin(; f, p = nothing)         -> Robin(f, p)
    Robin(a::Number, b::Number)     -> Robin((u, t, p) -> (oftype(u, a), oftype(u, b)))
    Robin(a₂::Function, b₂::Function, p = nothing) -> Robin((u, t, p) -> (a₂(u, t, p), b₂(u, t, p)), p)
"""
Base.@kwdef struct Robin{F,P} <: AbstractBoundaryCondition{F,P}
    f::F
    p::P = nothing
    Robin(f::F, p::P=nothing) where {F,P} = new{F,P}(f, p)
end
Robin(a::Number, b::Number) =
    let a = a, b = b
        Robin((u, t, p) -> (oftype(u, a), oftype(u, b)))
    end
Robin(a₂::Function, b₂::Function, p) =
    let f = a₂, g = b₂, p = p
        Robin((u, t, p) -> (f(u, t, p), g(u, t, p)), p)
    end

is_dirichlet(::AbstractBoundaryCondition) = false
is_dirichlet(::Dirichlet) = true
is_neumann(::AbstractBoundaryCondition) = false
is_neumann(::Neumann) = true
is_robin(::AbstractBoundaryCondition) = false
is_robin(::Robin) = true

"""
    BoundaryConditions{L<:Union{<:Dirichlet, <:Neumann},R<:Union{<:Dirichlet, <:Neumann},M<:Robin}

The boundary conditions for the `MBProblem`.
    
# Fields 
- `lhs::L`: The left-hand side boundary condition. Must be `Dirichlet` or `Neumann`.
- `rhs::R`: The right-hand side boundary condition. Must be `Dirichlet` or `Neumann`.
- `moving_boundary::M`: The moving boundary condition. Must be `Robin`.

See also [`Dirichlet`](@ref), [`Neumann`](@ref), and [`Robin`](@ref) for the types of 
boundary conditions you can construct.
"""
Base.@kwdef struct BoundaryConditions{L,R,M}
    lhs::L
    rhs::R
    moving_boundary::M
    function BoundaryConditions(lhs::L, rhs::R, moving_boundary::M) where {L,R,M}
        if !(is_dirichlet(lhs) || is_neumann(lhs))
            throw(MethodError(BoundaryConditions, (lhs, rhs, moving_boundary)))
        elseif !(is_dirichlet(rhs) || is_neumann(rhs))
            throw(MethodError(BoundaryConditions, (lhs, rhs, moving_boundary)))
        elseif !is_robin(moving_boundary)
            throw(MethodError(BoundaryConditions, (lhs, rhs, moving_boundary)))
        end
        new{L,R,M}(lhs, rhs, moving_boundary)
    end
end