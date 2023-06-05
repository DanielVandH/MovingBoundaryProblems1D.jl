using ..MovingBoundaryProblems1D
MB = MovingBoundaryProblems1D

f = (u, t, p) -> 1.01u^p[1] + p[2] + t
p = (1.0, 2.0)
dirc = Dirichlet(f, p)
@test MB.is_dirichlet(dirc)
@test !MB.is_neumann(dirc)
@test !MB.is_robin(dirc)
@test dirc(0.881, 0.5) == f(0.881, 0.5, p)
@inferred dirc(0.881, 0.5) 
@test Dirichlet(f) == Dirichlet(f, nothing)
dirc = Dirichlet(0.5)
@test dirc(0.5, 0.2) == 0.5
@inferred dirc(0.5, 0.2)
dirc = Dirichlet(0)
@test dirc(0,0.2) == 0 && dirc(0,0.2) isa Int
@test dirc(0.0,0.2) == 0.0 && dirc(0.0,0.2) isa Float64
@inferred dirc(0,0.2)
@inferred dirc(0.0,0.2)

neum = Neumann(f, p)
@test !MB.is_dirichlet(neum)
@test MB.is_neumann(neum)
@test !MB.is_robin(neum)
@test neum(0.5, 0.25) == f(0.5,0.25, p)
@test neum(0.881, 0.192) == f(0.881, 0.192, p)
@inferred neum(0.881, 0.192)
@test Neumann(f) == Neumann(f, nothing)
neum = Neumann(0.5)
@test neum(0.5,0.2) == 0.5
@inferred neum(0.5,0.2)
neum = Neumann(0)
@test neum(0,0.2) == 0 && neum(0,0.2) isa Int
@test neum(0.0,0.2) == 0.0 && neum(0.0,0.2) isa Float64
@inferred neum(0,0.2)

robin = Robin(f, p)
@test !MB.is_dirichlet(robin)
@test !MB.is_neumann(robin)
@test MB.is_robin(robin)
robin = Robin(0.5, 0.2)
@test robin(0.3, 0.5) == (0.5, 0.2)
@inferred robin(0.3, 0.5)
f1 = (u, t, p) -> u*t + p[1]
f2 = (u, t, p) -> u^2 + t - p[2] 
robin = Robin(f1, f2, (1.3, 2.7))
@test robin(0.3, 0.5) == (0.3 * 0.5 + 1.3, 0.3^2 + 0.5 - 2.7)
robin = Robin(; f, p)
@test robin.f == f 
@test robin.p == p

bc = BoundaryConditions(dirc, neum, robin)
@test bc.lhs == dirc 
@test bc.rhs == neum
@test bc.moving_boundary == robin
@test_throws MethodError BoundaryConditions(robin, neum, robin)
@test_throws MethodError BoundaryConditions(dirc, robin, robin)
@test_throws MethodError BoundaryConditions(dirc, neum, dirc)
@test_throws MethodError BoundaryConditions(dirc, neum, neum)