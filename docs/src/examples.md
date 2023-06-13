```@meta
CurrentModule = MovingBoundaryProblems1D
```

# Examples 

This section gives some examples for how the package can be used. 

## Example I: Heat equation

We start with a simple problem that has an exact solution, as discussed in e.g. [this paper](https://doi.org/10.1016/S0377-0427(97)00034-4). The problem is a heat equation with a moving boundary (or a one-phase classical Stefan problem):

```math
\begin{align*}
\begin{array}{rcll}
\dfrac{\partial u}{\partial t} & = & \dfrac{\partial^2u}{\partial x^2} & 0 < x < L(t), \\[9pt]
\dfrac{\partial u}{\partial t} & = & -\exp(t) & x = 0,\, t > 0, \\[9pt]
u & = & 0 & x = L(t),\,t>0, \\[9pt]
\dfrac{\mathrm dL}{\mathrm dt} & = & -\dfrac{\partial u}{\partial x} & x = L(t),\,t>0, \\[9pt] 
L(0) & = & 0 \\[9pt]
u(x, 0) & = & 0 & 0 \leq x \leq L(0) = 0.
\end{array}
\end{align*}
```

The exact solution is $u(x, t) = \exp(t - x) - 1$, $L(t) = t$ for $0 \leq x \leq L(t)$ and $0 < t < 1$. Since we need to have $L(0) > 0$, and the domain is just $\{0\}$ at $t = 0$, we solve the problem starting at $t = 0.1$ rather than $t = 0$. We setup the problem as follows:

```julia
using MovingBoundaryProblems1D 

# Define the exact solutions for getting the initial data
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
```

The problem is then solved with `solve`, as is done in DifferentialEquations.jl.

```julia
using OrdinaryDiffEq, LinearSolve
sol = solve(prob, TRBDF2(linsolve=KLUFactorization()), saveat=0.1)
```

Now, for plotting the solution, we cannot simply use `mesh_points` since those are defined on $[0, 1]$. To obtain the scaled mesh points at each time, we use `scaled_mesh_points`:

```julia
julia> scaled_mesh = scaled_mesh_points(prob, sol)
5-element Vector{Vector{Float64}}:
 [0.0, 0.0002004008016032064, 0.0004008016032064128, 0.0006012024048096193, 0.0008016032064128256, 0.001002004008016032, 0.0012024048096192386, 0.001402805611222445, 0.0016032064128256513, 0.001803607214428858  …  0.09819639278557114, 0.09839679358717435, 0.09859719438877756, 0.09879759519038077, 0.09899799599198397, 0.09919839679358718, 0.09939879759519038, 0.0995991983967936, 0.0997995991983968, 0.1]
 [0.0, 0.00040080095829667935, 0.0008016019165933587, 0.0012024028748900382, 0.0016032038331867174, 0.002004004791483397, 0.0024048057497800764, 0.002805606708076756, 0.003206407666373435, 0.003607208624670115  …  0.1963924695653729, 0.19679327052366957, 0.19719407148196627, 0.19759487244026294, 0.1979956733985596, 0.1983964743568563, 0.19879727531515298, 0.19919807627344965, 0.19959887723174635, 0.19999967819004302]
 [0.0, 0.0006012011247570208, 0.0012024022495140416, 0.0018036033742710623, 0.002404804499028083, 0.003006005623785104, 0.0036072067485421245, 0.004208407873299146, 0.004809608998056166, 0.005410810122813187  …  0.2945885511309402, 0.2951897522556972, 0.29579095338045425, 0.2963921545052112, 0.29699335562996826, 0.2975945567547253, 0.2981957578794823, 0.2987969590042393, 0.29939816012899634, 0.29999936125375337]
 [0.0, 0.0008016010620834612, 0.0016032021241669224, 0.0024048031862503837, 0.0032064042483338447, 0.004008005310417306, 0.004809606372500767, 0.005611207434584229, 0.0064128084966676895, 0.007214409558751151  …  0.392784520420896, 0.39358612148297945, 0.39438772254506294, 0.39518932360714637, 0.39599092466922986, 0.3967925257313133, 0.3975941267933968, 0.3983957278554802, 0.3991973289175637, 0.39999892997964714]
 [0.0, 0.001002000880061787, 0.002004001760123574, 0.0030060026401853607, 0.004008003520247148, 0.005010004400308935, 0.006012005280370721, 0.007014006160432509, 0.008016007040494296, 0.009018007920556083  …  0.4909804312302756, 0.4919824321103374, 0.4929844329903992, 0.493986433870461, 0.49498843475052273, 0.49599043563058454, 0.49699243651064634, 0.49799443739070814, 0.4989964382707699, 0.4999984391508317]
```

Let us now make the plot. Notice that in the code we need to use `u[begin:(end-1)]` to get the values of the solution rather than simply `u`, since `u[end]` is where we store the value of $L$.

```julia
fig = Figure(fontsize=35)
ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"u(x, t)", width=600, height=300)
colors = [:red, :blue, :black, :magenta, :darkgreen]
[lines!(ax, x, u[begin:(end-1)], color=clr, linewidth=1.5) for (x, u, clr) in zip(scaled_mesh, sol.u, colors)]
[lines!(ax, x, exact_u.(x, t), color=clr, linestyle=:dash, linewidth=3) for (x, t, clr) in zip(scaled_mesh, sol.t, colors)]
resize_to_layout!(fig)
```

```@raw html
<figure>
    <img src='../figures/heat_equation.png', alt'Solution to the heat equation problem'><br>
</figure>
```

## Example II: One-phase Stefan problem and steady states

We now consider a second example. This problem is similar to the previous example, except we do not start with $L(0) = L_0$, we do not have an exact solution, and we will examine steady states. The problem is:

```math
\begin{array}{rcll}
\dfrac{\partial u}{\partial t} & = & D\dfrac{\partial^2 u}{\partial x^2} & 0 < x < L(t),\, t> 0, \\[9pt]
\dfrac{\partial u}{\partial x} & = & 0 & x = 0,\,t>0,\\[9pt]
u & = & 0 & x = L(t),\,t>0,\\[9pt]
\dfrac{\mathrm dL}{\mathrm dt} & = & -\kappa\dfrac{\partial u}{\partial x} & x = L(t).
\end{array}
```

As discussed in e.g. Eq. 1.6. of [this paper](https://doi.org/10.1098/rspa.2019.0378), $u(x, t) \to 0$ and $L(t) \to L_e$ as $t \to \infty$, where

```math
L_e = L(0) + \dfrac{\kappa}{D}\int_0^{L(0)} u(x, 0)\,\mathrm dx.
```

Out of interest, let us demonstrate how this is derived. Let us write

```math
\begin{align*}
N(t) &= \int_0^{L(t)} u(x, t)\,\mathrm dx \\[9pt]
\dfrac{\mathrm dN}{\mathrm dt} &= \dfrac{\mathrm dL}{\mathrm dt}\underbrace{u(L(t), t)}_0 + \int_0^{L(t)} \dfrac{\partial u}{\partial t}\,\mathrm dx \\[9pt]
&= \int_0^{L(t)} D\dfrac{\partial^2 u}{\partial x^2}\,\mathrm dx \\[9pt]
&= D\left[\dfrac{\partial u(L(t), t)}{\partial x} - \dfrac{\partial u(0, t)}{\partial x}\right] \\[9pt]
&= -\dfrac{D}{\kappa}\dfrac{\mathrm dL}{\mathrm dt}.
\end{align*}
```

Now integrating, we obtain $N(t) = -(D/\kappa)L(t) + C$ for some constant $C$. To evaluate this constant $C$, we integrate the initial condition:

```math
\begin{align*}
\int_0^{L(0)} u(x, 0)\,\mathrm dx &= -\dfrac{D}{\kappa}L(t) + C\\[9pt]
C &= \dfrac{D}{\kappa}L(0) + \int_0^{L(0)} u(x, 0)\,\mathrm dx.
\end{align*}
```

Thus,

```math
N(t) = \dfrac{D}{\kappa}\left(L(0) - L(t)\right) + \int_0^{L(0)} u(x, 0)\,\mathrm dx.
```

Now take $t \to \infty$ and note that $N(t) \to 0$, so that we get

```math
\begin{align*}
0 &= \dfrac{D}{\kappa}\left(L(0) - L_e\right) + \int_0^{L(0)} u(x, 0)\,\mathrm dx \\
0 &= L(0) - L_e + \dfrac{\kappa}{D}\int_0^{L(0)} u(x, 0)\,\mathrm dx \\
L_e &= L(0) + \dfrac{\kappa}{D}\int_0^{L(0)} u(x, 0)\,\mathrm dx.
\end{align*}
```

We will explore $L_e$ later. Let us start by now solving the problem. The initial condition we use is 

```math
u(x, 0) = \begin{cases} \alpha & x < \beta, \\ 0 & x \geq \beta,\end{cases}
```

where $\alpha = 1/2$ and $\beta = 1$. We start with $L(0) = \beta$. With these values, and $D = 0.1$ and $\kappa = 0.1$, we find

```math
L_e = 1 + \dfrac{0.1}{0.1}\int_0^1 \alpha \,\mathrm dx,
```

or $L_e = 3/2$. We can verify this estimate by solving the problem:

```julia
using MovingBoundaryProblems1D

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

using OrdinaryDiffEq, LinearSolve 
sol = solve(prob, TRBDF2(linsolve=KLUFactorization()), saveat=final_time / 250)

using CairoMakie 
let x = (stack ∘ scaled_mesh_points)(prob, sol), t = repeat(sol.t, inner=length(mesh_points)), u = sol[begin:(end-1), :]
    fig = Figure(fontsize=33)
    ax = Axis(fig[1, 1], width=600, height=300, xlabel=L"x", ylabel=L"t")
    tricontourf!(ax, vec(x), t, vec(u))
    resize_to_layout!(fig)
end
```

```@raw html
<figure>
    <img src='../figures/heat_stefan.png', alt'Solution to the Stefan problem'><br>
</figure>
```

We can indeed see that the solution is zero for large time, and the highest value of $x$ is indeed around $L_e = 3/2$. 

We provide a method for computing the steady state of a moving boundary problem. By wrapping the `prob` in `SteadyMBProblem`, another `solve` will give us the steady state:

```julia
using SteadyStateDiffEq 
sprob = SteadyMBProblem(prob)
```

```julia-repl
julia> ssol = solve(sprob, DynamicSS(TRBDF2()))
u: 501-element Vector{Float64}:
 3.8102357269537885e-8
 3.8102170476792114e-8
 3.8101610100268955e-8
 3.810067614511076e-8
 3.8099368619888047e-8
 3.80976875365995e-8
 3.8095632910671864e-8
 ⋮
 6.020854039360681e-10
 4.816570994416517e-10
 3.61232957797933e-10
 2.4081442647367477e-10
 1.2040295600361738e-10
 0.0
 1.49942609487085
```

Indeed, the first $500$ components approximately constitute the zero vector, and the last component (where $L$ is stored) is essentially $L_e \approx 3/2$.

## Example III: Fisher-Stefan model, exploring the spreading-extinction dichotomy

This example considers the Fisher-Stefan model (see e.g. [this paper](https://doi.org/10.1098/rspa.2019.0378)). This model is given by (in non-dimensional form)

```math
\begin{array}{rcll}
\dfrac{\partial u}{\partial t} & = & \dfrac{\partial^2 u}{\partial x^2} + u(1-u) & 0 < x < L(t), \, t > 0,\\[9pt]
\dfrac{\partial u}{\partial x} & = & 0 & x = 0,\,t>0, \\[9pt]
u & = & 0 & x = L(t),\,t>0, \\[9pt]
\dfrac{\mathrm dL}{\mathrm dt} & = & -\kappa\dfrac{\partial u}{\partial x} & x = L(t),\, t>0, \\[9pt]
u(x, 0) & = & \begin{cases} \alpha & x < \beta, \\ 0 & x \geq \beta, \end{cases} & 0 \leq x \leq L(0), \\[9pt]
L(0) &=& \beta,
\end{array}
```

A function that solves this for given $(\kappa, \alpha, \beta, T)$, where $T$ is the final time, and also returns the steady state is given below:

```julia
using MovingBoundaryProblems1D
using OrdinaryDiffEq, LinearSolve 
using SteadyStateDiffEq 

function construct_problem(𝝟, α, β, T, n=500)
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
```

Using this function, let us look at the solutions for $\kappa=20$ and $\kappa=0.45$, taking $\alpha=1/2$ and $\beta=1$.

```julia
prob1, x1, u1, t1, L1, su1, sL1, sx1, mL1, msol1 = solve_problem(20.0, 0.5, 1.0, 20.0);
prob2, x2, u2, t2, L2, su2, sL2, sx2, mL2, msol2 = solve_problem(0.45, 0.5, 1.0, 20.0);

using CairoMakie
fig = Figure(fontsize=33)
idx_rng = (1, 50, 100, 150, 200, 250)
colors = (:red, :blue, :black, :magenta, :darkgreen, :orange)
ax1 = Axis(fig[1, 1], width = 600, height = 300, xlabel=L"x",ylabel=L"u(x, t)", title=L"(a): $\kappa = 20$", titlealign=:left)
ax2 = Axis(fig[1, 2], width = 600, height = 300, xlabel=L"x",ylabel=L"u(x, t)", title=L"(b): $\kappa = 0.45$", titlealign=:left)
ax3 = Axis(fig[2, 1], width = 600, height = 300, xlabel=L"x",ylabel=L"t", title=L"(c): $\kappa = 20$", titlealign=:left)
ax4 = Axis(fig[2, 2], width = 600, height = 300, xlabel=L"x",ylabel=L"t", title=L"(d): $\kappa = 20$", titlealign=:left)
[lines!(ax1, x1[:, idx], u1[:, idx], color=clr, linewidth=1.5) for (idx, clr) in zip(idx_rng, colors)]
[lines!(ax2, x2[:, idx], u2[:, idx], color=clr, linewidth=1.5) for (idx, clr) in zip(idx_rng, colors)]
tricontourf!(ax3, vec(x1), repeat(t1, inner=size(x1, 1)), vec(u1), levels = 100)
tricontourf!(ax4, vec(x2), repeat(t2, inner=size(x2, 1)), vec(u2), levels = 20)
resize_to_layout!(fig)
```

```@raw html
<figure>
    <img src='../figures/fisher_stefan_1.png', alt'Solution to the Fisher-Stefan problem'><br>
</figure>
```

We see that with $\kappa=0.45$ the population goes extinct, but for $\kappa = 20$ the solution evolves to a travelling wave with some speed $c$. We can estimate the wave speed $c$ for $\kappa=20$ as follows:

```julia
using StatsBase
function estimate_wave_speed(sol)
    Ls = @views sol[end, (end ÷ 2) : end] # only take late time 
    ts = @views sol.t[(end ÷ 2) : end] 
    cs = diff(Ls) ./ diff(ts)
    return mean(cs)
end
```

```julia-repl
julia> c = estimate_wave_speed(msol1)
1.1860676450332943
```

This estimate of $c$ agrees with $c \approx 1.2$ reported by [El-Hachem et al. (2019)](https://doi.org/10.1098/rspa.2019.0378).

Let us now give a further study of the relationship between $c$ and $\kappa$. We start by reproducing Figure 5 in [El-Hachem et al. (2019)](https://doi.org/10.1098/rspa.2019.0378). 

```julia
function compute_wave_speed(κ, α, β, T)
    prob, x, u, t, L, su, sL, sx, mL, sol = solve_problem(κ, α, β, T, 2_000)
    if sL > 200.0 # the steady state solution makes no sense in this case, so we have a travelling wave
        c = estimate_wave_speed(sol)
        return c
    else
        return NaN
    end
end
```

Now we can reproduce the figure.

```julia
κ = LinRange(0.1, 2.0, 50)
c = compute_wave_speed.(κ, 0.5, 1.0, 5.0)
fig = Figure(fontsize=33)
ax = Axis(fig[1, 1], xlabel=L"\kappa", ylabel=L"c", width=600, height=300)
lines!(ax, κ, c, linewidth=3, color=:red)
resize_to_layout!(fig)
xlims!(ax, 0, 2)
ylims!(ax, 0, 1/2)
```

```@raw html
<figure>
    <img src='../figures/fisher_stefan_2.png', alt'c-k plot'><br>
</figure>
```

We see that there is a cutoff value of $\kappa$ called $\kappa_c$, such that the population goes extinct for $\kappa < \kappa_c$ and evolves to a travelling wave for $\kappa > \kappa_c$. We can estimate this $\kappa_c$:

```julia-repl
julia> κc = κ[findfirst(!isnan, c)]
0.48775510204081635
```

So, $\kappa_c \approx 0.48$, which agrees with the value in [El-Hachem et al. (2019)](https://doi.org/10.1098/rspa.2019.0378). 

We could go further and explore the crtiical length, verifying the results of [Du and Lin (2010)](http://dx.doi.org/doi:10.1137/090771089) and [El-Hachem et al. (2019)](https://doi.org/10.1098/rspa.2019.0378) to check that:

1. If $L(t) > \pi/2 - \kappa\int_0^{L(t)} u(x, t)\,\mathrm dx$ for at least one value of $t > 0$, the solution evolves to a travelling wave.
2. If $L(0) < \pi/2$, the population goes extinct and $L(0) < \lim_{t \to \infty} L(t) < \pi/2$.

We do not do this here, though.

## Example IV: Dirichlet boundary conditions and a parabolic initial condition 

We now consider

```math
\begin{align*}
\begin{array}{rcll}
\dfrac{\partial u}{\partial t} & = & \dfrac{\partial^2u}{\partial x^2} + u(1-u) & 0 < x < L(t), \, t>0, \\[9pt]
u(0, t) & = & 0 & t > 0, \\[9pt]
u(L(t), t) & = & 0 & t > 0, \\[9pt]
\dfrac{\mathrm dL(t)}{\mathrm dt} & = & -\dfrac{1}{2}\dfrac{\partial u(L(t), t)}{\partial t} & t > 0, \\[9pt]
u(x, 0) & = & 2x(1-x) & 0 \leq x \leq L(0), \\[9pt] 
L(0) & = & 1.
\end{array}
\end{align*}
```

This problem can be solved as follows:

```julia
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

let x = (stack ∘ scaled_mesh_points)(prob, sol), t = repeat(sol.t, inner=length(mesh_points)), u = sol[begin:(end-1), :]
    fig = Figure(fontsize=33)
    ax = Axis(fig[1, 1], width=600, height=300, xlabel=L"x", ylabel=L"t")
    tricontourf!(ax, vec(x), t, vec(u), levels=0.0:0.05:0.5, extendlow=:auto)
    resize_to_layout!(fig)
end
```

```@raw html
<figure>
    <img src='../figures/parabolic.png', alt'Solution to the parabolic problem'><br>
</figure>
```

## Example V: Epithelial dynamics 

We now consider a free boundary model of epithelial dynamics, following [Baker et al. (2019)](https://doi.org/10.1016/j.jtbi.2018.12.025) - refer to this paper for more detail about this problem. 

The PDE we take is

```math
\begin{align*}
\begin{array}{rcll}
\dfrac{\partial q}{\partial t} & = & \dfrac{\partial}{\partial x}\left(D(q)\dfrac{\partial q}{\partial x}\right) + qG\left(\dfrac{1}{q}\right) & 0 < x < L(t),\,t>0,\\[9pt]
\dfrac{\partial q}{\partial x} & = & 0 & x = 0,\,t>0,\\[9pt]
\dfrac{1}{\eta}F\left(\dfrac{1}{q}\right) + \dfrac{D(q)}{2q}\dfrac{\partial q}{\partial x} & = & 0 & x = L(t),\,t>0,\\[9pt] 
\dfrac{\mathrm dL}{\mathrm dt} & = & -\dfrac{D(q)}{q}\dfrac{\partial q}{\partial x} & x = L(t),\,t>0, \\[9pt]
q(x, 0) & = & q_0(x) & 0 \leq x \leq L(0).
\end{array}
\end{align*}
```

For the parameters of this problem, we take (these parameters follow a similar example in [Murphy et al. (2020)](https://doi.org/10.1007/s11538-020-00807-x)):

1. _Linear force law_: $F(q) = k(s-q)$ 
2. _Diffusion function_: $D(q) = k/(\eta q^2)$;
3. _Proliferation law_: $G(q) = \beta$;
4. _Spring constant_: $k = 10$;
5. _Resting spring length_: $s = 1$;
6. _Viscosity coefficient_: $\eta = 1$;
7. _Intrinsic proliferation rate_: $\beta = 0.00577$. (This value is chosen so that $\int_0^{400} q(x, 400)\,\mathrm dx \approx 400$.)

To define $q_0(x)$, we take a Gaussian initial density with variance three centred at $x = L_0/2$, where $L_0 = 10$, and we then scale it such that $N(0) = 40$, where $N(t) = \int_0^{L(t)} q(x, t)\,\mathrm dx$ is the number of cells at the time $t$. To find this $q_0(x)$, write

```math
q_0(x) = \dfrac{A}{\sqrt{2\pi \sigma^2}}\exp\left\{-\dfrac{1}{2}\left(\dfrac{x - L_0/2}{\sigma}\right)^2\right\},
```

where $\sigma^2 = 3$. Integrating this, we find that 

```math 
\int_0^{L_0} q_0(x)\,\mathrm dx = A\text{erf}\left(\dfrac{L_0}{4\sigma}\sqrt{2}\right).
``` 

So, setting $N(0) = A\text{erf}(L_0\sqrt{2}/4\sigma)$, we find

```math
q_0(x) = \dfrac{N(0)}{\text{erf}\left(\dfrac{L_0}{4\sigma}\sqrt{2}\right)\sqrt{2\pi\sigma^2}}\exp\left\{-\dfrac{\left(2x - L_0\right)^2}{8\sigma^2}\right\}.
```

The last thing to do before we define our problem is to rearrange the boundary condition at $x = L(t)$ into the form $\partial q/\partial x = a_2(q(L(t), t), t)$:

```math
\dfrac{\partial q(L(t), t)}{\partial x} = -\dfrac{2q}{\eta D(q)}F\left(\dfrac{1}{q}\right).
```

Let's now construct and solve the problem.

```julia
using MovingBoundaryProblems1D, SpecialFunctions 

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
using OrdinaryDiffEq, LinearSolve
sol = solve(prob, TRBDF2(linsolve=KLUFactorization()))
```

To visualise the solution, we will also look at the cell number $N(t)$. This cell number is computed using the following:

```julia
using DataInterpolations
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
```

We now plot.

```julia
using CairoMakie 

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
```

```@raw html
<figure>
    <img src='../figures/epithelial.png', alt'Solution to the epithelial problem'><br>
</figure>
```

