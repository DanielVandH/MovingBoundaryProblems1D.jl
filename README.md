# MovingBoundaryProblems1D

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://DanielVandH.github.io/MovingBoundaryProblems1D.jl/dev/)
[![Build Status](https://github.com/DanielVandH/MovingBoundaryProblems1D.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DanielVandH/MovingBoundaryProblems1D.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/DanielVandH/MovingBoundaryProblems1D.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/DanielVandH/MovingBoundaryProblems1D.jl)

This is a package for solving one-dimensional single phase moving boundary problems using the finite volume method. The problems that can be solved using this package take the form

```math
\begin{align*}
\begin{array}{rcll}
\dfrac{\partial u(x, t)}{\partial t} & = & \dfrac{\partial}{\partial x}\left(D\left(u(x, t), x, t\right)\dfrac{\partial u(x, t)}{\partial x}\right) + R\left(u(x, t), x, t\right) & 0 < x < L(t),\, t>0, \\[9pt]
\dfrac{\partial u(0, t)}{\partial x} & = & a_0\left(u(0, t), t\right) & x = 0,\,t>0, \\[9pt]
\dfrac{\partial u(L(t), t)}{\partial x} & = & a_1\left(u(L(t), t)\right) & x=L(t),\,t>0, \\[9pt]
\dfrac{\mathrm dL}{\mathrm dt} & = & a_2\left(u(L(t), t), t\right) + b_2\left(u(L(t), t), t\right)\dfrac{\partial u(L(t), t)}{\partial x} & x = L(t),\,t>0, \\[9pt]
u(x, 0) & = & u_0(x) & 0 \leq x \leq L(0), \\[9pt]
L(0) & = & L_0.
\end{array}
\end{align*}
```

You can also use Dirichlet boundary conditions rather than Neumann boundary conditions for $u(0, t)$ and $u(L(t), t)$. 

For examples on how to use it, please see the docs. 