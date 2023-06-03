```@meta
CurrentModule = MovingBoundaryProblems1D
```

# Mathematical Details 

In this section, we provide some of the mathematical details for discretising the PDEs we consider. Recall that the PDEs we consider take the form:

```math
\begin{align*}
\begin{array}{rcll}
\dfrac{\partial u}{\partial t} & = & \dfrac{\partial}{\partial x}\left(D\left(u(x, t), x, t\right)\dfrac{\partial u}{\partial x}\right) + R(u) & 0 < x < L(t),\, t>0, \\[9pt]
a_0\left(u(0, t), t\right) + b_0\left(u(0, t), t\right)\dfrac{\partial u}{\partial x} & = & 0 & x = 0,\,t>0, \\[9pt]
a_1\left(u(L(t), t)\right) + b_1\left(u(L(t), t)\right)\dfrac{\partial u}{\partial x} & = & 0 & x=L(t),\,t>0, \\[9pt]
\dfrac{\mathrm dL}{\mathrm dt} & = & a_2\left(u(L(t), t)\right) + b_2\left(u(L(t), t)\right)\dfrac{\partial u}{\partial x} & x = L(t),\,t>0, \\[9pt]
u(x, 0) & = & u_0(x) & 0 \leq x \leq L(0), \\[9pt]
L(0) & = & L_0.
\end{array}
\end{align*}
```

It is assumed that $L(t) \neq 0$ for $t \geq 0$. (Note that, in this discussion, $t \geq 0$, but the package allows for generic timespans.)