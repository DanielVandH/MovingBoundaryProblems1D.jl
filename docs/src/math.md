```@meta
CurrentModule = MovingBoundaryProblems1D
```

# Mathematical Details 

In this section, we provide some of the mathematical details for discretising the PDEs we consider. Recall that the PDEs we consider take the form:

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

It is assumed that $L(t) > 0$ for $t \geq 0$. (Note that, in this discussion, $t \geq 0$, but the package allows for generic timespans.) If you need other types of problems, e.g. $L(t)$ is negative somewhere, you might need to apply some transformations to the domain first.

## Transforming onto a Fixed Domain

We apply the Landau transform $(\xi, \tau) \mapsto (x/L(t), t)$ to transform the moving boundary problem onto a fixed domain. Writing $q(\xi, \tau) = u(x, t)$, the chain rule gives

```math
\begin{align*}
\dfrac{\partial u}{\partial t} &= \dfrac{\partial q}{\partial \xi}\dfrac{\partial \xi}{\partial t} + \dfrac{\partial q}{\partial\tau}\dfrac{\mathrm d\tau}{\mathrm dt} \\
&= -\frac{x}{L^2}\dfrac{\mathrm dL}{\mathrm dt}\dfrac{\partial q}{\partial \xi} + \dfrac{\partial q}{\partial \tau} \\
\dfrac{\partial q}{\partial\tau} &= \dfrac{\partial u}{\partial t} + \frac{\xi}{L}\dfrac{\mathrm dL}{\mathrm dt}\dfrac{\partial q}{\partial \xi} \\
&= \dfrac{\partial}{\partial x}\left(D\dfrac{\partial u}{\partial x}\right) + R + \frac{\xi}{L}\dfrac{\mathrm dL}{\mathrm dt}\dfrac{\partial q}{\partial \xi} \\
&= \frac{\xi}{L}\dfrac{\mathrm dL}{\mathrm dt}\dfrac{\partial q}{\partial \xi} + \frac{1}{L^2}\dfrac{\partial}{\partial \xi}\left(D\dfrac{\partial q}{\partial \xi}\right).
\end{align*}
```

Similarly, $\partial u(0, t)/\partial x = (1/L)\partial q(0, t)/\partial\xi$ so $\partial q(0, t)/\partial\xi = La_0$ and $\partial q(1, t)/\partial\xi = La_1$. Lastly, $\mathrm dL/\mathrm dt = a_2 + (b_2/L)\partial q(1, t)/\partial\xi$. Thus, the problem to discretise is (replacing $\tau$ with $t$):


```math
\begin{align*}
\begin{array}{rcll}
\dfrac{\partial q}{\partial t} & = & \dfrac{\xi}{L}\dfrac{\mathrm dL}{\mathrm dt}\dfrac{\partial q}{\partial \xi} + \dfrac{1}{L^2}\dfrac{\partial}{\partial \xi}\left(D\left(q, \xi L, t\right)\dfrac{\partial q}{\partial \xi}\right) + R\left(q, \xi L, t\right) & 0 < \xi < 1,\, t>0, \\[9pt]
\dfrac{\partial q}{\partial \xi} & = & La_0(q, t) & \xi = 0,\, t>0, \\[9pt]
\dfrac{\partial q}{\partial \xi} & = & La_1(q, t) & \xi = 1,\, t>0, \\[9pt]
\dfrac{\mathrm dL}{\mathrm dt} & = & a_2 + \dfrac{b_2}{L}\dfrac{\partial q}{\partial \xi} & \xi = 1,\, t>0, \\[9pt]
q(\xi, 0) & = & u_0\left(\xi L_0\right) & 0 \leq \xi \leq 1,\\[9pt] 
L(0) & = & L_0.
\end{array}
\end{align*}
```

Note that we could replace $\partial q/\partial \xi$ in the $\mathrm dL/\mathrm dt$ condition with $La_1(q, t)$, but in case we have a Dirichlet boundary condition for $q(1, t)$ rather than a Neumann boundary condition, a finite difference would be needed, so we delay this point.

## Interior Discretisation

We now proceed by discretising the fixed boundary problem. We define some grid $\xi_1,\ldots,\xi_n$ over $[0, 1]$, where $0 = \xi_1 < \xi_2 < \cdots < \xi_n = 1$, and let

```math
\begin{align*}
w_i &= \begin{cases} \xi_1 & i=1, \\[6pt] \dfrac12\left(\xi_{i-1} + \xi_i\right) & i=2,\ldots,n, \end{cases} \\[6pt]
e_i &= \begin{cases} \dfrac12\left(\xi_i + \xi_{i+1}\right) & i=1,\ldots,n-1, \\[6pt] \xi_n & i = n. \end{cases}
\end{align*}
```

Next, we define $V_i = e_i - w_i$, $i = 1,\ldots,n$, and $h_i = \xi_{i+1} - \xi_i$, $i=1,\ldots,n-1$. Now integrate the PDE over a control volume $[w_i, e_i]$:

```math
\begin{align*}
\int_{w_i}^{e_i} \dfrac{\partial q}{\partial t}\,\mathrm{d}\xi &= \frac{1}{L}\dfrac{\mathrm dL}{\mathrm dt}\int_{w_i}^{e_i} \xi\dfrac{\partial q}{\partial \xi}\,\mathrm{d}\xi + \int_{w_i}^{e_i} R(q, \xi L, t)\,\mathrm{d}\xi \\
&+ \frac{1}{L^2}\left\{D\left(q(e_i, t), e_iL, t\right)\dfrac{\partial q(e_i, t)}{\partial \xi} - D\left(q(w_i, t), w_iL, t\right)\dfrac{\partial q(w_i, t)}{\partial \xi}\right\}. 
\end{align*}
```

This first integral can be handled using integration by parts,

```math
\begin{align*}
\int_{w_i}^{e_i} \xi\dfrac{\partial q}{\partial \xi}\,\mathrm{d}\xi &= \left[\xi q\right]_{w_i}^{e_i} - \int_{w_i}^{e_i} q\,\mathrm{d}\xi \\
&= e_iq(e_i, t) - w_iq(w_i, t) - \int_{w_i}^{e_i} q\,\mathrm{d}\xi.
\end{align*}
```

Now, we also let $\bar q_i = (1/V_i)\int_{w_i}^{e_i} q\,\mathrm d\xi$ and $\bar R_i = (1/V_i)\int_{w_i}^{e_i} R\,\mathrm d\xi$. With this, we obtain:

```math
\begin{align*}
V_i\dfrac{\mathrm d\bar q_i}{\mathrm dt} &= \dfrac{1}{L}\dfrac{\mathrm dL}{\mathrm dt}\left\{e_iq(e_i, t) - w_iq(w_i, t) - V_i\bar q_i\right\} +V_i \bar R_i \\
&+ \frac{1}{L^2}\left\{D\left(q(e_i, t), e_iL, t\right)\dfrac{\partial q(e_i, t)}{\partial \xi} - D\left(q(w_i, t), w_iL, t\right)\dfrac{\partial q(w_i, t)}{\partial \xi}\right\} \\
\dfrac{\mathrm d\bar q_i}{\mathrm dt} &= \dfrac{1}{V_iL}\dfrac{\mathrm dL}{\mathrm dt}\left[e_iq(e_i, t) - w_iq(w_i, t)\right] - \dfrac{1}{L}\dfrac{\mathrm dL}{\mathrm dt}\bar q_i + \bar R_i \\
&+ \frac{1}{V_iL^2}\left\{D\left(q(e_i, t), e_iL, t\right)\dfrac{\partial q(e_i, t)}{\partial \xi} - D\left(q(w_i, t), w_iL, t\right)\dfrac{\partial q(w_i, t)}{\partial \xi}\right\}. 
\end{align*}
```

To now proceed, we let $q_i = q(\xi_i, t)$, $D_i = D(q_i, \xi_i L, t)$, and $R_i = R(q_i, \xi_i L, t)$. We then make the approximations:

```math
\begin{align*}
\begin{array}{rcll}
\bar q_i &\approx &  q_i & i=1,\ldots,n, \\[9pt]
\bar R_i &\approx & R_i & i=1,\ldots,n,\\[9pt]
q(e_i, t) &\approx& \dfrac12\left(q_i + q_{i+1}\right) & i=1,\ldots,n-1,\\[9pt]
q(w_i, t) &\approx & \dfrac12\left(q_{i-1} + q_i\right) & i=2,\ldots,n, \\[9pt]
D\left(q(e_i, t), e_iL, t\right) &\approx& \dfrac12\left(D_i + D_{i+1}\right) \quad & i=1,\ldots,n-1,\\[9pt]
D\left(q(w_i, t), w_iL, t\right) &\approx& \dfrac12\left(D_{i-1} + D_i\right) & i=2,\ldots,n\\[9pt] 
\dfrac{\partial q(e_i, t)}{\partial \xi} &\approx& \dfrac{q_{i+1} - q_i}{h_i} & i=1,\ldots,n-1, \\[9pt]
\dfrac{\partial q(w_i, t)}{\partial \xi} &\approx& \dfrac{q_i - q_{i-1}}{h_{i-1}} & i=2,\ldots,n.
\end{array}
\end{align*}
```

With these approximations, our approximation of the PDE at the interior nodes becomes

```math
\begin{align*}
\dfrac{\mathrm dq_i}{\mathrm dt} &= \frac{1}{V_iL}\dfrac{\mathrm dL}{\mathrm dt}\left[e_i\left(\dfrac{q_i+q_{i+1}}{2}\right) - w_i\left(\dfrac{q_{i-1} + q_i}{2}\right)\right] - \dfrac{1}{L}\dfrac{\mathrm dL}{\mathrm dt}q_i + R_i \\
&+ \frac{1}{V_iL^2}\left[\left(\dfrac{D_i + D_{i+1}}{2}\right)\left(\dfrac{q_{i+1} - q_i}{h_i}\right) - \left(\dfrac{D_{i-1} + D_i}{2}\right)\left(\dfrac{q_i - q_{i-1}}{h_{i-1}}\right)\right]
\end{align*}
```

## Boundary Discretisation

Now we handle discretising the boundary conditions. The boundary condition at $\xi = 0$ gives $\partial q/\partial \xi = La_0(q_1, t)$, so, also using $w_1 = \xi_1 = 0$,

```math
\begin{align*}
\dfrac{\mathrm dq_1}{\mathrm dt} &= \dfrac{1}{V_1L}\dfrac{\mathrm dL}{\mathrm dt}e_1\left(\dfrac{q_1 + q_2}{2}\right) - \dfrac{1}{L}\dfrac{\mathrm dL}{\mathrm dt}q_1 + R_1 \\
&+ \frac{1}{V_1L^2}\left[\left(\dfrac{D_1 + D_2}{2}\right)\left(\dfrac{q_2 - q_1}{h_1}\right) - D_1La_0(q_1, t)\right].
\end{align*}
```

The boundary condition at $\xi = 1$ gives $\partial q/\partial \xi = La_1(q_n, t)$, so, also using $e_n = \xi_n = 1$,

```math
\begin{align*}
\dfrac{\mathrm dq_n}{\mathrm dt} &= \dfrac{1}{V_nL}\dfrac{\mathrm dL}{\mathrm dt}\left[q_n - w_n\left(\dfrac{q_{n-1} + q_n}{2}\right)\right] - \dfrac{1}{L}\dfrac{\mathrm dL}{\mathrm dt}q_n + R_n \\
&+ \dfrac{1}{V_nL^2}\left[D_nLa_1(q_n, t) - \left(\dfrac{D_{n-1} + D_n}{2}\right)\left(\dfrac{q_n - q_{n-1}}{h_{n-1}}\right)\right].
\end{align*}
```

Now we handle $\mathrm dL/\mathrm dt$. If we have a boundary condition $\partial q/\partial \xi = La_1(q_n, t)$ for $q(1, t)$, then we can write

```math
\begin{align*}
\dfrac{\mathrm dL}{\mathrm dt} &= a_2(q_n, t) + \frac{b_2(q_n, t)}{L}La_1(q_n, t) \\
&= a_2(q_n, t) + a_1(q_n, t)b_2(q_n, t).
\end{align*}
```

If instead we have a Dirichlet boundary condition on $q(1, t)$, we use a three-point finite difference, writing (see e.g. [this paper](http://www.m-hikari.com/ijma/ijma-password-2009/ijma-password17-20-2009/bhadauriaIJMA17-20-2009.pdf)).

```math
\dfrac{\partial q}{\partial \xi} \approx \dfrac{h_{n-1}}{h_{n-2}\left(h_{n-2} + h_{n-1}\right)}q_{n-2} - \dfrac{h_{n-2} + h_{n-1}}{h_{n-2}h_{n-1}}q_{n-1} + \dfrac{h_{n-2} + 2h_{n-1}}{h_{n-2}\left(h_{n-1} + h_{n-2}\right)}q_n,
```

so that

```math
\begin{align*}
\dfrac{\mathrm dL}{\mathrm dt} & = a_2(q_n, t) + \dfrac{b_2(q_n, t)}{L}\left(\dfrac{h_{n-1}}{h_{n-2}\left(h_{n-2} + h_{n-1}\right)}q_{n-2} - \dfrac{h_{n-2} + h_{n-1}}{h_{n-2}h_{n-1}}q_{n-1} + \dfrac{h_{n-2} + 2h_{n-1}}{h_{n-2}\left(h_{n-1} + h_{n-2}\right)}q_n\right).
\end{align*}
```

## The Complete Discretisation 

Putting all these results together, the complete system of ODEs is

```math
\begin{align*}
\dfrac{\mathrm dq_i}{\mathrm dt} &= \frac{1}{V_iL}\dfrac{\mathrm dL}{\mathrm dt}\left[e_i\left(\dfrac{q_i+q_{i+1}}{2}\right) - w_i\left(\dfrac{q_{i-1} + q_i}{2}\right)\right] - \dfrac{1}{L}\dfrac{\mathrm dL}{\mathrm dt}q_i + R_i \\
&+ \frac{1}{V_iL^2}\left[\left(\dfrac{D_i + D_{i+1}}{2}\right)\left(\dfrac{q_{i+1} - q_i}{h_i}\right) - \left(\dfrac{D_{i-1} + D_i}{2}\right)\left(\dfrac{q_i - q_{i-1}}{h_{i-1}}\right)\right], \\& \quad\quad \quad i=2,\ldots,n-1, \\[9pt]
\dfrac{\mathrm dq_1}{\mathrm dt} &= \dfrac{1}{V_1L}\dfrac{\mathrm dL}{\mathrm dt}e_1\left(\dfrac{q_1 + q_2}{2}\right) - \dfrac{1}{L}\dfrac{\mathrm dL}{\mathrm dt}q_1 + R_1 \\
&+ \frac{1}{V_1L^2}\left[\left(\dfrac{D_1 + D_2}{2}\right)\left(\dfrac{q_2 - q_1}{h_1}\right) - D_1La_0(q_1, t)\right], \\[9pt] 
\dfrac{\mathrm dq_n}{\mathrm dt} &= \dfrac{1}{V_nL}\dfrac{\mathrm dL}{\mathrm dt}\left[q_n - w_n\left(\dfrac{q_{n-1} + q_n}{2}\right)\right] - \dfrac{1}{L}\dfrac{\mathrm dL}{\mathrm dt}q_n + R_n \\
&+ \dfrac{1}{V_nL^2}\left[D_nLa_1(q_n, t) - \left(\dfrac{D_{n-1} + D_n}{2}\right)\left(\dfrac{q_n - q_{n-1}}{h_{n-1}}\right)\right], \\[9pt]
\dfrac{\mathrm dL}{\mathrm dt} &= a_2(q_n, t) + a_1(q_n, t)b_2(q_n, t).
\end{align*}
```

This system is solved, and then from $q(\xi_i, t_j)$ we get values of $u$ at $(x_iL(t_j), t_j)$.