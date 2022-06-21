# Modeling.jl

This is my Julia package for systems modeling

## State Space
---

This is a set of functions for State Space modeling techniques of the form

$$
\begin{align}
\dot x &= Ax + Bu \\
y &= Cx + Du
\end{align}
$$

that when discretized becomes

$$
\begin{align}
y_k &= Cx_k + Du_k \\
x_{k + 1} &= A_dx_k + B_dx_k
\end{align}
$$

where $A_d$ and $B_d$ are defined as follows

$$
\begin{align}
A_d &= e^{A\Delta t} \\
B_d &= A^{-1}\left(A_d - I\right)B \quad \mathrm{if }\det(A) \neq 0 \\
B_d &= \frac{\Delta t}{2}\left(I + A_d\right)B \quad \mathrm{if } \det(A) = 0
\end{align}
$$

Structure is defined as follows

```julia
struct StateSpace
    A::Matrix{<:Real}
    B::Vector{<:Real}
    C::Matrix{<:Real}
    D::Real
end
```
Where `D` is by default 0

Example for mass spring damper model

```julia
A = [0 1; -1 -1]
B = [1, 0]
C = [1 0]
ss = StateSpace(A, B, C)
```

#### Discretized
---
Full discretized example from spring-mass-damper system above

```julia
A = [0 1; -1 -1]
B = [1, 0]
C = [1 0]
ss = StateSpace(A, B, C)
dt = 0.01
t = 0:dt:10
dss = DiscreteSS(ss, dt)
x_i = [10., 0.]
x = [x_i]

for i in t[2:end]
    push!(x, next_x(dss, last(x), 0.))
end
```


## Modeling Techniques
---
These both work for State Space and other models

#### Euler's Method
---
Based on a linear approximation between time steps

`eulers` function in Julia is defined as follows

```julia
eulers(xdot::Function, x::Vector{<:Real}, dt::Real, params, u::Real)
```

#### 4th Order Runge-Kutta Method
---
Much more accurate than Euler's method, and maybe more accurate than discretized State Space

The `rk4` function would look like this

```julia
rk4(xdot::Function, x::Vector{<:Real}, dt::Real, params, u::Real)
```
where `xdot` should be defined as follows

```julia
xdot(x::Vector{<:Real}, u::Real)
```
and returns the derivative of the vector `x`

The `rk4` function would look like this for state space models
```julia
rk4(ss::StateSpace, x::Vector{<:Real}, dt::Real, u::Real = 0.)
```

Where `ss` is defined by using the `StateSpace` struct above (using Mass-spring-damper as example)
```julia
A = [0 1; -1 -1]
B = [1, 0]
C = [1 0]
ss = StateSpace(A, B, C)
x_next = rk4(ss, [10., 0.], 1.0)
```
