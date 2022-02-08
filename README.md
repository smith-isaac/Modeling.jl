# Modeling.jl

This is my Julia package for systems modeling

## State Space
---

Structure is defined as follows

```julia
struct StateSpace
    A::Matrix{Float64}
    B::Vector{Float64}
    C::Matrix{Float64}
    D::Float64
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
eulers(xdot::Function, x::Vector{Float64}, dt::Float64, params, u::Float64)
```

#### 4th Order Runge-Kutta Method
---
Much more accurate than Euler's method, and maybe more accurate than discretized State Space

The `rk4` function would look like this

```julia
rk4(xdot::Function, x::Vector{Float64}, dt::Float64, params, u::Float64)
```
where `xdot` should be defined as follows

```julia
xdot(x::Vector{Float64}, u::Float64)
```
and returns the derivative of the vector `x`

The `rk4` function would look like this for state space models
```julia
rk4(ss::StateSpace, x::Vector{Float64}, dt::Float64, u::Float64 = 0.)
```

Where `ss` is defined by using the `StateSpace` struct above (using Mass-spring-damper as example)
```julia
A = [0 1; -1 -1]
B = [1, 0]
C = [1 0]
ss = StateSpace(A, B, C)
x_next = rk4(ss, [10., 0.], 1.0)
```
