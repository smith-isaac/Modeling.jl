# Modeling.jl

This is my Julia package for systems modeling

## State Space
---
$$
\begin{gathered}
\dot x = Ax + Bu \\
y = Cx + Du \to y = Cx
\end{gathered}
$$

Struct is defined as follows

```julia
struct StateSpace
    A::Matrix{Float64}
    B::Vector{Float64}
    C::Matrix{Float64}
end
```


Example for mass spring damper model $\left(\ddot x + \dot x + x = F\right)$
$$
A = 
\begin{bmatrix}
0 & 1 \\ -1 & -1
\end{bmatrix}
,
B =
\begin{bmatrix}
1 \\ 0
\end{bmatrix}
,C =
\begin{bmatrix} 1 & 0
\end{bmatrix}
$$

```julia
A = [0 1; -1 -1]
B = [1, 0]
C = [1 0]
ss = StateSpace(A, B, C)
```

#### Discretized
---
This uses an exponential model in between each time-step instead of linear model
$$
\dot x = Ax \to x = e^{\int A dt}
$$

$$
\begin{gathered}
A_d = e^{A\Delta t} \\
B_d = A^{-1} \left(A_d - I\right)B \quad \textrm{or} \quad B_d = \frac{\Delta t}{2}\left(I + A_d\right) \quad \textrm{if}\quad \det A = 0
\end{gathered}
$$

Then iterating through the time steps
$$
\begin{gathered}
y_k = Cx_k \\
x_{k+1} = A_d x_k + B_du_k
\end{gathered}
$$

Full discretized example from spring-mass-damper system above $\left(\ddot x + \dot x + x = F\right)$

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
$$
\dot x = f(x, u) \quad\textrm{or}\quad \dot x = Ax + Bu
$$

#### Euler's Method
---
Based on a linear approximation between time steps
$$
\begin{gathered}
x_{k + 1} = \dot x_k \Delta t + x_k \\
x_{k + 1} = f(x_k, u_k)\Delta t + x_k \quad\textrm{or}\quad x_{k + 1} = (A x_k + Bu_k)\Delta t + x_k
\end{gathered}
$$

`eulers` function in Julia is defined as follows

```julia
eulers(xdot::Function, x::Vector{Float64}, params, u::Float64)
```

#### 4th Order Runge-Kutta Method
---
Much more accurate than Euler's method, and maybe more accurate than discretized State Space

Definition comes from the differential equation being of the form $\dot x = f(x, t)$ but can be expanded to fit our models.

$$
\begin{aligned}
\dot x &= f(x, t)\qquad
h = \Delta t \\
k_1 &= f(x_n, t_n) \\
k_2 &= f(x_n + k_1\frac{h}{2}, t + \frac{h}{2}) \\
k_3 &= f(x_n + k_2\frac{h}{2}, t + \frac{h}{2}) \\
k_4 &= f(x_n + k_3h, t + h) \\
x_{n + 1} &= x_n + \frac{h}{6} (k_1 + 2k_2 + 2k_3 + k_4)
\end{aligned}
$$

For differential models such as $(\dot x = f(x, u))$ this looks like the following
$$
\begin{aligned}
k_1 &= f(x_n, u_k) \\
k_2 &= f(x_n + k_1\frac{h}{2},u_k) \\
k_3 &= f(x_n + k_2\frac{h}{2},u_k) \\
k_4 &= f(x_n + k_3h,u_k) \\
x_{n + 1} &= x_n + \frac{h}{6} (k_1 + 2k_2 + 2k_3 + k_4)
\end{aligned}
$$

The function would look like this

```julia
rk4(xdot::Function, x::Vector{Float64}, params, u::Float64)
```

For State Space models $(\dot x = Ax + Bu)$ it would look like this

$$
\begin{aligned}
k_1 &= Ax_n + Bu_k \\
k_2 &= A(k_1\frac{h}{2}) + Bu_k \\
k_3 &= A(k_2\frac{h}{2}) + Bu_k \\
k_4 &= A(k_3h) + Bu_k \\
x_{n + 1} &= x_n + \frac{h}{6} (k_1 + 2k_2 + 2k_3 + k_4)
\end{aligned}
$$
The function would look like this
```julia
rk4(ss::StateSpace, x::Vector{Float64}, u::Float64)
```
