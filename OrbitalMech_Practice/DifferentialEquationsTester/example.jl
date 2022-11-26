# Testing DifferentialEquations package on Lorenz Attractor

cd(@__DIR__)
using Pkg
Pkg.activate(".")

# Original Lorenz attractor function
function lorenz(u, p, t)
    dx = 10.0 * (u[2] - u[1])
    dy = u[1] * (28.0 - u[3])
    dz = (u[1] * u[2]) - (8/3)*u[3]
    [dx, dy, dz]
end

using DifferentialEquations, BenchmarkTools

u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz, u0, tspan)

# Benchmark baseline
@benchmark solve(prob, Tsit5())
# Benchmark without saving the array every step
@benchmark solve(prob, Tsit5(), save_everystep = false)

# Altered Lorenz attractor function using a cache array, du. This makes the code non-allocating
function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3])
    du[3] = (u[1] * u[2]) - (8/3)*u[3]
end

prob1 = ODEProblem(lorenz!, u0, tspan)

@benchmark solve(prob1, Tsit5())
@benchmark solve(prob1, Tsit5(), save_everystep = false)