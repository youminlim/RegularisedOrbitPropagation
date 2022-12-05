using Base: undef_ref_alignment
# Calculate satellite trajectory about a central body using the 2 body problem

cd(@__DIR__)
using Pkg
Pkg.activate(".")

using DifferentialEquations
using LinearAlgebra
using Plots

# Define 2 body problem for a central body at [0, 0, 0]
function f(du, u, p, t)
    r̲ = u[1:3]
    du[1:3] = u[4:6]

    G = 6.667e-20
    μ = G * 5.972e24

    r = norm(r̲)

    du[4:6] = -μ * 1000 * r̲ / r^3
end

# Define function parameters
u0 = [8000, 0, 6000, 0, 5, 5]
tspan = (0, 14400)
p = 0

problem = ODEProblem(f, u0, tspan, p)
sol = solve(problem, reltol = 1e-8, saveat = 1)

plot(sol, idxs = (1, 2, 3))