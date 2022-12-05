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
    du[1:6] = u[7:12]

    G = 6.667e-20
    μ = G * 5.972e24

    r = norm(r̲)

    du[7:9] = -μ * 1000 * r̲ / r^3
    du[10:12] .= 0
end

# Define function parameters
G = 6.667e-20
μ = G * 5.972e24
u0 = [8000.0, 0.0, 6000.0, 0.0, 0.0, 0.0, 0.0, 5.0, 5.0, 0.0, 0.0, 0.0]
tspan = (0, 14400)
p = 0

problem = ODEProblem(f, u0, tspan, p)
solution = solve(problem, reltol = 1e-8, saveat = 0.1)

solutionCOE = []
for i in solution.u
    push!(solutionCOE, coe(i, μ))
end
angularMomentum = getindex.(solutionCOE, 1)
angularMomentum = angularMomentum/maximum(angularMomentum)
inclination = getindex.(solutionCOE, 2)
inclination = inclination/maximum(inclination)
raan = getindex.(solutionCOE, 3)
raan = raan/maximum(raan)
eccentricity = getindex.(solutionCOE, 4)
eccentricity = eccentricity/maximum(eccentricity)
aop = getindex.(solutionCOE, 5)
trueAnomaly = getindex.(solutionCOE, 6)


# Plotting
    # Plot trajectory        
    trajectoryPlot = plot(solution, idxs = (1, 2, 3), title = "Satellite Trajectory", xlabel = "x", ylabel = "y", zlabel = "z", label = "Keplerian", size = (1000, 500))

    # Plot COEs
    angularMomentumPlot = plot(solution.t, angularMomentum, title = "Angular Momentum")
    inclinationPlot = plot(solution.t, inclination, title = "Inclination")
    raanPlot = plot(solution.t, raan, title = "RAAN Ω")
    eccentricityPlot = plot(solution.t, eccentricity, title = "Eccentricity")
    aopPlot = plot(solution.t, aop, title = "Argument of Perigee ω")
    trueAnomalyPlot = plot(solution.t, trueAnomaly, title = "True Anomaly θ")
#

display(trajectoryPlot)
display(plot(angularMomentumPlot, inclinationPlot, raanPlot, eccentricityPlot, aopPlot, trueAnomalyPlot, layout = (3, 2), size = (1000, 750)))