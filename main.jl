""" Numerical Propagation for 2 Body Problem with Cowell's Method """

# Add Packages and Libraries
cd(@__DIR__)
using Pkg
Pkg.activate(".")

include("Planets.jl")
include("2Body3D.jl")
include("StateToCOE.jl")

using LinearAlgebra
using DifferentialEquations
using Plots

# Dictionary to store defaults perturbation settings
function defaultPerturbations()
    return Dict{String, Bool}(
                "J2" => false,
                "Drag" => false,
                "DragModel" => false,
                "SolarRadiationPressure" => false,
                "ThirdBody" => false,
                )
end

function defaultParameters()
    include("Planets.jl")
    return Dict(
        "satelliteMass" => 1000,
        "centralBody" => Earth(),
        "perturbations" => defaultPerturbations()
    )
end

function parameters(p = defaultParameters())
    # unpack parameters to use
    satelliteMass = p["satelliteMass"]
    centralBody = p["centralBody"]
    perturbations = p["perturbations"]

    # organise parameters
    input = [satelliteMass, centralBody["mass"], centralBody["radius"], centralBody["J2"]]
    perturbationKeys = keys(perturbations)
    selectedPerturbations = []
    for i in perturbationKeys
        if perturbations[i] == true
            push!(selectedPerturbations, i)
        end
    end
    input = [input; selectedPerturbations]

    return input
end

# 2 Body problem with Cowell's Method for a satellite orbit about a central body
function f!(du, u, p, t)
    """
    Inputs:
        state vector u - 6 states, [x₁, y₁, z₁, x₂, y₂, z₂, ẋ₁, ẏ₁, ż₁, ẋ₂, ẏ₂, ż₂], as motion of central body is not concerned
        p - parameters to pass to function, has structure : [satelliteMass, centralBody, perturbations]
            satelliteMass - contains a value representing the mass of the satellite

            centralBody - will have structure [mass, radius, J2] of the central body

            perturbations - will have structure [nothing], override by adding in the string name of each perturbation type, i.e. "J2" in order to apply the accelerations from J2 to the system
    Output:
        du - the derivative of the state vector i.e. u̇
    """

    # Housekeeping
    satelliteMass = p[1]
    centralBody = p[2:4]
    perturbations = p[4:length(p)]

    x, y, z = u[1:3]
    r̲ = [x, y, z]
    r = norm(r̲)
    du[1:3] = u[4:6]

    G = 6.67e-20
    μ = centralBody[1] * G

    # Calculate accelerations
    du[4:6] = - μ * r̲ * satelliteMass / r^3

    # Add perturbation effects
    if "J2" in perturbations
        # Housekeeping
        R = centralBody[2]
        J₂ = centralBody[3]
        ẋJ2 = (1 - (5*z^2)/r^2) * x
        ẏJ2 = (1 - (5*z^2)/r^2) * y
        żJ2 = (3 - (5*z^2)/r^2) * z
        a_J2 = (-(3/2)*(μ*J₂*R^2) / (r^5)) * [ẋJ2, ẏJ2, żJ2]

        du[4:6] += a_J2
    end

    if "Drag" in perturbations
        nothing
    end

end

# Define global constants
G = 6.67e-20
centralBody = Earth()
μ = G * centralBody["mass"]

# Define state vectors
u̲₀ = [8000.0, 6000.0, 0.0, 0.0, 5.0, 5.0]
tspan = (0.0, 1440.0)

# solve keplerian Trajectory
perturbations = defaultPerturbations()
p = defaultParameters()
input = parameters(p)

problem = ODEProblem(f!, u̲₀, tspan, input)
solution = solve(problem, dt = 0.1)

# Convert state vector data to COEs
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
    #aop = aop/maximum(aop)
    trueAnomaly = getindex.(solutionCOE, 6)
    #trueAnomaly = trueAnomaly/maximum(trueAnomaly)
#

# Plot data
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

# Set Perturbations
# perturbations["J2"] = true

# # Define parameter
# p["centralBody"] = Earth()
# p["perturbations"] = perturbations
# inputPerturbed  = parameters(p)

# problemPerturbed = ODEProblem(f!, u̲₀, tspan, inputPerturbed)
# solutionPerturbed = solve(problemPerturbed, reltol=1e-9, saveat = 0.1)

# # Plot data
#     plot!(trajectoryPlot, solutionPerturbed, idxs = (1, 2, 3), label = "J2 Perturbation")
# #
    
# Display Plots
    display(trajectoryPlot)
    display(plot(angularMomentumPlot, inclinationPlot, raanPlot, eccentricityPlot, aopPlot, trueAnomalyPlot, layout = (3, 2), size = (1000, 750)))