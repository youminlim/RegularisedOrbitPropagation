using Base: MainInclude, ForwardOrdering
""" Numerical Propagation for 2 Body Problem with Cowell's Method """

# Add Packages and Libraries
cd(@__DIR__)
using Pkg
Pkg.activate(".")

include("Planets.jl")
include("2Body3D.jl")
include("StateToCOE.jl")
include("DragPerturbation.jl")
include("atmosphere.jl")

using LinearAlgebra
using DifferentialEquations
using GeometryBasics
using GLMakie

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

function TwoBodyInitial(u)
    """
    Rewrites the state vector of a satellite to include the state vector of the central body as input to the 2 body solver
    """
    
    return [u[1:3]; 0; 0; 0; u[4:6]; 0; 0; 0] 
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
    du[1:6] = u[7:12]

    G = 6.67e-20
    μ₁ = centralBody[1] * G
    μ₂ = satelliteMass * G

    # Central Body

    # Calculate accelerations
    du[7:9] = - μ₁ * r̲ / r^3
    du[10:12] = μ₂ * r̲ / r^3

    # Add perturbation effects
    if "J2" in perturbations
        # Housekeeping
        R = centralBody[2]
        J₂ = centralBody[3]
        ẋJ2 = (1 - (5*z^2)/r^2) * x
        ẏJ2 = (1 - (5*z^2)/r^2) * y
        żJ2 = (3 - (5*z^2)/r^2) * z
        a_J2 = (-(3/2)*(μ*J₂*R^2) / (r^5)) * [ẋJ2, ẏJ2, żJ2]

        du[7:9] += a_J2
    end

    if "Drag" in perturbations
        BC = 1 / 82
        a_Drag = drag(BC, u)

        du[7:9] += a_Drag
    end

end

# Extract state vectors from solution of solve() to convert to COEs
function state2coe(solution, μ)
    solutionCOE = []

    for i in solution.u
        push!(solutionCOE, coe(i, μ))
    end

    return solutionCOE
end

# Extract each COE
function coeExtractor(solutionCOE)
    angularMomentum = getindex.(solutionCOE, 1)
    inclination = getindex.(solutionCOE, 2)
    raan = getindex.(solutionCOE, 3)
    eccentricity = getindex.(solutionCOE, 4)
    aop = getindex.(solutionCOE, 5)
    trueAnomaly = getindex.(solutionCOE, 6)

    return [angularMomentum, inclination, raan, eccentricity, aop, trueAnomaly]
end



# Define global constants
G = 6.67e-20
centralBody = Earth()
μ = G * centralBody["mass"]

# Define state vectors
u̲₀ = [8000.0, 0.0, 6000.0, 0.0, 5.0, 5.0]
u̲₀ = TwoBodyInitial(u̲₀)
tspan = (0.0, 160000.0)

# solve keplerian Trajectory
perturbations = defaultPerturbations()
p = defaultParameters()
input = parameters(p)

problem = ODEProblem(f!, u̲₀, tspan, input)
solution = solve(problem, RK4())

# Solve perturbation trajectory
p["perturbations"]["J2"] = true
p["perturbations"]["Drag"] = true
inputPerturbed = parameters(p)

problemPerturbed = ODEProblem(f!, u̲₀, tspan, inputPerturbed)
solutionPerturbed = solve(problemPerturbed, RK4())

# Convert state vector data to COEs
begin
    coes = coeExtractor(state2coe(solution, μ))
    coesPerturbed = coeExtractor(state2coe(solutionPerturbed, μ))

    time = solution.t / 3600 
    timePerturbed = solutionPerturbed.t / 3600
end
#

# Plot data
begin
    # Plot trajectory
    fig = Figure()
    ax = Axis3(
        fig[1, 1];
        title = "Satellite Trajectory"
    )

    x = solution[1,:]
    y = solution[2,:]
    z = solution[3,:]

    xₚ = solutionPerturbed[1,:]
    yₚ = solutionPerturbed[2,:]
    zₚ = solutionPerturbed[3,:]

    lines!(ax, x, y, z, label = "Keplerian Trajectory", color = "red")
    lines!(ax, xₚ, yₚ, zₚ, label = "Perturbed Trajectory", color = "green")
    mesh!(ax, Sphere(Point3f0(0), 6371))

    #trajectoryPlot = plot(solution, idxs = (1, 2, 3), title = "Satellite Trajectory", xlabel = "x", ylabel = "y", zlabel = "z", label = "Keplerian trajectory", size = (1000, 500))
    #plot!(trajectoryPlot, solutionPerturbed, idxs = (1, 2, 3), label = "Pertubed trajectory")
    fig[1, 2] = Legend(fig, ax, "Legend")
    fig
end

# Plot COEs
function CoePlotter(COEs, time, additionalCOEs::Vector{Vector{Float64}})
    # Plot COEs
    angularMomentumPlot = plot(time, COEs[1], title = "Angular Momentum")
    inclinationPlot = plot(time, COEs[2], title = "Inclination")
    raanPlot = plot(time, COEs[3], title = "RAAN Ω")
    eccentricityPlot = plot(time, COEs[4], title = "Eccentricity")
    aopPlot = plot(time, COEs[5], title = "Argument of Perigee ω")
    trueAnomalyPlot = plot(time, COEs[6], title = "True Anomaly θ")

    # Plot names 
    plotNames = [angularMomentumPlot, inclinationPlot, raanPlot, eccentricityPlot, aopPlot, trueAnomalyPlot]

    # Plot perturbed COEs
    for i in range(1, length(plotNames))
        plot!(plotNames[i], time, additionalCOEs[i])
    end

    display(plot(angularMomentumPlot, inclinationPlot, raanPlot, eccentricityPlot, aopPlot, trueAnomalyPlot, layout = (3, 2), size = (1250, 1150)))
end
 
# Setup Plot
begin
    # Create Figure
    fig = Figure(resolution = (1275, 1200))

    # Create Grids
        ga = fig[1:2, 1] = GridLayout()
        # gb = fig[2, 1] = GridLayout()
        gc = fig[1, 2] = GridLayout()
        gd = fig[2, 2] = GridLayout()

    # Labels
    labels = ["Keplerian", "Perturbed"]

    # For ga
        # Create Axis
        ax1 = Axis(ga[1, 1], title = "Angular Momentum", ylabel = "h in kgm²/s")
        ax2 = Axis(ga[2, 1], title = "Inclination", ylabel = "i in °")
        ax3 = Axis(ga[3, 1], title = "RAAN Ω", ylabel = "Ω in °")
        ax4 = Axis(ga[4, 1], title = "Eccentricity", ylabel = "e, dimensionless")
        ax5 = Axis(ga[5, 1], title = "AOP ω", ylabel = "ω in °")
        ax6 = Axis(ga[6, 1], title = "True Anomaly θ", xlabel = "Hours", ylabel = "θ in °")
        axis = [ax1, ax2, ax3, ax4, ax5, ax6]

        # Fill axes with data
        for i in range(1, length(axis))
            lines!(axis[i], time, coes[i], label = labels[1], color = "red")
            lines!(axis[i], timePerturbed, coesPerturbed[i], label = labels[2], color = "green")
        end
        for i in axis[1:5]
            hidexdecorations!(i, grid = false)
        end
    #
    # For gc
        # Create Axis
        axc = Axis3(gc[1, 1], title = "Satellite Trajectory")

        # Prepare data
        x = solution[1,:]
        y = solution[2,:]
        z = solution[3,:]

        xₚ = solutionPerturbed[1,:]
        yₚ = solutionPerturbed[2,:]
        zₚ = solutionPerturbed[3,:]

        # Fill axes with data
        lines!(axc, x, y, z, label = "Keplerian Trajectory", color = "red")
        lines!(axc, xₚ, yₚ, zₚ, label = "Perturbed Trajectory", color = "green")
        mesh!(axc, Sphere(Point3f0(0), 6371))
    #
    # For gd 
        # Create Axis
        ax = Axis(gd[1, 1], title = "x distance")
        ay = Axis(gd[2, 1], title = "y distance")
        az = Axis(gd[3, 1], title = "z distance")
        axisDistance = [ax, ay, az]

        # Prepare data
        keplerian = [x, y, z]
        perturbed = [xₚ, yₚ, zₚ]

        # Fill axes with data
        for i in range(1, length(axisDistance))
            lines!(axisDistance[i], time, keplerian[i], color = "red")
            lines!(axisDistance[i], timePerturbed, perturbed[i], color = "green")
        end
        for i in axisDistance[1:2]
            hidexdecorations!(i, grid = false)
        end
    #

    # leg = Legend(gd[1, 1], ax1)
    # leg.tellheight = true

    colsize!(fig.layout, 1, Auto(0.6))
    colsize!(gc, 1, Auto())
    fig
end