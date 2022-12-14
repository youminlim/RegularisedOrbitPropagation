""" This function calculates the accelerations for the 2 body problem """

# Add Pre-requisite functions for the solver
include("DragPerturbation.jl")

function twoBodySolver!(du, u, p, t)
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
    μ₁ = (G * centralBody[1])
    μ₂ = (satelliteMass * G)

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