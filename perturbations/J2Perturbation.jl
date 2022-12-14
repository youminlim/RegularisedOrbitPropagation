function J2coe(coe, centralBodyRadius, J₂)
    """
    Calculates the J2 Perturbation effect onto COEs of a body orbiting a planet with radius r
    Inputs:
        COE - Vector of COEs of the form [h, i, Ω, e, ω, θ]. Not all COEs are required but this is for simplicity
        centralBodyRadius - Radius of the central body to be orbited in km
        J₂ - Second zonal harmonic perturbation of the central body to be orbited
    Outputs:
        [Ω_dot, ω_dot] - rate of change in RAAN and AOP
    """

    # Unpack COE vector
    [h, i, Ω, e, ω, θ] = coe

    # Calculate miscellaneos COEs
    a = (h²/μ) * (1 / (1-e²)) 

    Ω_dot = -(3/2) * ((sqrt(mu)*J₂*centralBodyRadius^2) / (((1-e^2)^2) * a^(7/2))) * cosd(i)
    ω_dot = Ω_dot * ((5/2)*(sind(i)^2) - 2) / cosd(i)

    return [Ω_dot, ω_dot]
end

function J2state(x, y, z, centralBodyRadius = 6378, centralBodyJ₂ = 1.08262668e-3)
    """Calculates the state vector accelerations experiences by the satellite due to the J2 perturbation of the central body 
    Inputs:
        x, y, z - the state vector values expressing x, y and z dimensionless
        centralBodyRadius - Radius of the central body to be orbited in km
        J₂ - Second zonal harmonic perturbation of the central body to be orbited
    Outputs:
        [ẍ, ÿ, z̈] - acceleration in x, y and z directions
    """

    ẋJ2 = (1 - (5*z^2)/r^2) * x
    ẏJ2 = (1 - (5*z^2)/r^2) * y
    żJ2 = (3 - (5*z^2)/r^2) * z

    return (-(3/2)*(μ*centralBodyJ₂*centralBodyRadius^2) / (r^5)) * [ẋJ2, ẏJ2, żJ2]
end