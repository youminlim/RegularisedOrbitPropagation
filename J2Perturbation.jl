function J2(coe, R, J₂)
    """
    Calculates the J2 Perturbation effect onto COEs of a body orbiting a planet with radius r
    Inputs:
        COE - Vector of COEs of the form [h, i, Ω, e, ω, θ]. Not all COEs are required but this is for simplicity
        R - Radius of the planet to be orbited in km
        J₂ - Second zonal harmonic perturbation of the planet to be orbited
    Outputs:
        not really sure yet
    """

    # Unpack COE vector
    [h, i, Ω, e, ω, θ] = coe

    # Calculate miscellaneos COEs
    a = (h²/μ) * (1 / (1-e²)) 

    Ω_dot = -(3/2) * ((sqrt(mu)*J₂*R^2) / (((1-e^2)^2) * a^(7/2))) * cosd(i)
    ω_dot = Ω_dot * ((5/2)*(sind(i)^2) - 2) / cosd(i)

    return [Ω_dot, ω_dot]
end