function coe(state, μ)
    """
    This function will return the COE's given from the state vector provided
        package LinearAlgebra is REQUIRED

    Inputs:
        state - This will be the state vector of the body in the form (r, v) or [r, v], and must be 3 dimensional
        μ - This is the gravitational parameter, which changes with main orbiting body
    Outputs:
        [h, i, Ω, e, ω, θ] - This will be the vector of COE's returned 
    """

    # Unpack state vector
    r̲, v̲ = state[1:3], state[7:9]

    # Preliminary calculations
    r, v = norm(r̲), norm(v̲) # Distance and Speed
    vᵣ = (r̲ ⋅ v̲) / r         # Radial Velocity
    
    # Angular Momentum, h
    h̲ = cross(r̲, v̲)
    h = norm(h̲)

    # Inclination, i
    i = acosd(h̲[3] / h)

    # Node line, N
    N̲ = cross([0, 0, 1], h̲)
    N = norm(N̲)

    # Right ascension of the ascending node, Ω
    if N̲[2] >= 0
        Ω = acosd(N̲[1] / N)
    else 
        Ω = 360 - acosd(N̲[1] / N)
    end
    
        # Eccentricity, e
        e̲ = (1 / μ) * (((v^2 - (μ/r))*r̲) - (r * vᵣ * v̲))
        e = norm(e̲)
            # Checker for Eccentricity
            # if e != sqrt(1 + (((h/μ)^2) * ((v^2)- (2*μ/r))))
            #     println("Error in eccentricity")
            # end

        # Argument of Perigee, ω
        if e̲[3] >= 0
            ω = acosd((N̲ ⋅ e̲) / (N * e))
        else
            ω = 360 - acosd((N̲ ⋅ e̲) / (N * e))
        end

        # True anomaly, θ
        if vᵣ >= 0
            # θ = acosd((1/e)*(((h^2)/μ*r) - 1))
            θ = acosd((e̲ ⋅ r̲)/(e*r))
        else
            # θ = 360 - acosd((1/e)*(((h^2)/μ*r) - 1))
            θ = 360 - acosd((e̲ ⋅ r̲)/(e*r))
        end

    return [h, i, Ω, e, ω, θ]
end