"""This file will collect any functions used for coordinate transfers"""

# Add in pre-requisite files
using LinearAlgebra

function coe(state, μ::Float64 = 398332.4)
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
            θ = acosd((e̲ ⋅ r̲)/(e*r))
        else
            θ = 360 - acosd((e̲ ⋅ r̲)/(e*r))
        end

    return [h, i, Ω, e, ω, θ]
end

function coeCB(state, μ::Float64 = 398332.4)
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
    r̲, v̲ = state[4:6], state[10:12]

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

# Converts a single state vector to a line of COEs
function state(coe, μ::Float64 = 398332.4)
    """
    Converts COEs into state vector form
    Inputs:
        coe - vector of COEs in the form [h, i, Ω, e, ω, θ]
        μ - Gravtiational Parameter
    Outputs:
        [x y z ẋ ẏ ż] - vector of states in cartesian coordinates
    """
    
    # Unpack COE vector
    h, i, Ω, e, ω, θ = coe

    # Calculate miscellaneos COEs

    # Calculate r̲₀ and v̲₀
    r̲₀ = ((h^2)/μ) * (1 / (1 + (e * cosd(θ)))) .* [cosd(θ); sind(θ); 0]
    v̲₀ = (μ/h) .* [-sind(θ); e + cosd(θ); 0]

    # Calculate Rotation Matrix Q̲
    Q11 = -sind(Ω)*cosd(i)*sind(ω) + cosd(Ω)*cosd(ω)
    Q12 = -sind(Ω)*cosd(i)*cosd(ω) - cosd(Ω)*sind(ω)
    Q13 = sind(Ω)*sind(i)
    Q21 = cosd(Ω)*cosd(i)*sind(ω) + sind(Ω)*cosd(ω)
    Q22 = cosd(Ω)*cosd(i)*cosd(ω) - sind(Ω)*sind(ω)
    Q23 = -cosd(Ω)*sind(i)
    Q31 = sind(i)*sind(ω)
    Q32 = sind(i)*cosd(ω)
    Q33 = cosd(i)
    Q̲ = [Q11 Q12 Q13; Q21 Q22 Q23; Q31 Q32 Q33]

    # Calculate r̲ₓ and v̲ₓ
    r̲ₓ = Q̲ * r̲₀
    v̲ₓ = Q̲ * v̲₀

    return [r̲ₓ; v̲ₓ]
end

# Extract state vectors from solution of solve() to convert to COEs
function state2coe(solution, μ)
    solutionCOE = []

    for i in solution.u
        push!(solutionCOE, coe(i, μ))
    end

    return solutionCOE
end

# Extract state vectors from solution of solve() to convert to COEs for the central body
function state2coeCB(solutionCB, μCB)
    solutionCBCOE = []

    for i in solutionCB.u
        push!(solutionCBCOE, coeCB(i, μCB))        
    end

    return solutionCBCOE
end

# Extract each COE from a vector of COEs
function coeExtractor(solutionCOE)
    angularMomentum = getindex.(solutionCOE, 1)
    inclination = getindex.(solutionCOE, 2)
    raan = getindex.(solutionCOE, 3)
    eccentricity = getindex.(solutionCOE, 4)
    aop = getindex.(solutionCOE, 5)
    trueAnomaly = getindex.(solutionCOE, 6)

    return [angularMomentum, inclination, raan, eccentricity, aop, trueAnomaly]
end

# Test values
# e = 0.0004094
# i = 51.6417 # degrees
# Ω = 183.0231 # degrees
# ω = 142.6476 # degrees
# θ = 100 # degrees
# apogee = Rₑ + 420.0 # km
# perigee = Rₑ + 414.0 # km
# a = (apogee + perigee) / 2
# b = a * sqrt(1 - e^2)
# frequencyPerDay = 15.49860975 
# orbitalPeriod = 86400 / frequencyPerDay # seconds
# h = ((2 * π * a * b) / orbitalPeriod) * 10^-6
# satellite = state([h, i, Ω, ω, e, θ])