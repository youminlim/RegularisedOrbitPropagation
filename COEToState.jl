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