function drag(BC, u)
    """
    Calculates the perturbing acceleration from a drag force
    """
    
    r̲ = u[1:3]
    ρ = atmosphere(norm(r̲) - 6371)
    v̲ᵣ = u[7:9]
    vᵣ = norm(v̲ᵣ)

    return -(1/2) * ρ * vᵣ * BC * v̲ᵣ
end