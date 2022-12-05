function drag(BC, ρ, v̲)
    """
    Calculates the perturbing acceleration from a drag force
    """

    v̲ᵣ = v̲ - v̲ₐ
    vᵣ = norm(v̲ᵣ)
    
    return -(1/2) * ρ * vᵣ * BC * v̲ᵣ
end