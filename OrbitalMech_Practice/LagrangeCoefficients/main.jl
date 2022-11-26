using Base: stacktrace_contract_userdir
# Calculate r and v after Δθ from r̲₀ and v̲₀ using Lagrange Coefficients

# Import Packages
cd(@__DIR__)
using Pkg
Pkg.activate(".")

using LinearAlgebra
using Plots

# Compute Lagrange Coefficients
function lagrangeCoefficients(r̲₀, v̲₀, Δθ)
    # Convert Δθ
    Δθ *= π/180
    s = sin(Δθ)
    c = cos(Δθ)

    # Magnitude of r̲₀ and v̲₀
    r₀ = sqrt(dot(r̲₀, r̲₀))
    v₀ = sqrt(dot(v̲₀, v̲₀))

    # Radial component of v̲₀
    vᵣ₀ = dot(r̲₀, v̲₀) / r₀

    # Magnitude of constant angular momentum
    h = r₀ * sqrt(v₀^2 - vᵣ₀^2)

    # Calculate r
    r = (h^2 / μ) * (1 / (1 + ((h^2/(μ*r₀))-1)*cos(Δθ) - (h*vᵣ₀/μ)*s))

    # Calculate lagrange Coefficients
    f = 1 - ((μ*r/h^2)*(1 - c))
    g = (r*r₀/h) * s
    ḟ = (μ/h) * ((1 - c)/s) * ((μ/h^2)*(1 - c) - 1/r₀ - 1/r)
    ġ = 1 - (μ*r₀/h^2)*(1 - c)

    # Check if lagrange coefficients satisfy conservation of momentum
    condition = (f*ġ) - (ḟ*g)
    if ( 0.99999 > condition) && (condition > 1.0001)
        println("Error!")
        println(r̲₀)
    end

    # Calculate new position and velocity vectors
    r̲ = (f * r̲₀) + (g * v̲₀)
    v̲ = (ḟ * r̲₀) + (ġ * v̲₀)

    return r̲, v̲
end

# Global Conditions
G = 6.67e-20
M = 5.972e24
R = 6371e3
μ = G * M

# Initial Conditions
r̲₀ = [8182.4; -6865.9]
v̲₀ = [0.47572; 8.8116]

Δθ = 0.01
θ = range(0, 6*pi, step = Δθ)
r̲ = zeros(2, length(θ))
v̲ = zeros(2, length(θ))
r̲[:, 1] = r̲₀
v̲[:, 1] = v̲₀

# Calculate r̲ and v̲ for range of θ
for i in range(1, length(θ)-1)
    r̲[:, i+1], v̲[:, i+1] = lagrangeCoefficients(r̲[:, i], v̲[:, i], Δθ)
end

# Plot Results
    a = size(r̲)

    x = r̲[1, :]
    y = r̲[2, :]
    plot(x, y)

    savefig("plot.png")