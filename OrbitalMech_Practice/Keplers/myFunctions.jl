function solver(M, e)
    E = M
    error = 1

    while error >= 0.0001
        f = E - (e * sin(E)) - M
        df = 1 - (e * cos(E))

        # Newton-Raphson formula
        E_updated = E - (f / df)

        # Update error and eccentric anomaly
        error = E_updated - E
        E = E_updated
        
        return E
    end
end

function planet_radius(planetName)
    # Returns the planetName's radius in km
    if planetName == "earth"
        return 6379
    elseif planetName == "sun"
        return 696340
    end
end

function planetary_body_2Dpoints(r)
    # r - radius of the planetary body to plot
    theta = LinRange(0, 2*pi, 50)
    return r*cos.(theta), r*sin.(theta)
end

function plotter(x, y; planetName = "earth")

    plot(planetary_body_2Dpoints(planet_radius(planetName)),  title = "Satellite Trajectory", label = "Earth", aspect_ratio = :equal)
    plot!(x, y, label = "Satellite")
end