# Stores data for planets

function Earth()
    return Dict{String, Float64}(
        "mass" => 5.972e24,
        "radius" => 6379,
        "J2" => 1.08262668e-3
    )
end

function Moon()
    return Dict{String, Float64}(
        "mass" => 7.34767309e22,
        "radius" => 1737.2,
        "J2" => nothing
    )
end