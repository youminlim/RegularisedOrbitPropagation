function atmosphere(altitude)
    """
    Calculates the air density for a given altitude above Earth sea level
        Uses the simple exponential model
    Inputs:
        altitude - altitude above Earth sea level in km
    Outputs:
        ρ - density at specified altitude above Earth sea level
    """

    baseAltitude = [
    0, 25, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 
    180, 200, 250, 300, 350, 400, 460, 500, 600, 700, 800, 900, 1000
    ]

    nominalDensities = [
        1.225, 3.899e-2, 1.774e-2, 3.972e-3, 1.057e-3, 3.206e-4, 8.770e-5, 
        1.905e-5, 3.396e-6, 5.297e-7, 9.661e-8, 2.438e-8, 8.484e-9, 3.845e-9, 
        2.070e-9, 5.464e-10, 2.789e-10, 7.248e-11, 2.418e-11, 9.518e-12, 
        3.725e-12, 1.585e-12, 6.967e-13, 1.454e-13, 3.614e-14, 1.170e-14, 
        5.245e-15, 3.019e-15
    ]

    scaleHeight = [
        7.249, 6.349, 6.682, 7.554, 8.382, 7.714, 6.549, 5.799, 5.382, 5.877, 
        7.263, 9.473, 12.636, 16.149, 22.523, 29.740, 37.105, 45.546, 53.628, 
        53.298, 58.515, 60.828, 63.822, 71.835, 88.667, 124.64, 181.05, 268.00
    ]

    ρ = 0

    for i in baseAltitude
        if altitude < 0
            break
        elseif altitude < i
            j = findall(x -> x == i, baseAltitude)
            j = j[1] - 1
            ρ = nominalDensities[j] * exp(-(altitude - baseAltitude[j]) / scaleHeight[j])
            break
        end
    end

    return ρ
end