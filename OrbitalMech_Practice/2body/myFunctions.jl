# Functions to solve and visualise 2 Body Problem

function twoBody3D(t, u, r)
    r = eulerian(u[1:3], u[4:6])

    # Second derivatives
    d2x1 = G*m2*(u[4] - u[1])/r^3
    d2y1 = G*m2*(u[5] - u[2])/r^3
    d2z1 = G*m2*(u[6] - u[3])/r^3
    d2x2 = G*m2*(u[1] - u[4])/r^3
    d2y2 = G*m2*(u[2] - u[5])/r^3
    d2z2 = G*m2*(u[3] - u[6])/r^3

    [u[7], u[8], u[9], u[10], u[11], u[12], d2x1, d2y1, d2z1, d2x2, d2y2, d2z2]
end

function twoBody2D(t, u, r)
    r = eulerian(u[1:2], u[3:4])

    # Second derivatives
    d2x1 = G*m2*(u[3] - u[1])/r^3
    d2y1 = G*m2*(u[4] - u[2])/r^3
    d2x2 = G*m2*(u[1] - u[3])/r^3
    d2y2 = G*m2*(u[2] - u[4])/r^3

    [u[5], u[6], u[7], u[8], d2x1, d2y1, d2x2, d2y2]
end


function eulerian(vector1, vector2)
    # Calculates the Eulerian distance between two vectors

    distance = vector2 - vector1
    if length(distance) == 1
        return distance
    else
        squares = 0
        for i in distance
            squares += i^2
        end
        return sqrt(squares)
    end
end


function plot2Body(t, y)
    a = size(y)

    if a[2] == 8 
        body1 = y[:, 1:2]
        body2 = y[:, 3:4]
        plot(body1[:, 1], body1[:, 2], t, title = "2 Body Problem", aspect_ratio = :equal, label = "Mass 1")
        plot!(body2[:, 1], body2[:, 2], t, label = "Mass 2")

    elseif a[2] == 12
        body1 = y[:, 1:3]
        body2 = y[:, 4:6]
        plot(body1[:, 1], body1[:, 2], body1[:, 3], title = "2 Body Problem", aspect_ratio = :equal, label = "Mass 1")
        plot!(body2[:, 1], body2[:, 2], body2[:, 3], label = "Mass 2")
    end
end


function rkSolver(f, t0, tf, y0, n)
    # Declare arrays
    h = (tf - t0) / n
    t = range(t0, tf, n)
    y = zeros(n, length(y0))

    # Add in initial conditions
    y[1, :] = y0

    if length(y0) == 2
        r = eulerian(y[1, 1:2], y[1, 3:4])

        for i in range(1, n-1)
            k1 = f(t[i], y[i, :], r)
            k2 = f(t[i] + (0.5*h), y[i, :] + (0.5*h*k1), r)
            k3 = f(t[i] + (0.5*h), y[i, :] + (0.5*h*k2), r)
            k4 = f(t[i] + h, y[i, :] + (h*k3), r)
            y[i+1, :] = y[i, :] + (h * (k1/6 + k2/3 + k3/3 + k4/6))
            r = eulerian(y[i+1, 1:2], y[i+1, 3:4])
        end
    else
        r = eulerian(y[1, 1:3], y[1, 4:6])

        for i in range(1, n-1)
            k1 = f(t[i], y[i, :], r)
            k2 = f(t[i] + (0.5*h), y[i, :] + (0.5*h*k1), r)
            k3 = f(t[i] + (0.5*h), y[i, :] + (0.5*h*k2), r)
            k4 = f(t[i] + h, y[i, :] + (h*k3), r)
            y[i+1, :] = y[i, :] + (h * (k1/6 + k2/3 + k3/3 + k4/6))
            r = eulerian(y[i+1, 1:2], y[i+1, 3:4])
        end
    end

    return t, y
end


# Physics Checker
function gravity(G, m, y)
    # Calculates the gravitational force experienced by a body

    r = eulerian(y[1:3], y[4:6])
    return (G * m) / r^2
end


function accelChecker(G, m1, m2, y)
    # Checks if the acceleration forces experienced by the bodies are of correct magnitude
    
    dims = size(y)
    for i in range(1, length = dims[1])
        if gravity(G, m1, y[i, :]) != -(m1 / m2) * gravity(G, m2, y[i, :])
         println("Error!")
        end 
     end 
end