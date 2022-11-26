# Applies Runge Kutta method to the system

function rkSolver(rk, f, t0, tf, y0, n)
    # Declare arrays
    h = (tf - t0) / n
    t = range(t0, tf, n)
    y = zeros(n, length(y0))

    # Add in initial conditions
    y[1, :] = y0

    # RK order checker
    if rk == 1
        for i in range(1, n-1)
            k1 = f(t[i], y[i, :], r)
            y[i+1, :] = y[i, :] + (h * (k1/6))
        end
    elseif rk == 2
        for i in range(1, n-1)
            k1 = f(t[i], y[i, :])
            k2 = f(t[i] + h, y[i, :] + (h*k1))
            y[i+1, :] = y[i, :] + (h * (k1/2 + k2/2))
        end
    elseif rk == 3
        for i in range(1, n-1)
            k1 = f(t[i], y[i, :])
            k2 = f(t[i] + (0.5*h), y[i, :] + (0.5*h*k1))
            k3 = f(t[i] + h, y[i, :] + (h*(-k1+2*k2)))
            y[i+1, :] = y[i, :] + (h * (k1/6 + k2/3 + k3/6))
        end
    elseif rk == 4
        for i in range(1, n-1)
            k1 = f(t[i], y[i, :])
            k2 = f(t[i] + (0.5*h), y[i, :] + (0.5*h*k1))
            k3 = f(t[i] + (0.5*h), y[i, :] + (0.5*h*k2))
            k4 = f(t[i] + h, y[i, :] + (h*k3))
            y[i+1, :] = y[i, :] + (h * (k1/6 + k2/3 + k3/3 + k4/6))
        end
    else
        printIn("Order of Runge-Kutta method not defined")
    end

    return t, y
end