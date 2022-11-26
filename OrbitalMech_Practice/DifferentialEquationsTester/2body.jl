# Testing DifferentialEquations package on Lorenz Attractor

cd(@__DIR__)
using Pkg
Pkg.activate(".")

using DifferentialEquations, BenchmarkTools, Plots

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

G = 6.67e-20
m1 = 1e26
m2 = 1e26

p = [G, m1, m2]
tspan = (0.0, 480.0)
u0 = [0; 0; 0; 3000; 0; 0; 10; 20; 30; 0; 40; 0]

# Original Function
    function f(u, p, t)
        r = eulerian(u[1:3], u[4:6])

        # Second derivatives
        d2x1 = G*m2*(u[4] - u[1])/r^3
        d2y1 = G*m2*(u[5] - u[2])/r^3
        d2z1 = G*m2*(u[6] - u[3])/r^3
        d2x2 = G*m1*(u[1] - u[4])/r^3
        d2y2 = G*m1*(u[2] - u[5])/r^3
        d2z2 = G*m1*(u[3] - u[6])/r^3

        [u[7], u[8], u[9], u[10], u[11], u[12], d2x1, d2y1, d2z1, d2x2, d2y2, d2z2]
    end

    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(prob, reltol=1e-8, saveat = 0.1)

    @benchmark solve(prob, reltol=1e-8, saveat=0.1)
#

# Function with Cache Array
    function f!(du, u, p, t)
        r = eulerian(u[1:3], u[4:6])

        # Build cache array
        du[1:6] = [u[7], u[8], u[9], u[10], u[11], u[12]]

        # Second derivatives
        du[7] = p[1]*p[3]*(u[4] - u[1])/r^3
        du[8] = p[1]*p[3]*(u[5] - u[2])/r^3
        du[9] = p[1]*p[3]*(u[6] - u[3])/r^3
        du[10] = p[1]*p[2]*(u[1] - u[4])/r^3
        du[11] = p[1]*p[2]*(u[2] - u[5])/r^3
        du[12] = p[1]*p[2]*(u[3] - u[6])/r^3
    end

    prob1 = ODEProblem(f!, u0, tspan, p)
    sol1 = solve(prob1, reltol=1e-8, saveat = 0.1)

    @benchmark solve(prob1, reltol=1e-8, saveat = 0.1)
#    

# Function with Cache Array and Precalculations
    function h!(du, u, p, t)
        r³ = eulerian(u[1:3], u[4:6])^3
        body1 = p[1] * p[3]
        body2 = p[1] * p[2]

        # Build cache array
        du[1:6] = u[7:12]

        # Second derivatives
        du[7] = body1*(u[4] - u[1])/r³
        du[8] = body1*(u[5] - u[2])/r³
        du[9] = body1*(u[6] - u[3])/r³
        du[10] = body2*(u[1] - u[4])/r³
        du[11] = body2*(u[2] - u[5])/r³
        du[12] = body2*(u[3] - u[6])/r³
    end

    prob2 = ODEProblem(h!, u0, tspan, p)
    sol2 = solve(prob2, reltol=1e-8, saveat = 0.1)

    @benchmark solve(prob2, reltol=1e-8, saveat = 0.1)
#    

# Function with Cache Array, Precalculations and Vectorisation
    function g!(du, u, p, t)
        r³ = eulerian(u[1:3], u[4:6])^3

        # Build cache array
        du[1:6] = u[7:12]

        # Second derivatives
        body1 = p[1] * p[3]
        body2 = p[1] * p[2]
        du[7:9] = body1 * (u[4:6] - u[1:3]) / r³
        du[10:12] = body2 * (u[1:3] - u[4:6]) / r³
    end

    prob3 = ODEProblem(g!, u0, tspan, p)
    sol3 = solve(prob3, reltol=1e-8, saveat = 0.1)

    @benchmark solve(prob3, reltol=1e-8, saveat = 0.1)
#


# # Plot relative position
# plot(sol1, idxs = (1,2,3), title = "2 Body Problem", label = "Mass 1")
# plot!(sol1, idxs = (4,5,6), label = "Mass 2", show = true)

# # Plot velocities
# plot(sol1, idxs = (7,8,9), title = "velocities", label = "Mass 1")
# plot!(sol1, idxs = (10,11,12), label = "Mass 2", show = true)

x1f = [last(sol1.u)[1]]
y1f = [last(sol1.u)[2]]
z1f = [last(sol1.u)[3]]
x2f = [last(sol1.u)[4]]
y2f = [last(sol1.u)[5]]
z2f = [last(sol1.u)[6]]

p1 = plot(sol1, idxs = (1, 2, 3), colour = :blue, title = "Position", label = "Mass 1", xlabel = "x", ylabel = "y", zlabel = "z")
plot!(p1, sol1, idxs = (4, 5, 6), colour = :red, label = "Mass 2", cam = (30,30), xlims = (0, 6e6))
plot!(p1, x1f, y1f, z1f, colour = :blue, seriestype = :scatter, label = nothing, lw = 10)
plot!(p1, x2f, y2f, z2f, colour = :red, seriestype = :scatter, label = nothing, lw = 10)
p2 = plot(sol1, idxs = (7, 8, 9), colour = :blue, title = "Velocity", label = "Mass 1")
plot!(p2, sol1, idxs = (10, 11, 12), colour = :red, label = "Mass 2", xlims = (-6.5e4, 7.5e4))
display(plot(p1, p2, layout = (1, 2), size = (1000,500)))

# savefig("plot.png")