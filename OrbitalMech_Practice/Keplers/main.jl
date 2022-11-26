# Solver for Kepler's Equation

# Activate packages in environment
cd(@__DIR__)
using Pkg
Pkg.activate(".")

include("myFunctions.jl")

using Plots

# Define global variables
altitude = 2000 # At perigee
r_earth = 6379
G = 6.67*10^-11
M = 5.972*10^24

# Orbit parameters
e = 0.9
a = (altitude + r_earth) / (1 - e)
# i = 0
# Ω = 0
# ω = 0
M0 = 0

#initial = [e, a, i, Ω, ω, M0]
M0 = deg2rad(M0)

### Program 

mu = G * M * 10^-9
period = 2 * pi * sqrt(a^3 / mu)
tmin = 0
tmax = period
time = [tmin, tmax]
dt = LinRange(tmin, tmax, 1000)

# COE's
E = zeros(length(dt));
x = zeros(length(dt));
y = zeros(length(dt));
n = sqrt(mu / a^3);

for i in range(1, length(dt))
    M = M0 + (n * dt[i]);
    E[i] = solver(M, e);
    x[i] = a * (cos(E[i]) - e);
    y[i] = a * sin(E[i]) * sqrt(1 - e^2);
end

plotter(x,y)
# Save plot figure to folder
# savefig("plot.png")