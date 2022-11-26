# 2 Body Problem Solver

# Add Packages and Libraries
cd(@__DIR__)
using Pkg
Pkg.activate(".")

using Plots
include("myFunctions.jl")

# System parameters
G = 6.67e-11
m1 = 10e26
m2 = 10e26
#y0 = [0, 0, 0, 3000e3, 0, 0, 10e3, 20e3, 30e3, 0, 40e3, 0] # In the form y0 = [x01 y01 z01 x02 y02 z02 dot(x01 y01 z01 x02 y02 z02)]
y0 = [0, 0, 3000e3, 0, 10e3, 20e3, 0, 40e3]

# Time interval
t0 = 0      # Initial time
tf = 60   # Final time
n = 1000

t, y = rkSolver(twoBody2D, t0, tf, y0, n)

plot2Body(t, y)

savefig("plot.png")

#accelChecker(G, m1, m2, y)