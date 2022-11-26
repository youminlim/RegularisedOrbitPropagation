# Main file to test numerical methods stored in files

# Import numerical method to file
cd(@__DIR__)
using Pkg
Pkg.activate(".")
include("rk1_4.jl")

# Define a test function
f(x, y) = x * sin(y)

# Initialise solver parameters
rk = 3      # RK solevr order
x0 = 0      # Initial x value
xf = 100    # Final x value
n = 101     # Number of steps
y0 = 0      # Initial y value

# Apply solver to system
x, y = rkSolver(rk, f, x0, xf, y0, n)

# Print solver values to terminal
for i in range(0, length = length(x))
    a = x[i]
    b = y[i]
    println("$a ,  $b")
end