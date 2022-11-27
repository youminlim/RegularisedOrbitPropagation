function TwoBody3D!(du, u, p, t)
    """
    Function definition of the 2 Body Problem\
    Input:
        u - Initial state vector
        p - Parameters to pass through for calculation, re
            REQUIRED: [G, M₁, M₂]:
                G - Universal Gravitational Parameter
                M₁ - Mass of body 1
                M₂ - Mass of body 2
        t - Time variable

    Output:
        du - Cache array which stores the state vector values after an iteration
    """

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

# Define 2 body function with varargs
# function twobody3d(du, u::Vector{Float64}, p::Vector{Float64, N}, N::Int64 == 1, Bodies::Vararg{String, :})
#     r_cubed = eulerian(u[1:3], u[4:6])^3
#     body1 = p[1] * p[3]
#     body2 = p[1] * p[2]

#     # Build cache array
#     du[1:6] = u[7:12]

#     # Second derivatives
#     du[7] = body1*(u[4] - u[1])/r_cubed
#     du[8] = body1*(u[5] - u[2])/r_cubed
#     du[9] = body1*(u[6] - u[3])/r_cubed
#     du[10] = body2*(u[1] - u[4])/r_cubed
#     du[11] = body2*(u[2] - u[5])/r_cubed
#     du[12] = body2*(u[3] - u[6])/r_cubed
# end