function eulerian(a̲, b̲)
    """ Calculates the Eulerian distance between two points in arbitrary space - valid for 2D and 3D vectors
    input:
        a̲ - First vector
        b̲ - Second vector
        NOTE: Both vectors, a and b, have to be of the same dimensions
    Output:
        r - eulerian distance between two points
    """

    # Check dimensions of a and b match
    if length(a̲) == length(b̲)
        r̲ = a̲ - b̲
        r = norm(r)
    else
        println("Error: Vector dimensions do not match")
    end
    
    return r
end