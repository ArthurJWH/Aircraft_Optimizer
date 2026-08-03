"""
    read_dat(filename::String; header_lines::Int=1)

    Reads a .dat file containing airfoil coordinates and returns the points as a 2D array.

    Arguments
    ---------
    filename : String
        The path to the .dat file.
    header_lines : Int, optional
        The number of header lines to skip (default is 1).

    Returns
    -------
    points : Array{Float64, 2}
        A 2D array containing the airfoil coordinates.

    Example
    -------
    ```julia
    points = read_dat("path/to/airfoil.dat"; header_lines=1)
    ```
"""
function read_dat(filename::String; header_lines::Int=1)
    lines = readlines(filename)
    len = length(lines) - header_lines

    points = Array{Float64, 2}(undef, len, 2)

    for i in 1:len
        line = lines[header_lines + i]
        coords = split(line)
        points[i, 1] = parse(Float64, coords[1])
        points[i, 2] = parse(Float64, coords[2])
    end

    return points
end
