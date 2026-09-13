"""
    read_dat(filename::String; header_lines::Int=1)

Read a `.dat` file containing airfoil coordinates and return the points as a 2D array.

# Arguments
- `filename::String`: Path to the `.dat` file.
- `header_lines::Int`: Number of header lines to skip (default: `1`).

# Returns
- `Array{Float64, 2}`: A `(N, 2)` matrix of `(x, z)` airfoil coordinates, ordered from
  trailing edge over the upper surface through the leading edge and back to the trailing edge.

# Example
```julia
points = read_dat("path/to/airfoil.dat")
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
