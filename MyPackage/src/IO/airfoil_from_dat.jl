using ..Utils
using ..Geometry: Airfoil

"""
    airfoil_from_dat(datfile::String)

Loads an airfoil profile from a `.dat` coordinate file (Selig or Lednicer format), interpolates upper and lower surfaces with linear splines, and constructs an [`Airfoil`](@ref MyPackage.Geometry.Airfoil) object.

# Arguments

  - `datfile::String`: Path to the `.dat` file containing coordinate points ordered from trailing edge over the upper surface to leading edge and back across the lower surface to the trailing edge.

# Returns

  - `Airfoil`: Initialized airfoil object containing upper surface, lower surface, and mean camber line closures.

# Example

```julia
airfoil = airfoil_from_dat(\"assets/airfoils/Plain/Plain.dat\")
```
"""
function airfoil_from_dat(datfile::String)
    name = split(basename(datfile), ".")[1]
    println("Loading airfoil: \$name from \$datfile")
    airfoil_data = read_dat(datfile)
    top_surface, bottom_surface = calc_surfaces(airfoil_data)
    camber = calc_camber(top_surface, bottom_surface)
    return Airfoil(name, datfile, top_surface, bottom_surface, camber)
end

"""
    calc_surfaces(airfoil_data::AbstractMatrix{Float64})

Separates 2D airfoil coordinates into upper (suction) and lower (pressure) surfaces, enforcing leading-edge origin alignment at ``(0, 0)``, and returns linear splines for each.

# Returns

  - `(top_surface, bottom_surface)`: A tuple of splines mapping normalized chordwise station ``x \\in [0, 1]`` to ``z/c``.
"""
function calc_surfaces(airfoil_data)
    # find index of the leading-edge (minimum x)
    le_index = argmin(airfoil_data[:, 1])
    le_x = airfoil_data[le_index, 1]

    if le_x > 0.01
        # leading-edge not at x≈0 — prepend a (0,0) point
        top_data = vcat([(0.0, 0.0)], @views airfoil_data[le_index:-1:1, :])
        bottom_data = vcat([(0.0, 0.0)], @views airfoil_data[le_index:end, :])
    else
        top_data = @views airfoil_data[le_index:-1:1, :]
        bottom_data = @views airfoil_data[le_index:end, :]
    end

    # TODO: B-Spline top and bottom surfaces
    top_surface = LinearSpline(top_data[:, 1], top_data[:, 2])
    bottom_surface = LinearSpline(bottom_data[:, 1], bottom_data[:, 2])
    return top_surface, bottom_surface
end

"""
    calc_camber(top_surface, bottom_surface)
    calc_camber(datfile::String)

Computes the mean camber line function ``z_c(x) = \\frac{1}{2} [z_{\\text{top}}(x) + z_{\\text{bottom}}(x)]`` of an airfoil.
"""
function calc_camber(top_surface, bottom_surface)
    camber_func = x -> (top_surface(x) + bottom_surface(x)) / 2
    return camber_func
end

function calc_camber(datfile::String)
    airfoil_data = read_dat(datfile)
    top_surface, bottom_surface = calc_surfaces(airfoil_data)
    camber_func = x -> (top_surface(x) + bottom_surface(x)) / 2
    return camber_func
end
