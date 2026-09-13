"""
    Airfoil{T, B, C}

Represents a 2D airfoil geometry with top/bottom surface splines and a mean camber line closure.

# Fields
- `name::SubString{String}`: Name of the airfoil, extracted from the file basename.
- `datfile::String`: Path to the source `.dat` coordinate file.
- `top_surface::T`: Function/spline evaluating the upper (suction) surface coordinates ``z(x)`` for ``x \\in [0, 1]``.
- `bottom_surface::B`: Function/spline evaluating the lower (pressure) surface coordinates ``z(x)`` for ``x \\in [0, 1]``.
- `camber::C`: Function evaluating the mean camber line ``z_c(x) = \\frac{1}{2}[z_{\\text{top}}(x) + z_{\\text{bot}}(x)]``.

Construct via [`MyPackage.IO.airfoil_from_dat`](@ref).
"""
struct Airfoil{T, B, C}
    name::SubString{String}
    datfile::String
    top_surface::T
    bottom_surface::B
    camber::C
end
