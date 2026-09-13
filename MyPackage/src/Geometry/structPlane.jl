using ..PlaneInfo

"""
    Plane

Top-level representation of an aircraft assembly consisting of multiple lifting surfaces.

# Fields

  - `surfaces::Vector{<:Aerosurface}`: Collection of constituent aerodynamic surfaces ([`Aerosurface`](@ref)).
  - `coeffs::Coeffs`: Aerodynamic coefficients container.
  - `data::Data`: Aircraft performance and geometry data, including center of gravity `CG`.

# Constructors

```julia
Plane(surfaces::Vector{<:Aerosurface}; CG=(0.0, 0.0, 0.0))
```

If `CG` is left as `(0.0, 0.0, 0.0)`, it defaults automatically to the quarter-chord of the mean aerodynamic chord
of the primary surface: `(pos[1] + 0.25 * surfaces[1].MAC, pos[2], pos[3])`.
"""
struct Plane
    surfaces::Vector{<:Aerosurface}
    coeffs::Coeffs
    data::Data
end

function Plane(
    surfaces::Vector{<:Aerosurface}; CG::NTuple{3, Float64}=(0.0, 0.0, 0.0)
)
    coeffs = Coeffs()

    if CG === (0.0, 0.0, 0.0)
        s = surfaces[1]
        pos = s.pos
        CG = (pos[1] + 0.25 * s.MAC, pos[2], pos[3])
    end

    data = Data(; CG=CG)

    return Plane(surfaces, coeffs, data)
end
