using ..Utils
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
Plane(surfaces::Vector{<:Aerosurface})
Plane(surfaces::Vector{<:Aerosurface}, CG::NTuple{3, Float64})
```

When `CG` is omitted, it is computed automatically at the aerodynamic center (quarter-chord of the Mean Aerodynamic Chord) of the primary surface (`surfaces[1]`):
1. Finds the spanwise location ``y_\\text{MAC} \\in [0, 1]`` where ``c(y) = \\text{MAC}`` using [`bisection`](@ref MyPackage.Utils.bisection).
2. Computes the longitudinal sweep displacement to ``y_\\text{MAC}`` via Gauss-Legendre quadrature ([`IntegrateGLQ`](@ref MyPackage.Utils.IntegrateGLQ)): ``\\Delta x_\\text{sweep} = \\frac{b}{2} \\int_0^{y_\\text{MAC}} \\tan\\Lambda(t) \\, dt``.
3. Adds the reference chord offset: ``\\Delta x_\\text{offset} = \\text{sw\\_center} \\cdot [c(0) - c(y_\\text{MAC})]``.
4. Sets the final CG at ``x_\\text{CG} = x_\\text{pos} + \\Delta x_\\text{sweep} + \\Delta x_\\text{offset} + 0.25 \\cdot \\text{MAC}``, ``y_\\text{CG} = y_\\text{pos}``, ``z_\\text{CG} = z_\\text{pos}``.
"""
struct Plane
    surfaces::Vector{<:Aerosurface}
    coeffs::Coeffs
    data::Data
end

function Plane(
    surfaces::Vector{<:Aerosurface}
)
    coeffs = Coeffs()

    s = surfaces[1]
    pos = s.pos
    y_mac = bisection(y -> s.chord(y) - s.MAC, 0.0, 1.0; tol=1e-8)
    sweep_integral = IntegrateGLQ(t -> tand(s.sweep(t)); n=5)
    sweep_length = s.b * sweep_integral(0.0, y_mac) / 2
    offset = s.sw_center * (s.chord(0.0) - s.chord(y_mac))
    CG = (pos[1] + sweep_length + offset + 0.25 * s.MAC, pos[2], pos[3])

    data = Data(; CG=CG)

    return Plane(surfaces, coeffs, data)
end

function Plane(
    surfaces::Vector{<:Aerosurface}, CG::NTuple{3, Float64}
)
    coeffs = Coeffs()

    data = Data(; CG=CG)

    return Plane(surfaces, coeffs, data)
end
