using ..Utils

"""
    Aerosurface

Parametric representation of an aerodynamic lifting surface (main wing, horizontal stabilizer, vertical tail, canard, or winglet).

Geometric distributions (`chord`, `twist`, `sweep`, `dihedral`) are defined as functions of the non-dimensional semi-span station
``y \\in [0, 1]`` (where ``y = 0`` corresponds to the root and ``y = 1`` corresponds to the tip).

# Fields
- `name::String`: User-defined surface name (default: `"Aerosurface"`).
- `mirror_xz::Bool`: Whether the surface is mirrored across the XZ-plane (Y -> -Y) to produce a symmetric port semi-span (default: `true`).
- `vertical::Bool`: Whether the surface is oriented vertically (span along Z axis) (default: `false`).
- `pos::Tuple{Float64, Float64, Float64}`: Global 3D offset `(x, y, z)` in `m` of the root leading edge (default: `(0.0, 0.0, 0.0)`).
- `rot::Tuple{Float64, Float64, Float64}`: 3D rotation angles (reserved for future kinematic transformations).
- `b::Float64`: Total physical wingspan in `m` (inclusive of mirrored semi-span if `mirror_xz = true`).
- `S::Float64`: Reference planform area, in `m²`, computed via numerical quadrature ``S = b \\int_0^1 c(y) dy``.
- `AR::Float64`: Aspect ratio ``AR = b^2 / S``.
- `MGC::Float64`: Mean geometric chord in `m`, ``MGC = \\int_0^1 c(y) dy``.
- `MAC::Float64`: Mean aerodynamic chord in `m`, ``MAC = \\frac{1}{MGC} \\int_0^1 c(y)^2 dy``.
- `ys::Vector{Float64}`: Spanwise stations ``y_i \\in [0, 1]`` corresponding to the defined `airfoils`.
- `airfoils::Vector{<:Airfoil}`: Collection of [`Airfoil`](@ref) profiles placed at each station in `ys`.
- `chord::Function`: Chord distribution function ``c(y)`` returning physical chord length in `m`.
- `twist::Function`: Geometric twist distribution function ``\\theta(y)`` in `deg` (positive pitches section up).
- `tw_center::Float64`: Chordwise fraction of section about which twist rotation is applied (default: `0.25`).
- `sweep::Function`: Sweep angle function ``\\Lambda(y)`` in `deg` (positive sweeps aft +X).
- `sw_center::Float64`: Chordwise fraction reference line for sweep offset (default: `0.25`).
- `dihedral::Function`: Dihedral angle function ``\\Gamma(y)`` in `deg` (positive deflects tip upward +Z).

# Example
```julia
plain = airfoil_from_dat("assets/airfoils/Plain/Plain.dat")
wing = Aerosurface(
    name="MainWing",
    airfoils=[plain, plain],
    b=6.0,
    chord=y -> 2 * sqrt(1 - y^2),  # Elliptical chord distribution
    sw_center=0.5,
)
```
"""
mutable struct Aerosurface{chordF, twistF, sweepF, dihedralF}
    # y is the spanwise fraction coordinate
    # x is the chordwise fraction coordinate

    name::String
    mirror_xz::Bool
    vertical::Bool
    pos::Tuple{Float64, Float64, Float64}
    rot::Tuple{Float64, Float64, Float64}

    b::Float64
    S::Float64
    AR::Float64
    MGC::Float64
    MAC::Float64
    # Swet::Float64 future drag build up

    ys::Vector{Float64}
    airfoils::Vector{<:Airfoil}

    # ys_mesh::Vector{Float64}

    chord::chordF
    twist::twistF
    tw_center::Float64
    sweep::sweepF
    sw_center::Float64
    dihedral::dihedralF
end

function Aerosurface(;
    name::String                          = "Aerosurface",
    mirror_xz::Bool                       = true,
    vertical::Bool                        = false,
    pos::Tuple{Float64, Float64, Float64} = (0.0, 0.0, 0.0),
    rot::Tuple{Float64, Float64, Float64} = (0.0, 0.0, 0.0),
    b::Float64                            = 1.0,
    ys::Vector{Float64}                   = [0.0, 1.0],
    airfoils::Vector{<:Airfoil}           = [Airfoil("../../assets/airfoils/Plain/Plain.dat"), Airfoil("../../assets/airfoils/Plain/Plain.dat")],
    chord::chordF                         = y -> 1.0,
    twist::twistF                         = y -> 0.0,
    tw_center::Float64                    = 0.25,
    sweep::sweepF                         = y -> 0.0,
    sw_center::Float64                    = 0.25,
    dihedral::dihedralF                   = y -> 0.0,
) where {chordF, twistF, sweepF, dihedralF}
    MGC = IntegrateGLQ(chord)(0.0, 1.0)
    S = MGC * b
    AR = b^2 / S
    MAC = IntegrateGLQ(x -> chord(x)^2)(0.0, 1.0) / MGC # Assume wing is symmetric about the centerline
    return Aerosurface{chordF, twistF, sweepF, dihedralF}(
        name,
        mirror_xz,
        vertical,
        pos,
        rot,
        b,
        S,
        AR,
        MGC,
        MAC,
        ys,
        airfoils,
        chord,
        twist,
        tw_center,
        sweep,
        sw_center,
        dihedral,
    )
end
