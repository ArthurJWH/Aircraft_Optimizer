using ..Utils

"""
    Aerosurface

    A mutable struct that represents an aerodynamic surface, such as wing, horizontal stabilizer, or vertical stabilizer.

    Fields
    ------
    name : String, optional
        The name of the aerosurface. Default is "Aerosurface".
    mirror_xz : Bool, optional
        Indicates whether the aerosurface is mirrored across the xz-plane. Default is true.
    vertical : Bool, optional
        Indicates whether the aerosurface is vertical. Default is false.
    pos : Tuple{Float64, Float64, Float64}, optional
        The position of the aerosurface root leading edge in 3D space. Default is (0.0, 0.0, 0.0).
    rot : Tuple{Float64, Float64, Float64}, optional
        The rotation of the aerosurface in 3D space. Default is (0.0, 0.0, 0.0).
        Currently not used in the code, but can be used for future development.
    b : Float64, optional
        The span of the aerosurface. Default is 1.0.
    ys : Vector{Float64}, optional
        The spanwise fraction coordinates of the airfoil positions in the aerosurface. Default is [0.0, 1.0].
    airfoils : Vector{<:Airfoil}, optional
        The airfoils used along the span of the aerosurface. Default is a vector of two plain airfoils.
    chord : Function, optional
        A function that defines the chord length distribution in terms of span fractions of the aerosurface. Default is a constant function returning 1.0.
    twist : Function, optional
        A function that defines the twist distribution in terms of span fractions of the aerosurface. Default is a constant function returning 0.0.
    tw_center : Float64, optional
        The chordwise fraction of each section where the twist is applied. Default is 0.25.
    sweep : Function, optional
        A function that defines the sweep distribution in terms of span fractions of the aerosurface. Default is a constant function returning 0.0.
    sw_center : Float64, optional
        The chordwise fraction of each section where the sweep is applied. Default is 0.25.
    dihedral : Function, optional
        A function that defines the dihedral distribution in terms of span fractions of the aerosurface. Default is a constant function returning 0.0.
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
