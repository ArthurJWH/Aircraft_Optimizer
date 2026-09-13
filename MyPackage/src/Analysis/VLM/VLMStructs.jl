# Core value types used throughout the VLM solver: the Vec3 alias, panel/
# ring/ground-transform geometry primitives.

"""
    Vec3

Convenience alias for a 3-element static vector `StaticArrays.SVector{3, Float64}`, used for all 3D points, directions, and forces.
"""
const Vec3 = SVector{3, Float64}

"""
    VLMSurface

Represents structural metadata for one aerodynamic lifting surface in the VLM discretization.

# Fields
- `n_chord::Int`: Number of chordwise panels.
- `n_span::Int`: Number of spanwise panels per semi-span.
- `mirror_xz::Bool`: Whether the surface is mirrored across the XZ-plane.
- `range::UnitRange{Int}`: Index range of vortex rings in the global rings vector corresponding to this surface.
"""
struct VLMSurface
    n_chord::Int
    n_span::Int
    mirror_xz::Bool
    range::UnitRange{Int}
end

"""
    VortexRing

A discrete quadrilateral vortex ring panel forming the building block of the Vortex Lattice Method.

# Fields
- `corners::NTuple{4, Vec3}`: Coordinates of the four panel corners in counter-clockwise order ``(A, B, C, D)``:
  - ``A``: Front-left (leading quarter-chord, left)
  - ``B``: Front-right (leading quarter-chord, right)
  - ``C``: Rear-right (trailing quarter-chord, right)
  - ``D``: Rear-left (trailing quarter-chord, left)
- `colpt::Vec3`: 3D coordinates of the panel collocation point, located at the panel's 3/4-chord line and centered spanwise.
- `normal::Vec3`: Unit outward normal vector ``\\vec{n}`` defining the local surface orientation.
- `area::Float64`: Area of the quadrilateral panel, in `m²`.
- `surface_id::Int`: Identifier of the parent [`VLMSurface`](@ref) containing this ring.
"""
struct VortexRing
    corners::NTuple{4, Vec3}
    colpt::Vec3
    normal::Vec3
    area::Float64
    surface_id::Int
end

"""
    GroundTransform

Encapsulates coordinate transformation parameters used to evaluate ground-effect image vortex systems.

When `ground == true`, the flat ground plane at ``z = -h`` is modeled by mirroring all vortex rings and the wake
across the ground plane to enforce zero normal velocity (non-penetration boundary condition) at the ground.

# Fields
- `h::Float64`: Height of the aircraft center of gravity `CG` above the ground plane, in `m`.
- `c2a::Float64`: Precomputed ``\\cos(2\\alpha)``.
- `s2a::Float64`: Precomputed ``\\sin(2\\alpha)``.
- `shift_x::Float64`: Longitudinal ground image translation offset, in `m`.
- `shift_z::Float64`: Vertical ground image translation offset, in `m`.
- `x_cg::Float64`: Center of gravity X coordinate, in `m`.
- `z_cg::Float64`: Center of gravity Z coordinate, in `m`.
- `sa::Float64`: Precomputed ``\\sin(\\alpha)``.
"""
struct GroundTransform
    h::Float64
    c2a::Float64
    s2a::Float64
    shift_x::Float64
    shift_z::Float64
    x_cg::Float64
    z_cg::Float64
    sa::Float64
end

"""
    GroundTransform(rot::NTuple{2, Float64}, h::Float64, CG::Vec3)

Constructs a [`GroundTransform`](@ref) for flow rotation `rot = (alpha, beta)` and ground height `h`.
"""
function GroundTransform(rot::NTuple{2, Float64}, h::Float64, CG::Vec3)
    alpha = rot[1]
    c2a   = cosd(2 * alpha)
    s2a   = sind(2 * alpha)
    sa    = sind(alpha)
    return GroundTransform(
        h, c2a, s2a, 2 * h * s2a, 2 * h * c2a, CG[1], CG[3], sa
    )
end
