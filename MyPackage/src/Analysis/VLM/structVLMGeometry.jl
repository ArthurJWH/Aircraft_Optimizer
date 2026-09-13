# VLMGeometry: the mesh and vortex-ring geometry for one discretization of
# a Plane, plus the mesh-generation functions that build it. Depends on
# vlm_types.jl (VortexRing, VLMSurface) and vlm_polar.jl only through
# doc cross-references, not code.

"""
    VLMGeometry

The geometric and topological discretization of a [`Plane`](@ref MyPackage.Geometry.Plane) for Vortex Lattice Method analysis.

Contains the underlying surface meshes, discrete vortex rings, panel normals, collocation points, and reference aerodynamic quantities.
Decoupled from flow conditions (`V_inf`, `rho`, `ground`), allowing a single `VLMGeometry` instance to be shared across
multiple [`VLMSetup`](@ref) configurations without re-meshing.

# Fields

  - `meshes::Vector{VLMMesh}`: Per-surface panel meshes containing vertex coordinates.
  - `rings::Vector{VortexRing}`: Global flattened vector of all discrete quadrilateral vortex rings.
  - `n_panels::Int`: Total number of vortex panels across all surfaces.
  - `n_surfaces::Int`: Number of lifting surfaces in the plane.
  - `surfaces::Vector{VLMSurface}`: Metadata and index ranges for each surface.
  - `wake_map::Vector{Int}`: Global indices of trailing-edge panels shedding wake vortex filaments.
  - `CG::Vec3`: Center of gravity coordinates `(x, y, z)` in `m` about which moments are evaluated.
  - `Sref::Float64`: Reference planform area in m² (inherited from primary wing `plane.surfaces[1].S`).
  - `cref::Float64`: Reference chord (MAC) in `m` (inherited from primary wing `plane.surfaces[1].MAC`).
  - `bref::Float64`: Reference wingspan in `m` (inherited from primary wing `plane.surfaces[1].b`).

# Constructors

```julia
VLMGeometry(plane::Plane, n_chordxspan::Vector{NTuple{2, Int}})
```

`n_chordxspan` specifies `(n_chord, n_span)` panel counts for each surface of the plane.
"""
struct VLMGeometry
    meshes::Vector{VLMMesh}
    rings::Vector{VortexRing}
    n_panels::Int
    n_surfaces::Int
    surfaces::Vector{VLMSurface}
    wake_map::Vector{Int}
    CG::Vec3
    Sref::Float64
    cref::Float64
    bref::Float64
end

function VLMGeometry(
    plane::Plane, n_chordxspan::Vector{NTuple{2, Int}}=[(0, 0)]
)
    meshes = VLMMesh(plane, n_chordxspan)
    rings, n_panels, n_surfaces, surfaces, wake_map = _gen_vortex_geom(meshes)
    CG = Vec3(plane.data.CG...)
    Sref = plane.surfaces[1].S
    cref = plane.surfaces[1].MAC
    bref = plane.surfaces[1].b
    return VLMGeometry(
        meshes,
        rings,
        n_panels,
        n_surfaces,
        surfaces,
        wake_map,
        CG,
        Sref,
        cref,
        bref,
    )
end

function _gen_vortex_geom(meshes::AbstractVector{<:VLMMesh})
    n_meshes = length(meshes)
    surfaces = Vector{VLMSurface}(undef, n_meshes)
    rings    = VortexRing[]
    wake_map = Int[]

    sizehint!(rings, 1024) #TODO: improve sizehint
    sizehint!(wake_map, 128)

    i_start = 0
    for (surface_id, mesh) in enumerate(meshes)
        vertices = mesh.vertices
        mirror_xz = mesh.mirror_xz
        sz = size(vertices)
        # @assert sz[1] == 3 "First dimension must be 3 (x,y,z)"
        n_span = sz[2] - 1
        n_chord = sz[3] - 1

        _gen_mesh_rings!(
            rings, vertices, n_span, n_chord, surface_id, Val(mirror_xz)
        )

        i_end = i_start + (1 + mirror_xz) * n_span * n_chord

        surfaces[surface_id] = VLMSurface(
            n_chord, n_span, mirror_xz, (i_start + 1):i_end
        )
        append!(wake_map, (i_end - (1 + mirror_xz) * n_span + 1):i_end)

        i_start = i_end
    end

    n_panels = i_start
    return rings, n_panels, length(meshes), surfaces, wake_map
end

function _gen_mesh_rings!(
    rings::Vector{VortexRing},
    vertices::Array{Float64, 3},
    n_span::Int,
    n_chord::Int,
    surface_id::Int,
    ::Val{false},
)::Nothing
    corners = Matrix{Vec3}(undef, n_span + 1, n_chord + 1)
    col_ref = Matrix{Vec3}(undef, n_span + 1, n_chord + 1)

    for i_chord in 1:n_chord
        @inbounds for i_span in 1:(n_span + 1)
            le = Vec3(
                vertices[1, i_span, i_chord],
                vertices[2, i_span, i_chord],
                vertices[3, i_span, i_chord],
            )
            te = Vec3(
                vertices[1, i_span, i_chord + 1],
                vertices[2, i_span, i_chord + 1],
                vertices[3, i_span, i_chord + 1],
            )
            corners[i_span, i_chord] = le + 0.25 * (te - le)
            col_ref[i_span, i_chord] = le + 0.75 * (te - le)
        end
    end
    @inbounds for i_span in 1:(n_span + 1)
        corners[i_span, n_chord + 1] = Vec3(
            vertices[1, i_span, n_chord + 1],
            vertices[2, i_span, n_chord + 1],
            vertices[3, i_span, n_chord + 1],
        )
    end

    @inbounds for i_chord in 1:n_chord
        for i_span in 1:n_span
            A = corners[i_span, i_chord]
            B = corners[i_span + 1, i_chord]
            C = corners[i_span + 1, i_chord + 1]
            D = corners[i_span, i_chord + 1]

            colpt =
                0.5 * (col_ref[i_span, i_chord] + col_ref[i_span + 1, i_chord])
            n, area = _panel_normal_area(A, B, C, D)

            push!(rings, VortexRing((A, B, C, D), colpt, n, area, surface_id))
        end
    end

    return nothing
end

function _gen_mesh_rings!(
    rings::Vector{VortexRing},
    vertices::Array{Float64, 3},
    n_span::Int,
    n_chord::Int,
    surface_id::Int,
    ::Val{true},
)::Nothing
    corners = Matrix{Vec3}(undef, n_span + 1, n_chord + 1)
    col_ref = Matrix{Vec3}(undef, n_span + 1, n_chord + 1)
    rings_helper = Vector{VortexRing}(undef, 2 * n_span * n_chord)

    for i_chord in 1:n_chord
        @inbounds for i_span in 1:(n_span + 1)
            le = Vec3(
                vertices[1, i_span, i_chord],
                vertices[2, i_span, i_chord],
                vertices[3, i_span, i_chord],
            )
            te = Vec3(
                vertices[1, i_span, i_chord + 1],
                vertices[2, i_span, i_chord + 1],
                vertices[3, i_span, i_chord + 1],
            )
            corners[i_span, i_chord] = le + 0.25 * (te - le)
            col_ref[i_span, i_chord] = le + 0.75 * (te - le)
        end
    end
    @inbounds for i_span in 1:(n_span + 1)
        corners[i_span, n_chord + 1] = Vec3(
            vertices[1, i_span, n_chord + 1],
            vertices[2, i_span, n_chord + 1],
            vertices[3, i_span, n_chord + 1],
        )
    end

    @inbounds for i_chord in 1:n_chord
        for i_span in 1:n_span
            A = corners[i_span, i_chord]
            B = corners[i_span + 1, i_chord]
            C = corners[i_span + 1, i_chord + 1]
            D = corners[i_span, i_chord + 1]
            colpt =
                0.5 * (col_ref[i_span, i_chord] + col_ref[i_span + 1, i_chord])
            n, area = _panel_normal_area(A, B, C, D)

            rings_helper[(i_chord - 1) * 2 * n_span + n_span + i_span] = VortexRing(
                (A, B, C, D), colpt, n, area, surface_id
            )

            A_m = (A[1], -A[2], A[3])
            B_m = (B[1], -B[2], B[3])
            C_m = (C[1], -C[2], C[3])
            D_m = (D[1], -D[2], D[3])
            colpt_m = Vec3(colpt[1], -colpt[2], colpt[3])
            n_m = Vec3(n[1], -n[2], n[3])

            rings_helper[(i_chord - 1) * 2 * n_span + n_span - i_span + 1] = VortexRing(
                (B_m, A_m, D_m, C_m), colpt_m, n_m, area, surface_id
            )
        end
    end

    append!(rings, rings_helper)

    return nothing
end

@inline function _panel_normal_area(
    A::Vec3, B::Vec3, C::Vec3, D::Vec3
)::Tuple{Vec3, Float64}
    d1 = C - A   # diagonal 1
    d2 = D - B   # diagonal 2
    cr = cross(d2, d1)
    normal = norm(cr)
    n = cr / normal
    area = 0.5 * normal
    return n, area
end
