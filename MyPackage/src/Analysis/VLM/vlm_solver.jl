using LinearAlgebra
using StaticArrays
using Base.Threads

using ..Geometry

const Vec3 = SVector{3, Float64}

struct VLMSurface
    n_chord::Int
    n_span::Int
    mirror_xz::Bool
    range::UnitRange{Int}
end

struct VortexRing
    corners::NTuple{4, Vec3}
    colpt::Vec3
    normal::Vec3
    area::Float64
    surface_id::Int
end

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

function GroundTransform(rot::NTuple{2, Float64}, h::Float64, CG::Vec3)
    alpha = rot[1]
    c2a   = cosd(2 * alpha)
    s2a   = sind(2 * alpha)
    sa    = sind(alpha)
    return GroundTransform(
        h, c2a, s2a, 2 * h * s2a, 2 * h * c2a, CG[1], CG[3], sa
    )
end

struct VLMSetup
    initialized::Bool
    AIC_rings::Matrix{Float64}
    n_panels::Int
    surfaces::Vector{VLMSurface}
    wake_map::Vector{Int}
    panel_rings::Vector{VortexRing}
    ground::Bool
    h::Float64
end

# function VLMSolver(
#     plane::Plane,
#     V_inf::Float64;
#     alpha::AbstractVector{<:Float64}=[0.0],
#     beta::AbstractVector{<:Float64}=[0.0],
#     ground::Bool=false,
#     h::Float64=0.0,
#     rho::Float64=1.225,
#     epsilon2::Float64=1e-10,
#     n_chord::Int=10,
#     n_span::Int=10,
#     wake_length::Float64=3.0
#     )
#     for b in beta
#         for a in alpha
#             VLMSolver(plane, V_inf, (a, b), ground=ground, h=h, rho=rho, epsilon2=epsilon2, n_chord=n_chord, n_span=n_span, wake_length=wake_length)
#         end
#     end
# end

function VLMSolver(
    plane::Plane,
    V_inf::Float64,
    rot::NTuple{2, Float64};
    ground::Bool=false,
    h::Float64=0.0,
    rho::Float64=1.225,
    epsilon2::Float64=1e-10,
    n_chordxspan::Vector{NTuple{2, Int}}=[(0, 0)],
)
    CG = Vec3(plane.data.CG...)

    meshes = VLMMesh(plane, n_chordxspan)

    rings, n_panels, n_surfaces, surfaces, wake_map = _gen_vortex_geom(meshes)

    if ground
        AIC, AIC_rings, AIC_vx, AIC_vy, AIC_vz, RHS, V_dir = _assemble_sys(
            rings, n_panels, rot, wake_map, V_inf, h, CG, epsilon2
        )
    else
        AIC, AIC_rings, AIC_vx, AIC_vy, AIC_vz, RHS, V_dir = _assemble_sys(
            rings, n_panels, rot, wake_map, V_inf, epsilon2
        )
    end

    setup = VLMSetup(
        true, AIC_rings, n_panels, surfaces, wake_map, rings, ground, h
    )

    F = lu(AIC)
    gamma = F \ RHS

    forces = _calc_forces(
        gamma,
        n_panels,
        n_surfaces,
        surfaces,
        rings,
        rho,
        V_dir,
        V_inf,
        AIC_vx,
        AIC_vy,
        AIC_vz,
    )

    loads = _calc_loads(forces, rings, CG, n_surfaces, rot)

    # NEW: Trefftz-plane induced drag and total lift, appended as extra
    # return values.
    L_trefftz, L_dist_trefftz, D_trefftz, D_dist_trefftz = _calc_trefftz_loads(
        gamma, rings, surfaces, rho, V_dir, V_inf, rot, epsilon2
    )

    return (loads..., L_trefftz, L_dist_trefftz, D_trefftz, D_dist_trefftz)
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

            rings[(i_chord - 1) * n_span + i_span] = VortexRing(
                (A, B, C, D), colpt, n, area, surface_id
            )
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
    cr = cross(d1, d2)
    normal = norm(cr)
    n = cr / normal
    area = 0.5 * normal
    return n, area
end

function _assemble_sys(
    rings::Vector{VortexRing},
    n_panels::Int,
    rot::NTuple{2, Float64},
    wake_map::Vector{Int},
    V_inf::Float64,
    epsilon2::Float64,
)
    V_dir = _wind_dir(rot)
    inv_wake_map = _build_inv_wake_map(wake_map, n_panels)

    AIC, AIC_rings, AIC_vx, AIC_vy, AIC_vz = _assemble_AIC(
        rings, n_panels, inv_wake_map, V_dir, epsilon2
    )

    RHS = _assemble_RHS(rings, n_panels, V_dir, V_inf)

    return AIC, AIC_rings, AIC_vx, AIC_vy, AIC_vz, RHS, V_dir
end

function _assemble_sys(
    rings::Vector{VortexRing},
    n_panels::Int,
    rot::NTuple{2, Float64},
    wake_map::Vector{Int},
    V_inf::Float64,
    h::Float64,
    CG::Vec3,
    epsilon2::Float64,
)
    V_dir = _wind_dir(rot)
    inv_wake_map = _build_inv_wake_map(wake_map, n_panels)

    gt       = GroundTransform(rot, h, CG)
    ring_img = _precompute_image_corners(rings, gt)
    dir_img  = _apply_ground_transform_dir(V_dir, gt)

    AIC, AIC_rings, AIC_vx, AIC_vy, AIC_vz = _assemble_AIC(
        rings, n_panels, inv_wake_map, V_dir, ring_img, dir_img, epsilon2
    )

    RHS = _assemble_RHS(rings, n_panels, V_dir, V_inf)

    return AIC, AIC_rings, AIC_vx, AIC_vy, AIC_vz, RHS, V_dir
end

@inline function _wind_dir(rot::NTuple{2, Float64})::Vec3
    alpha, beta = rot
    ca, sa = cosd(alpha), sind(alpha)
    cb, sb = cosd(beta), sind(beta)
    return Vec3(ca * cb, -sb, sa * cb)
end

function _build_inv_wake_map(wake_map::Vector{Int}, n_panels::Int)::Vector{Int}
    inv = zeros(Int, n_panels)
    @inbounds for k in eachindex(wake_map)
        inv[wake_map[k]] = k
    end
    return inv
end

function _precompute_image_corners(
    panels::Vector{VortexRing}, gt::GroundTransform
)::Vector{NTuple{4, Vec3}}
    n   = length(panels)
    img = Vector{NTuple{4, Vec3}}(undef, n)
    @inbounds for j in 1:n
        img[j] = _apply_ground_transform(panels[j].corners, gt)
    end
    return img
end

@inline function _apply_ground_transform(
    corners::NTuple{4, Vec3}, gt::GroundTransform
)::NTuple{4, Vec3}
    ntuple(Val(4)) do i
        x, y, z = corners[i]
        x_ref, z_ref = x - gt.x_cg, gt.z_cg - z
        x_t = gt.c2a * x_ref - gt.s2a * z_ref + gt.shift_x + gt.x_cg
        z_t = gt.s2a * x_ref + gt.c2a * z_ref - gt.shift_z + gt.z_cg
        # @assert z_t < z "Mirrored point is above the original point, check h and CG values."
        @assert gt.h > gt.sa * x_ref - gt.ca * z_ref "Plane intersects the ground, increase h or adjust CG."
        Vec3(x_t, y, z_t)
    end
end

@inline function _apply_ground_transform_dir(
    dir::Vec3, gt::GroundTransform
)::Vec3
    x, y, z = dir
    x_t = gt.c2a * x - gt.s2a * z
    z_t = gt.s2a * x + gt.c2a * z
    return Vec3(x_t, y, z_t)
end

function _assemble_AIC(
    rings::Vector{VortexRing},
    n_panels::Int,
    inv_wake_map::Vector{Int},
    V_dir::Vec3,
    epsilon2::Float64,
)
    AIC = Matrix{Float64}(undef, n_panels, n_panels)
    AIC_rings = Matrix{Float64}(undef, n_panels, n_panels)
    AIC_vx = Matrix{Float64}(undef, n_panels, n_panels)
    AIC_vy = Matrix{Float64}(undef, n_panels, n_panels)
    AIC_vz = Matrix{Float64}(undef, n_panels, n_panels)

    # Parallelise over rows (collocation points)
    @threads for i in 1:n_panels
        ri = rings[i]
        colpt = ri.colpt
        ni = ri.normal

        @inbounds for j in 1:n_panels
            corners = rings[j].corners
            vel = _ring_induced_v(colpt, corners, 1.0, epsilon2)

            AIC_rings[i, j] = dot(ni, vel)

            k = inv_wake_map[j]
            if k > 0
                vel += _horseshoe_induced_v(
                    colpt, corners[4], corners[3], V_dir, 1.0, epsilon2
                )
            end

            AIC_vx[i, j] = vel[1]
            AIC_vy[i, j] = vel[2]
            AIC_vz[i, j] = vel[3]

            AIC[i, j] = dot(ni, vel)
        end
    end

    return AIC, AIC_rings, AIC_vx, AIC_vy, AIC_vz
end

function _assemble_AIC(
    rings::Vector{VortexRing},
    n_panels::Int,
    inv_wake_map::Vector{Int},
    V_dir::Vec3,
    ring_img::Vector{NTuple{4, Vec3}},
    dir_img::Vec3,
    epsilon2::Float64,
)
    AIC = Matrix{Float64}(undef, n_panels, n_panels)
    AIC_rings = Matrix{Float64}(undef, n_panels, n_panels)
    AIC_vx = Matrix{Float64}(undef, n_panels, n_panels)
    AIC_vy = Matrix{Float64}(undef, n_panels, n_panels)
    AIC_vz = Matrix{Float64}(undef, n_panels, n_panels)

    # Parallelise over rows (collocation points)
    @threads for i in 1:n_panels
        ri = rings[i]
        colpt = ri.colpt
        ni = ri.normal

        @inbounds for j in 1:n_panels
            rj = rings[j]
            corners = rj.corners
            img_corners = ring_img[j]

            vel =
                _ring_induced_v(colpt, corners, 1.0, epsilon2) -
                _ring_induced_v(colpt, img_corners, 1.0, epsilon2)

            AIC_rings[i, j] = dot(ni, vel)

            k = inv_wake_map[j]
            if k > 0
                vel +=
                    _horseshoe_induced_v(
                        colpt, corners[4], corners[3], V_dir, 1.0, epsilon2
                    ) - _horseshoe_induced_v(
                        colpt,
                        img_corners[4],
                        img_corners[3],
                        dir_img,
                        1.0,
                        epsilon2,
                    )
            end

            AIC_vx[i, j] = vel[1]
            AIC_vy[i, j] = vel[2]
            AIC_vz[i, j] = vel[3]

            AIC[i, j] = dot(ni, vel)
        end
    end

    return AIC, AIC_rings, AIC_vx, AIC_vy, AIC_vz
end

@inline function _segment_induced_v(
    P::Vec3, A::Vec3, B::Vec3, Gamma::Float64, epsilon2::Float64
)::Vec3
    r1 = A - P
    r2 = B - P
    r0 = B - A
    r0_2 = dot(r0, r0)

    # degenerate (zero-length) segment guard, independent of P.
    if r0_2 < 1e-15 # or 1e-28
        return Vec3(0.0, 0.0, 0.0)
    end

    # epsilon2 is a dimensionless core-radius fraction of the
    # local segment length^2 ([core_len2] = m^2), so the regularization
    # scales with panel size instead of being a fixed, geometry-independent
    # absolute number.
    core2 = epsilon2 * r0_2

    cr = cross(r1, r2)
    cr2 = dot(cr, cr) + core2 * core2

    norm1 = sqrt(dot(r1, r1) + core2)
    norm2 = sqrt(dot(r2, r2) + core2)

    factor =
        Gamma / (4 * pi * cr2) * (dot(r0, r1) / norm1 - dot(r0, r2) / norm2)
    return factor * cr
end

@inline function _semi_infinite_induced_v(
    P::Vec3,
    A::Vec3,
    dir::Vec3,
    Gamma::Float64,
    ref_len2::Float64,
    epsilon2::Float64,
)::Vec3
    if ref_len2 < 1e-15 # or 1e-28
        return Vec3(0.0, 0.0, 0.0)
    end

    r1 = A - P
    r1_2 = dot(r1, r1)
    core2 = epsilon2 * ref_len2

    cr = cross(r1, dir)
    cr2 = dot(cr, cr) + core2 * core2

    norm1 = sqrt(r1_2 + core2)
    cos_th = dot(r1, dir) / norm1

    factor = Gamma / (4 * pi * cr2) * (cos_th - 1.0)
    return factor * cr
end

function _ring_induced_v(
    P::Vec3,
    A::Vec3,
    B::Vec3,
    C::Vec3,
    D::Vec3,
    Gamma::Float64,
    epsilon2::Float64,
)::Vec3
    v1 = _segment_induced_v(P, A, B, Gamma, epsilon2)
    v2 = _segment_induced_v(P, B, C, Gamma, epsilon2)
    v3 = _segment_induced_v(P, C, D, Gamma, epsilon2)
    v4 = _segment_induced_v(P, D, A, Gamma, epsilon2)

    return v1 + v2 + v3 + v4
end

@inline function _ring_induced_v(
    P::Vec3, corners::NTuple{4, Vec3}, Gamma::Float64, epsilon2::Float64
)::Vec3
    return _ring_induced_v(
        P, corners[1], corners[2], corners[3], corners[4], Gamma, epsilon2
    )
end

@inline function _horseshoe_induced_v(
    P::Vec3, A::Vec3, B::Vec3, dir::Vec3, Gamma::Float64, epsilon2::Float64
)::Vec3
    ref_len2 = dot(B - A, B - A)
    return _segment_induced_v(P, A, B, Gamma, epsilon2) +
           _semi_infinite_induced_v(P, B, dir, Gamma, ref_len2, epsilon2) -
           _semi_infinite_induced_v(P, A, dir, Gamma, ref_len2, epsilon2)
end

function _assemble_RHS(
    rings::Vector{VortexRing}, n_panels::Int, V_dir::Vec3, V_inf::Float64
)::Vector{Float64}
    RHS = Vector{Float64}(undef, n_panels)

    @threads for i in 1:n_panels
        RHS[i] = -dot(rings[i].normal, V_dir * V_inf)
    end

    return RHS
end

function _calc_forces(
    gamma::Vector{Float64},
    n_panels::Int,
    n_surfaces::Int,
    surfaces::Vector{VLMSurface},
    rings::Vector{VortexRing},
    rho::Float64,
    V_dir::Vec3,
    V_inf::Float64,
    AIC_vx::Matrix{Float64},
    AIC_vy::Matrix{Float64},
    AIC_vz::Matrix{Float64},
)
    V_free = V_dir * V_inf
    forces = Vector{Matrix{Vec3}}(undef, n_surfaces)

    Vx = AIC_vx * gamma
    Vy = AIC_vy * gamma
    Vz = AIC_vz * gamma

    for i in 1:n_surfaces
        surface = surfaces[i]
        n_span = surface.n_span
        n_chord = surface.n_chord
        mirror_xz = surface.mirror_xz
        range = surface.range
        start = range.start

        forces[i] = Matrix{Vec3}(undef, (mirror_xz + 1) * n_span, n_chord)

        for j in range
            rj = rings[j]
            span_vec = rj.corners[2] - rj.corners[1]
            chord_vec = rj.corners[3] - rj.corners[2]

            row_len = (1 + mirror_xz) * n_span

            delta_gamma_c, delta_gamma_s = _calc_delta_gamma(
                gamma, j, start, row_len
            )

            V_total = V_free + Vec3(Vx[j], Vy[j], Vz[j])
            vector = delta_gamma_c * span_vec + delta_gamma_s * chord_vec

            if (j - start) % row_len == 0
                vector += gamma[j] * (rj.corners[4] - rj.corners[1])
            end

            i_span, i_chord = _j1dto2d(j, start, n_span, Val(mirror_xz))
            forces[i][i_span, i_chord] = rho * cross(V_total, vector)
        end
    end

    return forces
end

function _calc_delta_gamma(
    gamma::Vector{Float64}, j::Int, start::Int, row_len::Int
)
    delta_gamma_s =
        (((j - start + 1) % row_len == 0) ? 0.0 : gamma[j + 1]) - gamma[j]
    delta_gamma_c =
        gamma[j] - ((j - start) < row_len ? 0.0 : gamma[j - row_len])
    return delta_gamma_c, delta_gamma_s
end

function _j1dto2d(j::Int, start::Int, n_span::Int, ::Val{false})
    i_span = (j - start) % n_span + 1
    i_chord = div(j - start, n_span) + 1

    return i_span, i_chord
end

function _j1dto2d(j::Int, start::Int, n_span::Int, ::Val{true})
    i_span = (j - start) % (2 * n_span) + 1
    i_chord = div(j - start, 2 * n_span) + 1

    return i_span, i_chord
end

function _calc_loads(
    forces::Vector{Matrix{Vec3}},
    rings::Vector{VortexRing},
    CG::Vec3,
    n_surfaces::Int,
    rot::NTuple{2, Float64},
)
    alpha = rot[1]
    beta = rot[2]
    sa = sind(alpha)
    ca = cosd(alpha)
    sb = sind(beta)
    cb = cosd(beta)

    FX = fill(0.0, n_surfaces) # Force in x direction
    FY = fill(0.0, n_surfaces) # Force in y direction
    FZ = fill(0.0, n_surfaces) # Force in z direction

    FX_dist = Vector{Vector{Float64}}(undef, n_surfaces)
    FY_dist = Vector{Vector{Float64}}(undef, n_surfaces)
    FZ_dist = Vector{Vector{Float64}}(undef, n_surfaces)

    L = fill(0.0, n_surfaces) # Lift
    D = fill(0.0, n_surfaces) # Drag

    L_dist = Vector{Vector{Float64}}(undef, n_surfaces)
    D_dist = Vector{Vector{Float64}}(undef, n_surfaces)

    M = fill(0.0, n_surfaces) # Pitching moment
    Ml = fill(0.0, n_surfaces) # Roll moment
    N = fill(0.0, n_surfaces) # Yaw moment

    M_dist = Vector{Vector{Float64}}(undef, n_surfaces)
    Ml_dist = Vector{Vector{Float64}}(undef, n_surfaces)
    N_dist = Vector{Vector{Float64}}(undef, n_surfaces)

    i1d = 1

    for i in 1:n_surfaces
        surface_forces = forces[i]
        n_span, n_chord = size(surface_forces)

        FX_dist[i] = fill(0.0, n_span)
        FY_dist[i] = fill(0.0, n_span)
        FZ_dist[i] = fill(0.0, n_span)

        L_dist[i] = fill(0.0, n_span)
        D_dist[i] = fill(0.0, n_span)

        M_dist[i] = fill(0.0, n_span)
        Ml_dist[i] = fill(0.0, n_span)
        N_dist[i] = fill(0.0, n_span)

        for j in 1:n_span
            for k in 1:n_chord
                fx = surface_forces[j, k][1]
                fy = surface_forces[j, k][2]
                fz = surface_forces[j, k][3]

                FX[i] += fx
                FY[i] += fy
                FZ[i] += fz

                FX_dist[i][j] += fx
                FY_dist[i][j] += fy
                FZ_dist[i][j] += fz

                lift = -fx * sa + fz * ca
                drag = fx * ca * cb - fy * sb + fz * sa * cb

                L[i] += lift
                D[i] += drag

                L_dist[i][j] += lift
                D_dist[i][j] += drag

                rj = rings[i1d]
                i1d += 1
                moment_arm = rj.colpt - CG

                (ml, m, n) = cross(moment_arm, surface_forces[j, k])

                M[i] += m
                Ml[i] += ml
                N[i] += n

                M_dist[i][j] += m
                Ml_dist[i][j] += ml
                N_dist[i][j] += n
            end
        end
    end

    return (
        FX,
        FX_dist,
        FY,
        FY_dist,
        FZ,
        FZ_dist,
        L,
        L_dist,
        D,
        D_dist,
        M,
        M_dist,
        Ml,
        Ml_dist,
        N,
        N_dist,
    )
end

@inline function _trefftz_filament_induced_v(
    P::Vec3,
    A::Vec3,
    Gamma::Float64,
    V_dir::Vec3,
    ref_len2::Float64,
    epsilon2::Float64,
)::Vec3
    r1 = A - P
    r1 -= dot(r1, V_dir) * V_dir
    cr = cross(r1, V_dir)
    h2 = dot(cr, cr) + epsilon2 * ref_len2
    return (-Gamma / (4 * pi * h2)) * cr
end

function _assemble_trefftz_wake(
    gamma::Vector{Float64},
    rings::Vector{VortexRing},
    surfaces::Vector{VLMSurface},
)::Tuple{Vector{Vec3}, Vector{Float64}, Vector{Float64}}
    n_fil = 0
    @inbounds for surf in surfaces
        n_fil += (1 + surf.mirror_xz) * surf.n_span + 1
    end

    fil_pos = Vector{Vec3}(undef, n_fil)
    fil_strength = Vector{Float64}(undef, n_fil)
    fil_ref2 = Vector{Float64}(undef, n_fil)

    m = 0
    for surf in surfaces
        row_len = (1 + surf.mirror_xz) * surf.n_span
        te_start = surf.range.stop - row_len + 1

        prev_g = 0.0
        @inbounds for k in 1:row_len
            j = te_start + k - 1
            g = gamma[j]
            rj = rings[j]
            edge2 = dot(
                rj.corners[3] - rj.corners[4], rj.corners[3] - rj.corners[4]
            )

            m += 1
            fil_pos[m] = rj.corners[4]
            fil_strength[m] = g - prev_g
            fil_ref2[m] = edge2
            prev_g = g

            if k == row_len
                m += 1
                fil_pos[m] = rj.corners[3]
                fil_strength[m] = -g
                fil_ref2[m] = edge2
            end
        end
    end

    return fil_pos, fil_strength, fil_ref2
end

function _assemble_trefftz_te(
    surfaces::Vector{VLMSurface}
)::Tuple{Vector{Int}, Vector{Int}, Vector{Int}, Vector{Int}}
    n_surfaces = length(surfaces)
    row_lens = Vector{Int}(undef, n_surfaces)
    n_te = 0
    @inbounds for (i, surf) in enumerate(surfaces)
        row_len = (1 + surf.mirror_xz) * surf.n_span
        row_lens[i] = row_len
        n_te += row_len
    end

    te_j = Vector{Int}(undef, n_te)
    te_surf = Vector{Int}(undef, n_te)
    te_k = Vector{Int}(undef, n_te)

    t = 0
    for (i, surf) in enumerate(surfaces)
        row_len = row_lens[i]
        te_start = surf.range.stop - row_len + 1
        @inbounds for k in 1:row_len
            t += 1
            te_j[t] = te_start + k - 1
            te_surf[t] = i
            te_k[t] = k
        end
    end

    return te_j, te_surf, te_k, row_lens
end

function _calc_trefftz_loads(
    gamma::Vector{Float64},
    rings::Vector{VortexRing},
    surfaces::Vector{VLMSurface},
    rho::Float64,
    V_dir::Vec3,
    V_inf::Float64,
    rot::NTuple{2, Float64},
    epsilon2::Float64,
)
    n_surfaces = length(surfaces)

    fil_pos, fil_strength, fil_ref2 = _assemble_trefftz_wake(
        gamma, rings, surfaces
    )
    n_fil = length(fil_pos)

    te_j, te_surf, te_k, row_lens = _assemble_trefftz_te(surfaces)
    n_te = length(te_j)

    V_free = V_dir * V_inf
    sa, ca = sind(rot[1]), cosd(rot[1])

    d_te = Vector{Float64}(undef, n_te)
    l_te = Vector{Float64}(undef, n_te)

    # Parallelise over trailing-edge panels (mirrors _assemble_AIC's
    # per-collocation-point threading; the O(n_te * n_fil) self-induction
    # sum below is the expensive part).
    @threads for t in 1:n_te
        j = te_j[t]
        rj = rings[j]
        span_vec = rj.corners[3] - rj.corners[4]
        midpt = 0.5 * (rj.corners[3] + rj.corners[4])
        span_vec_tp = span_vec - dot(span_vec, V_dir) * V_dir

        w = Vec3(0.0, 0.0, 0.0)
        @inbounds for m in 1:n_fil
            w += _trefftz_filament_induced_v(
                midpt, fil_pos[m], fil_strength[m], V_dir, fil_ref2[m], epsilon2
            )
        end

        f = rho * cross(V_free, gamma[j] * span_vec_tp)

        d_te[t] = dot(rho * cross(w, gamma[j] * span_vec_tp), V_dir)
        l_te[t] = -f[1] * sa + f[3] * ca
    end

    D = fill(0.0, n_surfaces)
    L = fill(0.0, n_surfaces)
    D_dist = Vector{Vector{Float64}}(undef, n_surfaces)
    L_dist = Vector{Vector{Float64}}(undef, n_surfaces)

    for i in 1:n_surfaces
        D_dist[i] = fill(0.0, row_lens[i])
        L_dist[i] = fill(0.0, row_lens[i])
    end

    @inbounds for t in 1:n_te
        i = te_surf[t]
        k = te_k[t]
        D_dist[i][k] = d_te[t]
        L_dist[i][k] = l_te[t]
        D[i] += d_te[t]
        L[i] += l_te[t]
    end

    return L, L_dist, D, D_dist
end
