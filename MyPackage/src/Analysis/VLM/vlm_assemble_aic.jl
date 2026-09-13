# The angle-independent + angle-dependent AIC assembly: wind direction,
# ground-image transforms, the Biot-Savart ring/horseshoe induced-velocity
# evaluations, and RHS assembly. This is the numerical core -- the part
# that actually implements the vortex lattice method. Pure functions only,
# no structs defined here, so it has no ordering dependency on the files
# around it beyond needing vlm_types.jl's types to exist by call time.

function _assemble_sys(
    rings::Vector{VortexRing},
    n_panels::Int,
    rot::NTuple{2, Float64},
    wake_map::Vector{Int},
    V_inf::Float64,
    ring_geom::NTuple{3, Matrix{Float64}},
    epsilon2::Float64,
)
    V_dir = _wind_dir(rot)
    inv_wake_map = _build_inv_wake_map(wake_map, n_panels)

    AIC, AIC_vx, AIC_vy, AIC_vz = _assemble_AIC(
        rings, n_panels, inv_wake_map, V_dir, ring_geom..., epsilon2
    )

    RHS = _assemble_RHS(rings, n_panels, V_dir, V_inf)

    return AIC, AIC_vx, AIC_vy, AIC_vz, RHS, V_dir
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

    return AIC, AIC_vx, AIC_vy, AIC_vz, RHS, V_dir
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

# The angle-independent half of AIC assembly: the pure ring-induced
# velocity at every (collocation point, ring) pair, with no wake/horseshoe
# term (that part depends on wind direction, computed separately per
# condition in _assemble_AIC below). Depends only on ring geometry and
# epsilon2 -- safe to compute once per VLMSetup and reuse for every
# (alpha, beta) solved under it (see _ring_geom!). NOT valid to reuse
# across different `ground` configurations (irrelevant here, this path is
# ground=false only) or different epsilon2 (already enforced by living on
# one VLMSetup, which fixes epsilon2).
function _assemble_ring_geom(
    rings::Vector{VortexRing}, n_panels::Int, epsilon2::Float64
)::NTuple{3, Matrix{Float64}}
    Vx_rings = Matrix{Float64}(undef, n_panels, n_panels)
    Vy_rings = Matrix{Float64}(undef, n_panels, n_panels)
    Vz_rings = Matrix{Float64}(undef, n_panels, n_panels)

    colpts = Vector{Vec3}(undef, n_panels)
    @inbounds for i in 1:n_panels
        colpts[i] = rings[i].colpt
    end

    @threads for j in 1:n_panels
        corners = rings[j].corners
        ring_inv = _ring_invariants(corners, epsilon2)
        @inbounds for i in 1:n_panels
            vel = _ring_induced_v(colpts[i], corners, ring_inv, 1.0)
            Vx_rings[i, j] = vel[1]
            Vy_rings[i, j] = vel[2]
            Vz_rings[i, j] = vel[3]
        end
    end

    return (Vx_rings, Vy_rings, Vz_rings)
end

# Angle-dependent AIC assembly given a precomputed, angle-independent
# ring_geom (see _assemble_ring_geom): adds only the wake/horseshoe term,
# which is the part that actually depends on V_dir, instead of
# recomputing the O(n_panels^2) ring-induced Biot-Savart evaluation on
# every call. Mathematically identical to (in fact, the same floating-point
# operations in the same order as) computing everything from scratch: since
# dot() is linear, dot(ni, vel_rings + vel_hs) == dot(ni, vel_rings) +
# dot(ni, vel_hs), so nothing is approximated here.
function _assemble_AIC(
    rings::Vector{VortexRing},
    n_panels::Int,
    inv_wake_map::Vector{Int},
    V_dir::Vec3,
    Vx_rings::Matrix{Float64},
    Vy_rings::Matrix{Float64},
    Vz_rings::Matrix{Float64},
    epsilon2::Float64,
)
    AIC = Matrix{Float64}(undef, n_panels, n_panels)
    AIC_vx = Matrix{Float64}(undef, n_panels, n_panels)
    AIC_vy = Matrix{Float64}(undef, n_panels, n_panels)
    AIC_vz = Matrix{Float64}(undef, n_panels, n_panels)

    colpts = Vector{Vec3}(undef, n_panels)
    normals = Vector{Vec3}(undef, n_panels)
    @inbounds for i in 1:n_panels
        colpts[i] = rings[i].colpt
        normals[i] = rings[i].normal
    end

    @threads for j in 1:n_panels
        corners = rings[j].corners
        k = inv_wake_map[j]

        if k > 0
            hs_inv = _horseshoe_invariants(corners[4], corners[3], epsilon2)
            @inbounds for i in 1:n_panels
                ni = normals[i]
                vel_hs = _horseshoe_induced_v(
                    colpts[i], corners[4], corners[3], V_dir, hs_inv, 1.0
                )
                vx = Vx_rings[i, j] + vel_hs[1]
                vy = Vy_rings[i, j] + vel_hs[2]
                vz = Vz_rings[i, j] + vel_hs[3]
                AIC_vx[i, j] = vx
                AIC_vy[i, j] = vy
                AIC_vz[i, j] = vz
                AIC[i, j] = dot(ni, Vec3(vx, vy, vz))
            end
        else
            @inbounds for i in 1:n_panels
                ni = normals[i]
                vx, vy, vz = Vx_rings[i, j], Vy_rings[i, j], Vz_rings[i, j]
                AIC_vx[i, j] = vx
                AIC_vy[i, j] = vy
                AIC_vz[i, j] = vz
                AIC[i, j] = dot(ni, Vec3(vx, vy, vz))
            end
        end
    end

    return AIC, AIC_vx, AIC_vy, AIC_vz
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

    colpts = Vector{Vec3}(undef, n_panels)
    normals = Vector{Vec3}(undef, n_panels)
    @inbounds for i in 1:n_panels
        colpts[i] = rings[i].colpt
        normals[i] = rings[i].normal
    end

    @threads for j in 1:n_panels
        rj = rings[j]
        corners = rj.corners
        img_corners = ring_img[j]
        ring_inv = _ring_invariants(corners, epsilon2)
        img_ring_inv = _ring_invariants(img_corners, epsilon2)
        k = inv_wake_map[j]

        if k > 0
            hs_inv = _horseshoe_invariants(corners[4], corners[3], epsilon2)
            img_hs_inv = _horseshoe_invariants(
                img_corners[4], img_corners[3], epsilon2
            )
            @inbounds for i in 1:n_panels
                P = colpts[i]
                ni = normals[i]

                vel =
                    _ring_induced_v(P, corners, ring_inv, 1.0) -
                    _ring_induced_v(P, img_corners, img_ring_inv, 1.0)
                AIC_rings[i, j] = dot(ni, vel)

                vel +=
                    _horseshoe_induced_v(
                        P, corners[4], corners[3], V_dir, hs_inv, 1.0
                    ) - _horseshoe_induced_v(
                        P,
                        img_corners[4],
                        img_corners[3],
                        dir_img,
                        img_hs_inv,
                        1.0,
                    )

                AIC_vx[i, j] = vel[1]
                AIC_vy[i, j] = vel[2]
                AIC_vz[i, j] = vel[3]
                AIC[i, j] = dot(ni, vel)
            end
        else
            @inbounds for i in 1:n_panels
                P = colpts[i]
                ni = normals[i]

                vel =
                    _ring_induced_v(P, corners, ring_inv, 1.0) -
                    _ring_induced_v(P, img_corners, img_ring_inv, 1.0)
                AIC_rings[i, j] = dot(ni, vel)

                AIC_vx[i, j] = vel[1]
                AIC_vy[i, j] = vel[2]
                AIC_vz[i, j] = vel[3]
                AIC[i, j] = dot(ni, vel)
            end
        end
    end

    return AIC, AIC_rings, AIC_vx, AIC_vy, AIC_vz
end

@inline function _segment_invariants(
    A::Vec3, B::Vec3, epsilon2::Float64
)::Tuple{Vec3, Float64, Float64}
    r0 = B - A
    r0_2 = dot(r0, r0)
    # epsilon2 is a dimensionless core-radius fraction of the local segment
    # length^2 ([core2] = m^2), so the regularization scales with panel
    # size instead of being a fixed, geometry-independent absolute number.
    core2 = epsilon2 * r0_2
    return r0, r0_2, core2
end

@inline function _ring_invariants(corners::NTuple{4, Vec3}, epsilon2::Float64)
    A, B, C, D = corners
    e1 = _segment_invariants(A, B, epsilon2)
    e2 = _segment_invariants(B, C, epsilon2)
    e3 = _segment_invariants(C, D, epsilon2)
    e4 = _segment_invariants(D, A, epsilon2)
    return (
        r0=(e1[1], e2[1], e3[1], e4[1]),
        r0_2=(e1[2], e2[2], e3[2], e4[2]),
        core2=(e1[3], e2[3], e3[3], e4[3]),
    )
end

@inline function _horseshoe_invariants(A::Vec3, B::Vec3, epsilon2::Float64)
    r0, r0_2, core2 = _segment_invariants(A, B, epsilon2)
    return (r0=r0, r0_2=r0_2, core2=core2) # r0_2 doubles as ref_len2
end

@inline function _segment_induced_v(
    P::Vec3,
    A::Vec3,
    B::Vec3,
    r0::Vec3,
    r0_2::Float64,
    core2::Float64,
    Gamma::Float64,
)::Vec3
    # degenerate (zero-length) segment guard, independent of P.
    if r0_2 < 1e-15 # or 1e-28
        return Vec3(0.0, 0.0, 0.0)
    end

    r1 = P - A
    r2 = P - B

    cr = cross(r1, r2)
    cr2 = dot(cr, cr) + core2 * r0_2

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
    core2::Float64,
)::Vec3
    if ref_len2 < 1e-15 # or 1e-28
        return Vec3(0.0, 0.0, 0.0)
    end

    r1 = P - A
    r1_2 = dot(r1, r1)

    cr = cross(dir, r1)
    cr2 = dot(cr, cr) + core2

    norm1 = sqrt(r1_2 + core2)
    cos_th = dot(dir, r1) / norm1

    factor = Gamma / (4 * pi * cr2) * (1 + cos_th)
    return factor * cr
end

@inline function _ring_induced_v(
    P::Vec3, corners::NTuple{4, Vec3}, inv, Gamma::Float64
)::Vec3
    A, B, C, D = corners
    v1 = _segment_induced_v(
        P, A, B, inv.r0[1], inv.r0_2[1], inv.core2[1], Gamma
    )
    v2 = _segment_induced_v(
        P, B, C, inv.r0[2], inv.r0_2[2], inv.core2[2], Gamma
    )
    v3 = _segment_induced_v(
        P, C, D, inv.r0[3], inv.r0_2[3], inv.core2[3], Gamma
    )
    v4 = _segment_induced_v(
        P, D, A, inv.r0[4], inv.r0_2[4], inv.core2[4], Gamma
    )

    return v1 + v2 + v3 + v4
end

@inline function _horseshoe_induced_v(
    P::Vec3, A::Vec3, B::Vec3, dir::Vec3, inv, Gamma::Float64
)::Vec3
    ref_len2 = inv.r0_2
    return _segment_induced_v(P, A, B, inv.r0, inv.r0_2, inv.core2, Gamma) +
           _semi_infinite_induced_v(P, B, dir, Gamma, ref_len2, inv.core2) -
           _semi_infinite_induced_v(P, A, dir, Gamma, ref_len2, inv.core2)
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
