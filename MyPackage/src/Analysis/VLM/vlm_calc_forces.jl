# Post-processing: turning circulation (gamma) and induced velocities into
# panel forces, then into per-surface loads (_calc_loads) and Trefftz-plane
# induced drag (_calc_trefftz_drag). Consumed by _solve_point in
# vlm_setup.jl; doesn't itself call into vlm_aic.jl.

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
    span_segments = Vector{Matrix{Float64}}(undef, n_surfaces)

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
        row_len = (1 + mirror_xz) * n_span

        forces[i] = Matrix{Vec3}(undef, row_len, n_chord)
        span_segments[i] = Matrix{Float64}(undef, row_len, n_chord)

        for j in range
            rj = rings[j]
            span_vec = rj.corners[2] - rj.corners[1]
            chord_vec = rj.corners[4] - rj.corners[1]

            delta_gamma_c, delta_gamma_s = _calc_delta_gamma(
                gamma, j, start, row_len
            )

            V_total = V_free + Vec3(Vx[j], Vy[j], Vz[j])
            vector = delta_gamma_c * span_vec + delta_gamma_s * chord_vec

            if (j - start + 1) % row_len == 0
                vector += -gamma[j] * (rj.corners[3] - rj.corners[2])
            end

            i_span, i_chord = _j1dto2d(j, start, n_span, Val(mirror_xz))
            span_segments[i][i_span, i_chord] = span_vec[2]
            forces[i][i_span, i_chord] = rho * cross(V_total, vector)
        end
    end

    return forces, span_segments
end

function _calc_delta_gamma(
    gamma::Vector{Float64}, j::Int, start::Int, row_len::Int
)
    delta_gamma_s =
        gamma[j] - (((j - start) % row_len == 0) ? 0.0 : gamma[j - 1])
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
    span_segments::Vector{Matrix{Float64}},
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
        surface_segments = span_segments[i]
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
                span_segment = surface_segments[j, k]

                FX[i] += fx
                FY[i] += fy
                FZ[i] += fz

                FX_dist[i][j] += fx / span_segment
                FY_dist[i][j] += fy / span_segment
                FZ_dist[i][j] += fz / span_segment

                lift = -fx * sa + fz * ca
                drag = fx * ca * cb - fy * sb + fz * sa * cb

                L[i] += lift
                D[i] += drag

                L_dist[i][j] += lift / span_segment
                D_dist[i][j] += drag / span_segment

                rj = rings[i1d]
                i1d += 1
                moment_arm = rj.colpt - CG

                (ml, m, n) = cross(moment_arm, surface_forces[j, k])

                M[i] += m
                Ml[i] += ml
                N[i] += n

                M_dist[i][j] += m / span_segment
                Ml_dist[i][j] += ml / span_segment
                N_dist[i][j] += n / span_segment
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

@inline function _trefftz_plane_basis(V_dir::Vec3)::Tuple{Vec3, Vec3}
    # Any unit seed not parallel to V_dir works; z fails only when V_dir
    # is itself near-vertical, so fall back to x then.
    seed = abs(V_dir[3]) < 0.9 ? Vec3(0.0, 0.0, 1.0) : Vec3(1.0, 0.0, 0.0)
    plane_u_dir = cross(seed, V_dir)
    plane_u_dir = plane_u_dir / sqrt(dot(plane_u_dir, plane_u_dir))
    plane_v_dir = cross(V_dir, plane_u_dir)
    return plane_u_dir, plane_v_dir
end

# Trefftz-plane induced drag. Kept separate from _calc_loads: its inputs
# (circulation and wake geometry) are entirely different from what
# _calc_loads needs (per-panel forces, moment arm, rotation angles), and
# it needs a global filament list assembled across ALL surfaces before any
# per-panel work can start (a wing's trailing vortices induce velocity at
# a downstream tail's TE panels), so it can't share _calc_loads' per-
# surface loop structure either.
function _calc_trefftz_drag(
    gamma::Vector{Float64},
    rings::Vector{VortexRing},
    surfaces::Vector{VLMSurface},
    rho::Float64,
    V_dir::Vec3,
    epsilon2::Float64,
)
    n_surfaces = length(surfaces)

    # --- Assemble the trailing-vortex wake filaments (one per spanwise TE
    # panel boundary, across all surfaces): filament strength is the
    # spanwise difference in TE circulation, with the two wingtips closing
    # the sheet back to zero circulation. ---
    n_wake_filaments = 0
    @inbounds for surface in surfaces
        n_wake_filaments += (1 + surface.mirror_xz) * surface.n_span + 1
    end

    wake_filament_pos = Vector{Vec3}(undef, n_wake_filaments)
    wake_filament_strength = Vector{Float64}(undef, n_wake_filaments)
    wake_filament_ref_len2 = Vector{Float64}(undef, n_wake_filaments)

    filament_idx = 0
    for surface in surfaces
        n_te_panels_surface = (1 + surface.mirror_xz) * surface.n_span
        te_ring_start = surface.range.stop - n_te_panels_surface + 1

        prev_gamma = 0.0
        @inbounds for k in 1:n_te_panels_surface
            ring_idx = te_ring_start + k - 1
            ring_gamma = gamma[ring_idx]
            te_ring = rings[ring_idx]
            te_edge_len2 = dot(
                te_ring.corners[3] - te_ring.corners[4],
                te_ring.corners[3] - te_ring.corners[4],
            )

            filament_idx += 1
            wake_filament_pos[filament_idx] = te_ring.corners[4]
            wake_filament_strength[filament_idx] = ring_gamma - prev_gamma
            wake_filament_ref_len2[filament_idx] = te_edge_len2
            prev_gamma = ring_gamma

            if k == n_te_panels_surface
                filament_idx += 1
                wake_filament_pos[filament_idx] = te_ring.corners[3]
                wake_filament_strength[filament_idx] = -ring_gamma
                wake_filament_ref_len2[filament_idx] = te_edge_len2
            end
        end
    end

    # --- Index the trailing-edge panels themselves (across all surfaces),
    # tracking which surface/spanwise-station each belongs to so the
    # per-panel drag can be scattered back into per-surface results. ---
    n_te_panels_per_surface = Vector{Int}(undef, n_surfaces)
    n_te_panels = 0
    @inbounds for (i, surface) in enumerate(surfaces)
        row_len = (1 + surface.mirror_xz) * surface.n_span
        n_te_panels_per_surface[i] = row_len
        n_te_panels += row_len
    end

    te_panel_ring_idx = Vector{Int}(undef, n_te_panels)
    te_panel_surface_idx = Vector{Int}(undef, n_te_panels)
    te_panel_span_idx = Vector{Int}(undef, n_te_panels)

    t = 0
    for (i, surface) in enumerate(surfaces)
        n_te_panels_surface = n_te_panels_per_surface[i]
        te_ring_start = surface.range.stop - n_te_panels_surface + 1
        @inbounds for k in 1:n_te_panels_surface
            t += 1
            te_panel_ring_idx[t] = te_ring_start + k - 1
            te_panel_surface_idx[t] = i
            te_panel_span_idx[t] = k
        end
    end

    # --- Project onto the Trefftz plane (perpendicular to V_dir) once, up
    # front: self-induction there is a 2D problem, so the O(n_te_panels *
    # n_wake_filaments) sum below is plain scalar 2D arithmetic instead of
    # a Vec3 Biot-Savart call per pair. ---
    plane_u_dir, plane_v_dir = _trefftz_plane_basis(V_dir)

    wake_filament_u = Vector{Float64}(undef, n_wake_filaments)
    wake_filament_v = Vector{Float64}(undef, n_wake_filaments)
    @inbounds for m in 1:n_wake_filaments
        wake_filament_u[m] = dot(wake_filament_pos[m], plane_u_dir)
        wake_filament_v[m] = dot(wake_filament_pos[m], plane_v_dir)
    end

    te_panel_u = Vector{Float64}(undef, n_te_panels)
    te_panel_v = Vector{Float64}(undef, n_te_panels)
    te_span_vec_u = Vector{Float64}(undef, n_te_panels)
    te_span_vec_v = Vector{Float64}(undef, n_te_panels)
    @inbounds for t in 1:n_te_panels
        te_ring = rings[te_panel_ring_idx[t]]
        span_vec = te_ring.corners[3] - te_ring.corners[4]
        midpt = 0.5 * (te_ring.corners[3] + te_ring.corners[4])
        te_panel_u[t] = dot(midpt, plane_u_dir)
        te_panel_v[t] = dot(midpt, plane_v_dir)
        te_span_vec_u[t] = dot(span_vec, plane_u_dir)
        te_span_vec_v[t] = dot(span_vec, plane_v_dir)
    end

    te_panel_drag = Vector{Float64}(undef, n_te_panels)

    # Parallelise over trailing-edge panels (mirrors _assemble_AIC's
    # per-collocation-point threading; the O(n_te_panels *
    # n_wake_filaments) self-induction sum below is the expensive part).
    @threads for t in 1:n_te_panels
        panel_u = te_panel_u[t]
        panel_v = te_panel_v[t]
        induced_u = 0.0
        induced_v = 0.0
        @inbounds for m in 1:n_wake_filaments
            delta_u = panel_u - wake_filament_u[m]
            delta_v = panel_v - wake_filament_v[m]
            dist2 = delta_u^2 + delta_v^2 + epsilon2 * wake_filament_ref_len2[m]
            vortex_coeff = wake_filament_strength[m] / (4 * pi * dist2)
            induced_u += vortex_coeff * delta_v
            induced_v -= vortex_coeff * delta_u
        end
        ring_idx = te_panel_ring_idx[t]
        te_panel_drag[t] =
            rho *
            gamma[ring_idx] *
            (induced_u * te_span_vec_v[t] - induced_v * te_span_vec_u[t])
    end

    D_trefftz = fill(0.0, n_surfaces)
    D_trefftz_dist = Vector{Vector{Float64}}(undef, n_surfaces)

    for i in 1:n_surfaces
        D_trefftz_dist[i] = fill(0.0, n_te_panels_per_surface[i])
    end

    @inbounds for t in 1:n_te_panels
        i = te_panel_surface_idx[t]
        k = te_panel_span_idx[t]
        D_trefftz_dist[i][k] = te_panel_drag[t]
        D_trefftz[i] += te_panel_drag[t]
    end

    return D_trefftz, D_trefftz_dist
end
