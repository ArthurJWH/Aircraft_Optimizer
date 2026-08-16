using ..Geometry
using JSON3

"""
PlaneToJSON.jl

Resolves a `Plane` into concrete, dimensioned 3D loft-section curves for a
solid outer-mold-line (OML), by reusing the *exact* transform math from
`structVLMMesh.jl::_generate_geom` / `_transform_section!` / `_transform_section_v!`,
extended from camber-only to full top+bottom thickness. Writes JSON for the
SolidWorks loft builder.

Geometry now matches your mesher exactly (previous version guessed at the
sweep/dihedral convention -- this one doesn't need to guess):

  - `surface.b` is the FULL span; each Aerosurface generates one semi-span
    from y=0 (root) to y=1 (tip), physically spanning b/2. mirror_xz
    produces the other semi-span by reflecting global Y -> -Y.
  - chord(y), twist(y) [deg] are evaluated directly at each spanwise
    station y in [0,1].
  - sweep(y), dihedral(y) [deg] are LOCAL angles; their tangents are
    integrated (Gauss-Legendre, matching IntegrateGLQ) from 0 to y and
    scaled by b/2 to get sweep_length(y) / dihedral_length(y) -- the
    physical offset of the sw_center-chord reference line.
  - The section (all chordwise points at a given span station) is:
      1. built in local unrotated (x0 = xi*chord, z0 = airfoil(xi)*chord)
         coordinates,
      2. rotated by twist about the pivot x = tw_center*chord (physical,
         i.e. tw_center fraction of THIS station's chord),
      3. translated in x by [sw_center*(root_chord - chord) + sweep_length]
         and in z by dihedral_length, then by surface.pos.
  - For `vertical == true` surfaces, span runs along z instead of y, and
    the "dihedral" offset applies to y instead of z (matches
    `_transform_section_v!`).
  - Extended to thickness: instead of one camber-based z per (station,
    chordwise-index), we carry top_surface AND bottom_surface, blended
    spanwise via linear interpolation between the airfoils at their
    defining `ys` (same technique your mesher uses for camber via
    LinearSpline), sampled at n_span+1 cosine-spaced stations -- so loft
    resolution is decoupled from how many airfoils you defined, exactly
    like your VLM mesh.
  - Also exports three guide curves per surface: leading edge (xi=0),
    trailing edge (xi=1), and the sw_center reference line (the actual
    sweep/quarter-chord line) -- useful as SolidWorks loft guide curves
    so the loft doesn't pinch/twist unexpectedly between stations.
"""

# ---------------------------------------------------------------------
# Gauss-Legendre quadrature (5-point, matches IntegrateGLQ(...; n=5) used
# in structVLMMesh.jl for the sweep/dihedral integrals).
# ---------------------------------------------------------------------
const _GL5_X = [
    -0.9061798459386640,
    -0.5384693101056831,
    0.0,
    0.5384693101056831,
    0.9061798459386640,
]
const _GL5_W = [
    0.2369268850561891,
    0.4786286704993665,
    0.5688888888888889,
    0.4786286704993665,
    0.2369268850561891,
]

function integrate_glq(f, a::Float64, b::Float64)
    a == b && return 0.0
    c1, c2 = (b - a) / 2, (b + a) / 2
    s = 0.0
    @inbounds for i in 1:5
        s += _GL5_W[i] * f(c2 + c1 * _GL5_X[i])
    end
    return s * c1
end

cosine01(n::Int) = (1 .- cos.(range(0; stop=pi, length=n + 1))) ./ 2

# ---------------------------------------------------------------------
# Standalone piecewise-linear spanwise blend between the given airfoils'
# top/bottom surfaces, evaluated at arbitrary y in [0,1]. Mirrors what
# structVLMMesh.jl does with LinearSpline(ys, camber_values) but keeps
# top and bottom separate instead of averaging into camber.
# ---------------------------------------------------------------------
function lerp_at(ys::Vector{Float64}, vals::Vector{Float64}, y::Float64)
    y <= ys[1] && return vals[1]
    y >= ys[end] && return vals[end]
    i = searchsortedlast(ys, y)
    i = clamp(i, 1, length(ys) - 1)
    t = (y - ys[i]) / (ys[i + 1] - ys[i])
    return vals[i] + t * (vals[i + 1] - vals[i])
end

"""
    build_thickness_splines(surface, n_chord)

For each of the n_chord+1 cosine-spaced chordwise stations xi, returns a
closure y -> (top(xi,y), bottom(xi,y)) via spanwise linear interpolation
across the surface's defined airfoils (at their `ys`), exactly mirroring
`splines_z` in structVLMMesh.jl but for top AND bottom.
"""
function build_thickness_splines(surface, n_chord::Int)
    xs = cosine01(n_chord)
    ys_defined = surface.ys
    top_at_xi = Vector{Vector{Float64}}(undef, length(xs))
    bot_at_xi = Vector{Vector{Float64}}(undef, length(xs))
    for (k, xi) in enumerate(xs)
        top_at_xi[k] = [af.top_surface(xi) for af in surface.airfoils]
        bot_at_xi[k] = [af.bottom_surface(xi) for af in surface.airfoils]
    end
    top_fn = (k, y) -> lerp_at(ys_defined, top_at_xi[k], y)
    bot_fn = (k, y) -> lerp_at(ys_defined, bot_at_xi[k], y)
    return xs, top_fn, bot_fn
end

"""
    generate_oml(surface; n_chord=50, n_span=40)

Reproduces `_generate_geom` + `_generate_vertices` + `_transform_section!`
(or `_v!` for vertical surfaces) from structVLMMesh.jl, but for a full
closed top+bottom loop per span station instead of a single camber point,
at n_span+1 cosine-spaced span stations (independent of how many airfoils
were defined -- pass n_span = length(surface.ys)-1 if you want stations
to land exactly on your defined airfoils instead).

Returns (stations, guides) where:
stations :: Vector of (y_frac, chord, twist_deg, points::Vector{[x,y,z]})
points form ONE closed loop: LE -> along top -> TE -> along
bottom (reversed) -> back to LE.
guides   :: Dict("leading_edge"=>pts, "trailing_edge"=>pts,
"reference_line"=>pts) each a Vector{[x,y,z]}, one
point per station, for use as SolidWorks loft guide curves.
"""
function generate_oml(surface; n_chord::Int=50, n_span::Int=40)
    b = surface.b
    ys_defined = surface.ys
    chord_fn, twist_fn = surface.chord, surface.twist
    sweep_fn, dihedral_fn = surface.sweep, surface.dihedral
    tw_center, sw_center = surface.tw_center, surface.sw_center
    pos = surface.pos
    vertical = surface.vertical

    xs, top_fn, bot_fn = build_thickness_splines(surface, n_chord)
    y = cosine01(n_span)

    chords = chord_fn.(y)
    root_chord = chord_fn(0.0)

    sweep_integrand(t) = tand(sweep_fn(t))
    dihedral_integrand(t) = tand(dihedral_fn(t))

    stations = Vector{Dict{String, Any}}()
    le_guide, te_guide, ref_guide = Vector{Any}(), Vector{Any}(), Vector{Any}()

    for (j, yj) in enumerate(y)
        c = chords[j]
        tw = twist_fn(yj)
        sweep_length = b * integrate_glq(sweep_integrand, 0.0, yj) / 2
        dihedral_length = b * integrate_glq(dihedral_integrand, 0.0, yj) / 2
        offset = sw_center * (root_chord - c)
        tw_center_phys = tw_center * c

        cd, sd = cosd(tw), sind(tw)

        # local (x0,z0) pairs, unrotated, chordwise LE->TE (top), TE->LE (bottom)
        top_pts = [(xi * c, top_fn(k, yj) * c) for (k, xi) in enumerate(xs)]
        bot_pts = [(xi * c, bot_fn(k, yj) * c) for (k, xi) in enumerate(xs)]
        loop_local = vcat(top_pts, reverse(bot_pts)[2:end])  # closed, no duplicate TE

        pts3d = Vector{Vector{Float64}}(undef, length(loop_local))
        for (i, (x0, z0)) in enumerate(loop_local)
            x = x0 - tw_center_phys
            new_x = x * cd + z0 * sd + tw_center_phys + offset + sweep_length
            new_off = -x * sd + z0 * cd + dihedral_length   # thickness-plane offset axis
            span_phys = yj * b / 2

            p = zeros(3)
            if !vertical
                p[1] = new_x + pos[1]
                p[2] = span_phys + pos[2]
                p[3] = new_off + pos[3]
            else
                p[1] = new_x + pos[1]
                p[2] = new_off + pos[2]
                p[3] = span_phys + pos[3]
            end
            pts3d[i] = p
        end

        push!(
            stations,
            Dict(
                "y_frac" => yj,
                "chord" => c,
                "twist_deg" => tw,
                "points" => pts3d,
            ),
        )

        # guide points: LE = xi=0 (index 1), TE = xi=1 (index n_chord+1),
        # reference/sweep line = the sw_center-chord point (constant x
        # in local terms == offset+sweep_length term derived analytically,
        # independent of thickness -- use camber-ish midline z at sw_center
        # via linear interpolation between nearest xi samples).
        le = pts3d[1]
        te = pts3d[n_chord + 1]
        # reference line point: x at sw_center fraction, z ~ average of
        # top/bottom at sw_center (put it ON the surface, not floating)
        ref_x0 = sw_center * c
        ref_top = top_fn(
            clamp(searchsortedlast(xs, sw_center), 1, n_chord + 1), yj
        )
        ref_bot = bot_fn(
            clamp(searchsortedlast(xs, sw_center), 1, n_chord + 1), yj
        )
        ref_z0 = 0.5 * (ref_top + ref_bot) * c
        rx = ref_x0 - tw_center_phys
        ref_new_x =
            rx * cd + ref_z0 * sd + tw_center_phys + offset + sweep_length
        ref_new_off = -rx * sd + ref_z0 * cd + dihedral_length
        span_phys = yj * b / 2
        refp = zeros(3)
        if !vertical
            refp[1] = ref_new_x + pos[1]
            refp[2] = span_phys + pos[2]
            refp[3] = ref_new_off + pos[3]
        else
            refp[1] = ref_new_x + pos[1]
            refp[2] = ref_new_off + pos[2]
            refp[3] = span_phys + pos[3]
        end

        push!(le_guide, le)
        push!(te_guide, te)
        push!(ref_guide, refp)
    end

    guides = Dict(
        "leading_edge" => le_guide,
        "trailing_edge" => te_guide,
        "reference_line" => ref_guide,
    )
    return stations, guides
end

function mirror_points(pts::Vector, vertical::Bool)
    # Reflect across the global XZ-plane (Y -> -Y), independent of
    # `vertical`, matching mirror_xz's stated meaning ("mirrored across
    # the xz-plane"). For vertical=true surfaces this still reflects the
    # y-coordinate (their lateral offset axis), e.g. to build both fins
    # of a V-tail.
    return [Any[p[1], -p[2], p[3]] for p in pts]
end

function resolve_surface(surface; n_chord::Int=50, n_span::Int=40)
    stations, guides = generate_oml(surface; n_chord=n_chord, n_span=n_span)

    result = Dict{String, Any}(
        "name" => surface.name,
        "vertical" => surface.vertical,
        "mirror_xz" => surface.mirror_xz,
        "stations" => stations,
        "guides" => guides,
    )

    if surface.mirror_xz
        mstations = [
            Dict(
                "y_frac" => st["y_frac"],
                "chord" => st["chord"],
                "twist_deg" => st["twist_deg"],
                "points" => mirror_points(st["points"], surface.vertical),
            ) for st in stations
        ]
        mguides = Dict(
            k => mirror_points(v, surface.vertical) for (k, v) in guides
        )
        result["mirrored_stations"] = mstations
        result["mirrored_guides"] = mguides
    end

    return result
end

"""
    export_plane_json(plane, filepath; n_chord=50, n_span=40)

Top-level entry point. `n_span` controls loft-station resolution
(independent of how many airfoils you defined -- pass
`n_span = length(surface.ys) - 1` per-surface if you'd rather loft
through exactly your defined airfoil stations and nothing in between).
"""
function export_plane_json(
    plane, filepath::String; n_chord::Int=50, n_span::Int=40
)
    surfaces_out = [
        resolve_surface(s; n_chord=n_chord, n_span=n_span) for
        s in plane.surfaces
    ]
    doc = Dict("units" => "meters", "surfaces" => surfaces_out)
    open(filepath, "w") do io
        JSON3.write(io, doc)
    end
    println(
        "Wrote OML geometry for $(length(surfaces_out)) surface(s) to $filepath"
    )
end

# ------------------------------------------------------------------
# Example usage:
#   include("PlaneToJSON.jl")
#   export_plane_json(my_plane, "plane_oml.json"; n_chord=60, n_span=40)
# ------------------------------------------------------------------
