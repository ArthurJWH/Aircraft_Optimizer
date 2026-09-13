# export_plane_json.jl
#
# Resolves a `Plane` into concrete, dimensioned 3D loft-section curves for a
# solid outer-mold-line (OML), reusing the transform math from
# `structVLMMesh.jl`, extended from camber-only to full top+bottom thickness.
# Writes JSON formatted for CAD / SolidWorks loft builders.

using JSON3

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
    length(ys) <= 1 && return vals[1]
    y <= ys[1] && return vals[1]
    y >= ys[end] && return vals[end]
    i = searchsortedlast(ys, y)
    i = clamp(i, 1, length(ys) - 1)
    t = (y - ys[i]) / (ys[i + 1] - ys[i])
    return vals[i] + t * (vals[i + 1] - vals[i])
end

"""
    build_thickness_splines(surface, n_chord)

Builds spanwise interpolation closures for the upper and lower airfoil thickness contours at `n_chord + 1` chordwise stations ``x_i \\in [0, 1]``.

# Arguments
- `surface::Aerosurface`: The aerodynamic surface definition.
- `n_chord::Int`: Number of chordwise discretization points.

# Returns
- `(xs, top_fn, bot_fn)`: Knot stations `xs`, upper contour interpolator `top_fn(k, y)`, and lower contour interpolator `bot_fn(k, y)`.
"""
function build_thickness_splines(surface, n_chord::Int)
    xs = collect(range(0.0, 1.0, length=n_chord + 1))
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
    generate_oml(surface; n_chord=50, n_span=40, min_thickness_rel=0.002)

Generates 3D loft stations and guide curves for an [`Aerosurface`](@ref MyPackage.Geometry.Aerosurface).

Spanwise sampling combines uniform discretization with all explicit stations defined in `surface.ys`.
Airfoil profile contours run continuously from lower trailing edge to leading edge to upper trailing edge.
Enforces a minimum aerodynamic thickness fraction `min_thickness_rel` for zero-thickness flat-plate profiles to ensure watertight CAD lofting.

# Arguments
- `surface::Aerosurface`: Aerodynamic lifting surface.
- `n_chord::Int`: Number of chordwise divisions (default: `50`).
- `n_span::Int`: Number of spanwise divisions (default: `40`).
- `min_thickness_rel::Float64`: Minimum relative thickness fraction (default: `0.002`).

# Returns
- `(stations, guides)`: Array of station curve dictionaries and guide curve coordinates.
"""
function generate_oml(
    surface;
    n_chord::Int=50,
    n_span::Int=40,
    min_thickness_rel::Float64=0.002,
)
    b = surface.b
    ys_defined = surface.ys
    chord_fn, twist_fn = surface.chord, surface.twist
    sweep_fn, dihedral_fn = surface.sweep, surface.dihedral
    tw_center, sw_center = surface.tw_center, surface.sw_center
    pos = surface.pos
    vertical = surface.vertical

    xs, top_fn, bot_fn = build_thickness_splines(surface, n_chord)
    # Uniform spanwise distribution plus exact user-defined airfoil stations
    y = sort(unique(vcat(collect(range(0.0, 1.0, length=n_span + 1)), ys_defined)))

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

        # Helper to transform 2D section point (x0, z0) into 3D global coordinate
        function to_3d(x0, z0)
            x = x0 - tw_center_phys
            new_x = x * cd + z0 * sd + tw_center_phys + offset + sweep_length
            new_off = -x * sd + z0 * cd + dihedral_length
            span_phys = yj * b / 2
            if !vertical
                return [new_x + pos[1], span_phys + pos[2], new_off + pos[3]]
            else
                return [new_x + pos[1], new_off + pos[2], span_phys + pos[3]]
            end
        end

        # Check maximum section thickness at this station
        max_th = maximum([top_fn(k, yj) - bot_fn(k, yj) for k in 1:(n_chord + 1)])
        apply_min_thickness = max_th < min_thickness_rel

        # Evaluate upper and lower surface profiles
        top_z = Vector{Float64}(undef, n_chord + 1)
        bot_z = Vector{Float64}(undef, n_chord + 1)
        for k in 1:(n_chord + 1)
            t_raw = top_fn(k, yj)
            b_raw = bot_fn(k, yj)
            if apply_min_thickness
                xi = xs[k]
                camber_k = 0.5 * (t_raw + b_raw)
                # Aerodynamic thickness envelope: zero gap at LE (xi=0), max thickness near mid-chord, thin finite TE (xi=1)
                shape_factor = 2.0 * sqrt(clamp(xi, 0.0, 1.0)) * (1.0 - 0.9 * xi)
                t_eff = max(t_raw - b_raw, min_thickness_rel * shape_factor)
                top_z[k] = (camber_k + 0.5 * t_eff) * c
                bot_z[k] = (camber_k - 0.5 * t_eff) * c
            else
                top_z[k] = t_raw * c
                bot_z[k] = b_raw * c
            end
        end

        # Lower surface points from TE (xi=1) down to LE (xi=0)
        pts_bot = [to_3d(xs[k] * c, bot_z[k]) for k in (n_chord + 1):-1:1]
        # Upper surface points from LE (xi=0) up to TE (xi=1)
        pts_top = [to_3d(xs[k] * c, top_z[k]) for k in 1:(n_chord + 1)]

        # Single continuous contour: Lower TE -> LE -> Upper TE (without duplicating LE)
        pts_contour = vcat(pts_bot, pts_top[2:end])

        push!(
            stations,
            Dict(
                "y_frac" => yj,
                "chord" => c,
                "twist_deg" => tw,
                "points" => pts_contour,
                "top_points" => pts_top,
                "bottom_points" => pts_bot,
            ),
        )

        # Guide points:
        # LE guide: xi=0
        le = pts_top[1]
        # TE guide: midpoint between top and bottom TE
        te_bot = pts_bot[1]
        te_top = pts_top[end]
        te = [0.5 * (te_top[1] + te_bot[1]), 0.5 * (te_top[2] + te_bot[2]), 0.5 * (te_top[3] + te_bot[3])]

        # Reference line guide at sw_center chord
        ref_x0 = sw_center * c
        idx = clamp(searchsortedlast(xs, sw_center), 1, n_chord + 1)
        ref_z0 = 0.5 * (top_z[idx] + bot_z[idx])
        refp = to_3d(ref_x0, ref_z0)

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

function resolve_surface(
    surface;
    n_chord::Int=50,
    n_span::Int=40,
    min_thickness_rel::Float64=0.002,
)
    stations, guides = generate_oml(
        surface;
        n_chord=n_chord,
        n_span=n_span,
        min_thickness_rel=min_thickness_rel,
    )

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
                "top_points" => mirror_points(st["top_points"], surface.vertical),
                "bottom_points" => mirror_points(st["bottom_points"], surface.vertical),
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
    export_plane_json(plane, filepath; n_chord=50, n_span=40, min_thickness_rel=0.002)

Exports the outer-mold-line (OML) loft sections and guide curves of all surfaces in a [`Plane`](@ref MyPackage.Geometry.Plane) to a JSON file.

The resulting JSON schema is structured for automated SolidWorks / CAD macro loft construction scripts.

# Arguments
- `plane::Plane`: The aircraft definition to export.
- `filepath::String`: Output path for the `.json` file.
- `n_chord::Int`: Number of chordwise discretization points per surface (default: `50`).
- `n_span::Int`: Number of spanwise loft stations per surface (default: `40`).
- `min_thickness_rel::Float64`: Minimum relative thickness fraction enforced for zero-thickness flat-plate airfoils (default: `0.002`).

# Example
```julia
using MyPackage.IO
export_plane_json(my_plane, "plane_oml.json"; n_chord=50, n_span=40)
```
"""
function export_plane_json(
    plane,
    filepath::String;
    n_chord::Int=50,
    n_span::Int=40,
    min_thickness_rel::Float64=0.002,
)
    surfaces_out = [
        resolve_surface(
            s;
            n_chord=n_chord,
            n_span=n_span,
            min_thickness_rel=min_thickness_rel,
        ) for s in plane.surfaces
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
