# The solver's results API: VLMPolarPoint (one solved condition),
# VLMPolar (the sparse store of them), and VLMLoadSlice (pulling a sorted
# alpha- or beta-polar out of that store for plotting). No dependency on
# mesh/AIC internals -- this file only ever sees plain numbers.

"""
    VLMLoadPoint

Dimensional aerodynamic loads and moments for exactly one solved `(alpha, beta)` flow condition.

Data is provided at three distinct levels of aggregation:
- `*_total::Float64`: The whole-aircraft scalar value obtained by summing across all surfaces.
- `*::Vector{Float64}`: Per-surface integrated totals of length `n_surfaces`.
- `*_dist::Vector{Vector{Float64}}`: Per-surface spanwise sectional distributions. `*_dist[i]` is a vector across the spanwise panel stations of surface `i` (normalized per unit span segment).

# Fields
- `FX`, `FX_dist`, `FX_total`: Body-axis aerodynamic force in the X direction (aft positive), in `N`.
- `FY`, `FY_dist`, `FY_total`: Body-axis aerodynamic side force in the Y direction (starboard positive), in `N`.
- `FZ`, `FZ_dist`, `FZ_total`: Body-axis aerodynamic force in the Z direction (upward positive), in `N`.
- `L`, `L_dist`, `L_total`: Wind-axis lift force perpendicular to the freestream direction, in `N`.
- `D`, `D_dist`, `D_total`: Wind-axis near-field drag force parallel to the freestream direction, in `N`.
- `M`, `M_dist`, `M_total`: Pitching moment about the aircraft center of gravity `CG` (pitch-up positive), in `N·m`.
- `Ml`, `Ml_dist`, `Ml_total`: Rolling moment about `CG` (right wing down positive), in `N·m`.
- `N`, `N_dist`, `N_total`: Yawing moment about `CG` (nose-right positive), in `N·m`.
- `D_trefftz`, `D_trefftz_dist`, `D_trefftz_total`: Induced drag evaluated via far-field Trefftz plane integration, in `N`.
"""
struct VLMLoadPoint
    FX::Vector{Float64}
    FX_dist::Vector{Vector{Float64}}
    FX_total::Float64
    FY::Vector{Float64}
    FY_dist::Vector{Vector{Float64}}
    FY_total::Float64
    FZ::Vector{Float64}
    FZ_dist::Vector{Vector{Float64}}
    FZ_total::Float64
    L::Vector{Float64}
    L_dist::Vector{Vector{Float64}}
    L_total::Float64
    D::Vector{Float64}
    D_dist::Vector{Vector{Float64}}
    D_total::Float64
    M::Vector{Float64}
    M_dist::Vector{Vector{Float64}}
    M_total::Float64
    Ml::Vector{Float64}
    Ml_dist::Vector{Vector{Float64}}
    Ml_total::Float64
    N::Vector{Float64}
    N_dist::Vector{Vector{Float64}}
    N_total::Float64
    D_trefftz::Vector{Float64}
    D_trefftz_dist::Vector{Vector{Float64}}
    D_trefftz_total::Float64
end

"""
    VLMLoad

Sparse repository of solved `(alpha, beta)` conditions stored within a `VLMSetup`.

# Fields
- `points::Dict{Tuple{Float64, Float64}, VLMLoadPoint}`: Dictionary mapping each solved `(alpha, beta)` pair (in `deg`) to its corresponding [`VLMLoadPoint`](@ref).

New conditions are inserted in ``O(1)`` time without reallocating existing solutions. Use [`VLMLoadSlice`](@ref) to extract ordered 1D polar curves for plotting.
"""
struct VLMLoad
    points::Dict{Tuple{Float64, Float64}, VLMLoadPoint}
end
VLMLoad() = VLMLoad(Dict{Tuple{Float64, Float64}, VLMLoadPoint}())

# Per-surface breakdown -- these hcat into (n_surfaces, n_points) matrices.
const _VLM_LOAD_FIELDS = (:FX, :FY, :FZ, :L, :D, :M, :Ml, :N, :D_trefftz)
# Whole-aircraft scalars -- these collect into (n_points,) vectors.
const _VLM_LOAD_TOTAL_FIELDS = (
    :FX_total, :FY_total, :FZ_total, :L_total, :D_total, :M_total, :Ml_total,
    :N_total, :D_trefftz_total,
)

function _slice_matrices(matches)
    surf = NamedTuple{_VLM_LOAD_FIELDS}(
        Tuple(reduce(hcat, getfield(pt, f) for (_, pt) in matches) for f in _VLM_LOAD_FIELDS)
    )
    tot = NamedTuple{_VLM_LOAD_TOTAL_FIELDS}(
        Tuple([getfield(pt, f) for (_, pt) in matches] for f in _VLM_LOAD_TOTAL_FIELDS)
    )
    return merge(surf, tot)
end

"""
    VLMLoadSlice(loads::VLMLoad; alpha=nothing, beta=nothing)

Extracts an ordered 1D polar slice from a [`VLMLoad`](@ref) storage container at a fixed `alpha` or `beta`.

Exactly one of `alpha` or `beta` must be provided as a keyword argument:
- `beta = b`: Returns points sorted by ascending angle of attack `alpha`.
- `alpha = a`: Returns points sorted by ascending sideslip angle `beta`.

# Returns
A `NamedTuple` containing:
- Coordinate vector (`alpha` or `beta`).
- Whole-aircraft loads as vectors: `L_total`, `D_total`, `M_total`, `Ml_total`, `N_total`, `FX_total`, `FY_total`, `FZ_total`, `D_trefftz_total`.
- Per-surface loads as `(n_surfaces, n_points)` matrices: `L`, `D`, `M`, `Ml`, `N`, `FX`, `FY`, `FZ`, `D_trefftz`.

# Example
```julia
# Extract an alpha polar at beta = 0.0 deg
polar = VLMLoadSlice(setup.loads; beta = 0.0)
plot(polar.alpha, polar.L_total, xlabel="alpha (deg)", ylabel="Lift (N)")
```
"""
function VLMLoadSlice(
    loads::VLMLoad;
    alpha::Union{Float64, Nothing}=nothing,
    beta::Union{Float64, Nothing}=nothing,
)
    (alpha === nothing) == (beta === nothing) &&
        throw(ArgumentError("give exactly one of alpha or beta"))

    if beta !== nothing
        matches = sort([(a, pt) for ((a, b), pt) in loads.points if b == beta]; by=first)
        isempty(matches) && throw(ArgumentError("no stored points at beta = $beta"))
        return merge((alpha=first.(matches),), _slice_matrices(matches))
    else
        matches = sort([(b, pt) for ((a, b), pt) in loads.points if a == alpha]; by=first)
        isempty(matches) && throw(ArgumentError("no stored points at alpha = $alpha"))
        return merge((beta=first.(matches),), _slice_matrices(matches))
    end
end

