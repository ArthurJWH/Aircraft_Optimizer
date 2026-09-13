# Non-dimensional aerodynamic coefficients (CL, CD, CM, ...) computed from
# a VLMLoadPoint's dimensional forces/moments, using VLMSetup's flow
# condition (V_inf, rho -> dynamic pressure q) and VLMGeometry's reference
# values (Sref, cref, bref). Pure post-processing of numbers already
# computed -- adds no solve-time cost.

"""
    VLMCoefficients

Non-dimensional aerodynamic coefficients corresponding to a single solved [`VLMLoadPoint`](@ref).

Values are normalized by dynamic pressure ``q = \\frac{1}{2} \\rho V_\\infty^2`` and aircraft reference geometry:

  - Force coefficients normalized by ``q S_{\\text{ref}}``:
      + Lift coefficient ``C_L = \\frac{L}{q S_{\\text{ref}}}``
      + Near-field drag coefficient ``C_D = \\frac{D}{q S_{\\text{ref}}}``
      + Trefftz plane induced drag coefficient ``C_{D,\\text{trefftz}} = \\frac{D_{\\text{trefftz}}}{q S_{\\text{ref}}}``
      + Side-force coefficient ``C_Y = \\frac{F_Y}{q S_{\\text{ref}}}``
  - Moment coefficients:
      + Pitching moment coefficient ``C_M = \\frac{M}{q S_{\\text{ref}} c_{\\text{ref}}}`` (normalized by MAC ``c_{\\text{ref}}``)
      + Rolling moment coefficient ``C_{Ml} = \\frac{M_l}{q S_{\\text{ref}} b_{\\text{ref}}}`` (normalized by wingspan ``b_{\\text{ref}}``)
      + Yawing moment coefficient ``C_N = \\frac{N}{q S_{\\text{ref}} b_{\\text{ref}}}`` (normalized by wingspan ``b_{\\text{ref}}``)

# Fields

Each coefficient is available at three levels of aggregation:

  - `CL`, `CL_dist`, `CL_total`: Lift coefficient.
  - `CD`, `CD_dist`, `CD_total`: Near-field drag coefficient.
  - `CD_trefftz`, `CD_trefftz_dist`, `CD_trefftz_total`: Trefftz plane induced drag coefficient.
  - `CY`, `CY_dist`, `CY_total`: Side-force coefficient.
  - `CM`, `CM_dist`, `CM_total`: Pitching moment coefficient about `CG` (pitch-up positive).
  - `CMl`, `CMl_dist`, `CMl_total`: Rolling moment coefficient about `CG` (starboard wing down positive).
  - `CN`, `CN_dist`, `CN_total`: Yawing moment coefficient about `CG` (nose-right positive).
"""
struct VLMCoefficients
    CL::Vector{Float64}
    CL_dist::Vector{Vector{Float64}}
    CL_total::Float64
    CD::Vector{Float64}
    CD_dist::Vector{Vector{Float64}}
    CD_total::Float64
    CD_trefftz::Vector{Float64}
    CD_trefftz_dist::Vector{Vector{Float64}}
    CD_trefftz_total::Float64
    CY::Vector{Float64}
    CY_dist::Vector{Vector{Float64}}
    CY_total::Float64
    CM::Vector{Float64}
    CM_dist::Vector{Vector{Float64}}
    CM_total::Float64
    CMl::Vector{Float64}
    CMl_dist::Vector{Vector{Float64}}
    CMl_total::Float64
    CN::Vector{Float64}
    CN_dist::Vector{Vector{Float64}}
    CN_total::Float64
end

_coef_vec(v::Vector{Float64}, denom::Float64) = v ./ denom
function _coef_dist(vd::Vector{Vector{Float64}}, denom::Float64)
    [d ./ denom for d in vd]
end

# Non-dimensionalizes one VLMLoadPoint's dimensional forces into
# coefficients -- raw values in, packaged tuple out.
function _calc_coefs(
    point::VLMLoadPoint, q::Float64, Sref::Float64, cref::Float64, bref::Float64
)
    qS  = q * Sref
    qSc = qS * cref
    qSb = qS * bref

    return (
        _coef_vec(point.L, qS),
        _coef_dist(point.L_dist, qS),
        point.L_total / qS,
        _coef_vec(point.D, qS),
        _coef_dist(point.D_dist, qS),
        point.D_total / qS,
        _coef_vec(point.D_trefftz, qS),
        _coef_dist(point.D_trefftz_dist, qS),
        point.D_trefftz_total / qS,
        _coef_vec(point.FY, qS),
        _coef_dist(point.FY_dist, qS),
        point.FY_total / qS,
        _coef_vec(point.M, qSc),
        _coef_dist(point.M_dist, qSc),
        point.M_total / qSc,
        _coef_vec(point.Ml, qSb),
        _coef_dist(point.Ml_dist, qSb),
        point.Ml_total / qSb,
        _coef_vec(point.N, qSb),
        _coef_dist(point.N_dist, qSb),
        point.N_total / qSb,
    )
end

"""
    VLMCoefficients(setup::VLMSetup, point::VLMLoadPoint)
    VLMCoefficients(setup::VLMSetup, alpha::Float64, beta::Float64)

Computes non-dimensional aerodynamic coefficients for a given [`VLMLoadPoint`](@ref) or looks up the condition `(alpha, beta)` from `setup.loads`.

Throws an `ArgumentError` if `(alpha, beta)` has not been solved yet (call `VLMSolver!` first).

# Example

```julia
# Single-point coefficients lookup
coefs = VLMCoefficients(setup, 4.0, 0.0)
println(
    \"CL = \",
    coefs.CL_total,
    \", CD = \",
    coefs.CD_total,
    \", CM = \",
    coefs.CM_total,
)
```
"""
function VLMCoefficients(setup::VLMSetup, point::VLMLoadPoint)::VLMCoefficients
    q = 0.5 * setup.rho * setup.V_inf^2
    geometry = setup.geometry
    return VLMCoefficients(
        _calc_coefs(point, q, geometry.Sref, geometry.cref, geometry.bref)...
    )
end

function VLMCoefficients(
    setup::VLMSetup, alpha::Float64, beta::Float64
)::VLMCoefficients
    point = get(setup.loads.points, (alpha, beta), nothing)
    point === nothing && throw(
        ArgumentError(
            "no solved point at (alpha, beta) = ($alpha, $beta) -- call VLMSolver! first",
        ),
    )
    return VLMCoefficients(setup, point)
end

"""
    VLMCoefficientsSlice(setup::VLMSetup; alpha=nothing, beta=nothing)

Extracts an ordered 1D polar slice of non-dimensional aerodynamic coefficients across a sweep of `alpha` or `beta`.

Exactly one of `alpha` or `beta` must be supplied:

  - `beta = b`: Returns a polar curve across sorted angles of attack `alpha`.
  - `alpha = a`: Returns a polar curve across sorted sideslip angles `beta`.

# Returns

A `NamedTuple` containing:

  - Coordinate vector (`alpha` or `beta`).
  - Total coefficients as vectors: `CL_total`, `CD_total`, `CD_trefftz_total`, `CY_total`, `Cm_total`, `Cl_total`, `Cn_total`.
  - Per-surface coefficients as `(n_surfaces, n_points)` matrices: `CL`, `CD`, `CD_trefftz`, `CY`, `CM`, `CMl`, `CN`.

# Example

```julia
# Generate CL vs alpha polar curve
polar = VLMCoefficientsSlice(setup; beta=0.0)
plot(
    polar.alpha,
    polar.CL_total;
    xlabel=\"alpha (deg)\",
    ylabel=\"CL\",
    label=\"Total Lift\",
)
```
"""
function VLMCoefficientsSlice(
    setup::VLMSetup;
    alpha::Union{Float64, Nothing}=nothing,
    beta::Union{Float64, Nothing}=nothing,
)
    slice = VLMLoadSlice(setup.loads; alpha=alpha, beta=beta)
    q = 0.5 * setup.rho * setup.V_inf^2
    geometry = setup.geometry
    qS = q * geometry.Sref
    qSc = qS * geometry.cref
    qSb = qS * geometry.bref

    coefs = (
        CL=slice.L ./ qS,
        CL_total=slice.L_total ./ qS,
        CD=slice.D ./ qS,
        CD_total=slice.D_total ./ qS,
        CD_trefftz=slice.D_trefftz ./ qS,
        CD_trefftz_total=slice.D_trefftz_total ./ qS,
        CY=slice.FY ./ qS,
        CY_total=slice.FY_total ./ qS,
        CM=slice.M ./ qSc,
        Cm_total=slice.M_total ./ qSc,
        CMl=slice.Ml ./ qSb,
        Cl_total=slice.Ml_total ./ qSb,
        CN=slice.N ./ qSb,
        Cn_total=slice.N_total ./ qSb,
    )

    return if alpha === nothing
        merge((alpha=slice.alpha,), coefs)
    else
        merge((beta=slice.beta,), coefs)
    end
end
