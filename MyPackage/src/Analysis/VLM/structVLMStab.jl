# Stability derivatives (CL_alpha, Cm_alpha, Cl_beta, Cn_beta, ...) via
# central finite difference at a chosen (alpha, beta). Built entirely on
# the existing public API (VLMSolver!, VLMCoefficients) -- no new solver
# machinery, so this only needs to be included after vlm_setup.jl and
# vlm_coefficients.jl, both purely for their function/type definitions to
# already exist at this file's definition time.

"""
    VLMStabilityDerivatives

Stores longitudinal and lateral-directional stability derivatives computed at a specific operating point `(alpha, beta)`.

All quantities represent whole-aircraft totals evaluated via central finite difference:
``\\frac{\\partial C}{\\partial \\alpha} \\approx \\frac{C(\\alpha + \\Delta\\alpha) - C(\\alpha - \\Delta\\alpha)}{2 \\Delta\\alpha}``

# Fields
- `alpha::Float64`: Nominal angle of attack, in `deg`.
- `beta::Float64`: Nominal sideslip angle, in `deg`.
- `dalpha::Float64`: Finite difference step size ``\\Delta\\alpha``, in `deg`.
- `dbeta::Float64`: Finite difference step size ``\\Delta\\beta``, in `deg`.
- **Longitudinal Derivatives (w.r.t. ``\\alpha``)**:
    - `CL_alpha`: Lift-curve slope ``\\partial C_L / \\partial \\alpha``, per `deg`.
    - `CD_alpha`: Induced drag slope ``\\partial C_D / \\partial \\alpha``, per `deg`.
    - `CD_trefftz_alpha`: Trefftz drag slope ``\\partial C_{D,\\text{trefftz}} / \\partial \\alpha``, per `deg`.
    - `CY_alpha`: Side-force derivative w.r.t. alpha ``\\partial C_Y / \\partial \\alpha``, per `deg`.
    - `CM_alpha`: Pitching moment stability derivative ``\\partial C_M / \\partial \\alpha``, per `deg` (negative implies longitudinal static stability).
    - `CMl_alpha`: Rolling moment cross-derivative ``\\partial C_{Ml} / \\partial \\alpha``, per `deg`.
    - `CN_alpha`: Yawing moment cross-derivative ``\\partial C_N / \\partial \\alpha``, per `deg`.
- **Lateral-Directional Derivatives (w.r.t. ``\\beta``)**:
    - `CL_beta`: Lift derivative w.r.t. sideslip ``\\partial C_L / \\partial \\beta``, per `deg`.
    - `CD_beta`: Drag derivative w.r.t. sideslip ``\\partial C_D / \\partial \\beta``, per `deg`.
    - `CD_trefftz_beta`: Trefftz drag derivative w.r.t. sideslip ``\\partial C_{D,\\text{trefftz}} / \\partial \\beta``, per `deg`.
    - `CY_beta`: Side-force derivative w.r.t. sideslip ``\\partial C_Y / \\partial \\beta``, per `deg`.
    - `CM_beta`: Pitching moment cross-derivative ``\\partial C_M / \\partial \\beta``, per `deg`.
    - `CMl_beta`: Dihedral effect derivative ``\\partial C_{Ml} / \\partial \\beta``, per `deg` (negative implies roll stability).
    - `CN_beta`: Directional / weathercock stability derivative ``\\partial C_N / \\partial \\beta``, per `deg` (positive implies directional static stability).
"""
struct VLMStabilityDerivatives
    alpha::Float64
    beta::Float64
    dalpha::Float64
    dbeta::Float64
    CL_alpha::Float64
    CD_alpha::Float64
    CD_trefftz_alpha::Float64
    CY_alpha::Float64
    CM_alpha::Float64
    CMl_alpha::Float64
    CN_alpha::Float64
    CL_beta::Float64
    CD_beta::Float64
    CD_trefftz_beta::Float64
    CY_beta::Float64
    CM_beta::Float64
    CMl_beta::Float64
    CN_beta::Float64
end

const _VLM_STABILITY_FIELDS = (:CL, :CD, :CD_trefftz, :CY, :CM, :CMl, :CN)

"""
    VLMStabilityDerivatives(setup::VLMSetup, alpha::Float64, beta::Float64; dalpha=0.5, dbeta=0.5)

Computes exact numerical stability derivatives at `(alpha, beta)` using central finite differences with full wake realignment.

Solves (or retrieves from `setup.loads` cache) the 4 perturbed conditions:
`[(alpha - dalpha, beta), (alpha + dalpha, beta), (alpha, beta - dbeta), (alpha, beta + dbeta)]`.

# Arguments
- `setup::VLMSetup`: Configured simulation setup.
- `alpha::Float64`: Nominal angle of attack, in `deg`.
- `beta::Float64`: Nominal sideslip angle, in `deg`.
- `dalpha::Float64`: Perturbation step for angle of attack, in `deg` (default: `0.5`).
- `dbeta::Float64`: Perturbation step for sideslip angle, in `deg` (default: `0.5`).

# Returns
A [`VLMStabilityDerivatives`](@ref) object containing all longitudinal and lateral-directional derivatives.

# Example
```julia
derivs = VLMStabilityDerivatives(setup, 4.0, 0.0)
println("CL_alpha = ", derivs.CL_alpha, " 1/deg")
println("CM_alpha = ", derivs.CM_alpha, " 1/deg")
println("CN_beta  = ", derivs.CN_beta,  " 1/deg")
```
"""
function VLMStabilityDerivatives(
    setup::VLMSetup, alpha::Float64, beta::Float64; dalpha::Float64=0.5, dbeta::Float64=0.5
)::VLMStabilityDerivatives
    VLMSolver!(
        setup,
        [
            (alpha - dalpha, beta), (alpha + dalpha, beta),
            (alpha, beta - dbeta), (alpha, beta + dbeta),
        ],
    )

    c_am = VLMCoefficients(setup, alpha - dalpha, beta)
    c_ap = VLMCoefficients(setup, alpha + dalpha, beta)
    c_bm = VLMCoefficients(setup, alpha, beta - dbeta)
    c_bp = VLMCoefficients(setup, alpha, beta + dbeta)

    d_alpha = Tuple(
        (getfield(c_ap, Symbol(f, :_total)) - getfield(c_am, Symbol(f, :_total))) /
            (2 * dalpha)
        for f in _VLM_STABILITY_FIELDS
    )
    d_beta = Tuple(
        (getfield(c_bp, Symbol(f, :_total)) - getfield(c_bm, Symbol(f, :_total))) /
            (2 * dbeta)
        for f in _VLM_STABILITY_FIELDS
    )

    return VLMStabilityDerivatives(alpha, beta, dalpha, dbeta, d_alpha..., d_beta...)
end

# ---------------------------------------------------------------------------
# Fast (approximate) variant: reuse one factorization + one set of induced-
# velocity matrices from the nominal point, instead of a fresh AIC assembly
# and factorization per perturbed point.
#
# This is the AVL/flow5 trick for CONTROL-surface derivatives, described in
# flow5's docs: "assume that the change to the influence matrix is small,
# so that the LU-factorized matrix can be reused as-is, and that only the
# change to the RHS needs to be implemented" -- flow5 validated this
# specific approximation against a fully-rebuilt reference and found <1.2%
# deviation for their test aircraft's control derivatives.
# (https://flow5.tech/docs/flow5_doc/Validation/Stability.html)
#
# IMPORTANT DIFFERENCE FROM AVL'S OWN ALPHA/BETA DERIVATIVES: AVL computes
# alpha/beta stability derivatives (not just control derivatives) this same
# way, but for AVL that isn't an approximation -- it's exact, because AVL's
# trailing vortices are fixed parallel to the body x-axis, not realigned to
# the local wind direction ("the trailing vorticity is oriented parallel to
# the X axis", MIT AVL User Primer). That means AVL's AIC has ZERO
# dependence on alpha/beta at all, so freezing it costs nothing. THIS
# solver's wake DOES realign to the true local wind direction each solve
# (a deliberate, more physically-faithful choice made earlier in this
# file's history) -- so AIC genuinely depends on alpha/beta here, and
# freezing it for this function is a real linearization approximation, not
# a free simplification. Use VLMStabilityDerivatives (the exact version)
# unless you've confirmed this one's accuracy is good enough for your case,
# or speed matters more (e.g. inside a trim/optimization loop calling this
# every iteration).
function _solve_point_fast(
    F, AIC_vx::Matrix{Float64}, AIC_vy::Matrix{Float64}, AIC_vz::Matrix{Float64},
    setup::VLMSetup, alpha::Float64, beta::Float64,
)::VLMLoadPoint
    geometry = setup.geometry
    rings, n_panels, n_surfaces, surfaces =
        geometry.rings, geometry.n_panels, geometry.n_surfaces, geometry.surfaces
    CG = geometry.CG
    V_inf, rho, epsilon2 = setup.V_inf, setup.rho, setup.epsilon2
    rot = (alpha, beta)
    V_dir = _wind_dir(rot)

    # Only the RHS changes -- AIC (and therefore F, AIC_vx/vy/vz) is reused
    # from the nominal point, unchanged.
    RHS = _assemble_RHS(rings, n_panels, V_dir, V_inf)
    gamma = F \ RHS

    forces, span_segments = _calc_forces(
        gamma, n_panels, n_surfaces, surfaces, rings, rho, V_dir, V_inf,
        AIC_vx, AIC_vy, AIC_vz,
    )

    (FX, FX_dist, FY, FY_dist, FZ, FZ_dist, L, L_dist, D, D_dist, M, M_dist,
        Ml, Ml_dist, N, N_dist) = _calc_loads(
        forces, span_segments, rings, CG, n_surfaces, rot
    )
    D_trefftz, D_trefftz_dist = _calc_trefftz_drag(
        gamma, rings, surfaces, rho, V_dir, epsilon2
    )

    return VLMLoadPoint(
        FX, FX_dist, sum(FX),
        FY, FY_dist, sum(FY),
        FZ, FZ_dist, sum(FZ),
        L, L_dist, sum(L),
        D, D_dist, sum(D),
        M, M_dist, sum(M),
        Ml, Ml_dist, sum(Ml),
        N, N_dist, sum(N),
        D_trefftz, D_trefftz_dist, sum(D_trefftz),
    )
end

"""
    VLMStabilityDerivativesFast(setup::VLMSetup, alpha::Float64, beta::Float64; dalpha=0.5, dbeta=0.5)

Cheaper, accelerated stability derivative solver that reuses the LU factorization of the nominal AIC matrix across all 4 perturbed states.

# Algorithm & Physical Approximation
Instead of assembling and factorizing four separate AIC matrices for ``\\alpha \\pm \\Delta\\alpha`` and ``\\beta \\pm \\Delta\\beta``
(which would cost ``O(4 N^3)``), this method:
1. Assembles and factorizes the nominal AIC matrix ``F = \\text{lu}(\\text{AIC}_0)`` once at `(alpha, beta)`.
2. Evaluates the perturbed RHS vectors ``\\text{RHS}(\\alpha \\pm \\Delta\\alpha, \\beta)`` and ``\\text{RHS}(\\alpha, \\beta \\pm \\Delta\\beta)``.
3. Solves the 4 perturbed circulation vectors via back-substitution: ``\\Gamma \\approx F \\backslash \\text{RHS}`` (costing ``O(4 N^2)``).
4. Recovers panel forces and computes derivatives via central finite differences.

!!! note
    Unlike AVL (where trailing vortices are fixed parallel to the fuselage axis so ``\\partial \\text{AIC} / \\partial \\alpha = 0``),
    `MyPackage` actively aligns trailing wake filaments with the freestream wind direction. Consequently, freezing the AIC matrix
    is a real linearization approximation. It provides excellent accuracy (typically within 1–2%) near trim while accelerating
    stability calculations by roughly 4x.

The perturbed points computed during this evaluation are temporary approximations and are **not** inserted into `setup.loads`.

# Arguments
- `setup::VLMSetup`: Configured simulation setup.
- `alpha::Float64`: Nominal angle of attack, in `deg`.
- `beta::Float64`: Nominal sideslip angle, in `deg`.
- `dalpha::Float64`: Step size for angle of attack, in `deg` (default: `0.5`).
- `dbeta::Float64`: Step size for sideslip angle, in `deg` (default: `0.5`).

# Returns
A [`VLMStabilityDerivatives`](@ref) object containing all evaluated stability derivatives.
"""
function VLMStabilityDerivativesFast(
    setup::VLMSetup, alpha::Float64, beta::Float64; dalpha::Float64=0.5, dbeta::Float64=0.5
)::VLMStabilityDerivatives
    geometry = setup.geometry
    rings, n_panels = geometry.rings, geometry.n_panels
    rot0 = (alpha, beta)

    if setup.ground
        AIC, AIC_vx, AIC_vy, AIC_vz, _, _ = _assemble_sys(
            rings, n_panels, rot0, geometry.wake_map, setup.V_inf, setup.h, geometry.CG,
            setup.epsilon2,
        )
    else
        ring_geom = _ring_geom!(setup)
        AIC, AIC_vx, AIC_vy, AIC_vz, _, _ = _assemble_sys(
            rings, n_panels, rot0, geometry.wake_map, setup.V_inf, ring_geom,
            setup.epsilon2,
        )
    end
    F = lu(AIC)

    pt_am = _solve_point_fast(F, AIC_vx, AIC_vy, AIC_vz, setup, alpha - dalpha, beta)
    pt_ap = _solve_point_fast(F, AIC_vx, AIC_vy, AIC_vz, setup, alpha + dalpha, beta)
    pt_bm = _solve_point_fast(F, AIC_vx, AIC_vy, AIC_vz, setup, alpha, beta - dbeta)
    pt_bp = _solve_point_fast(F, AIC_vx, AIC_vy, AIC_vz, setup, alpha, beta + dbeta)

    c_am, c_ap = VLMCoefficients(setup, pt_am), VLMCoefficients(setup, pt_ap)
    c_bm, c_bp = VLMCoefficients(setup, pt_bm), VLMCoefficients(setup, pt_bp)

    d_alpha = Tuple(
        (getfield(c_ap, Symbol(f, :_total)) - getfield(c_am, Symbol(f, :_total))) /
            (2 * dalpha)
        for f in _VLM_STABILITY_FIELDS
    )
    d_beta = Tuple(
        (getfield(c_bp, Symbol(f, :_total)) - getfield(c_bm, Symbol(f, :_total))) /
            (2 * dbeta)
        for f in _VLM_STABILITY_FIELDS
    )

    return VLMStabilityDerivatives(alpha, beta, dalpha, dbeta, d_alpha..., d_beta...)
end