# Post-Processing

After running [`VLMSolver!`](@ref), results are stored in `setup.loads`. This page
covers how to extract, inspect, and compare dimensional loads and non-dimensional
coefficients.

---

## Single-Point Inspection

Query coefficients at a specific solved condition with [`VLMCoefficients`](@ref):

```julia
coefs = VLMCoefficients(setup, 4.0, 0.0)   # alpha = 4°, beta = 0°

println("CL = ", coefs.CL_total)
println("CD = ", coefs.CD_total)
println("CM = ", coefs.CM_total)
```

!!! warning
    The requested `(alpha, beta)` must have been previously solved by `VLMSolver!`,
    otherwise an `ArgumentError` is thrown.

---

## Polar Curve Extraction

### Dimensional Loads — [`VLMLoadSlice`](@ref)

Extract a sorted 1D slice of dimensional forces and moments:

```julia
# All points at beta = 0°, sorted by alpha
polar = VLMLoadSlice(setup.loads; beta=0.0)

polar.alpha         # sorted angle of attack vector
    polar.L_total       # total lift, in N
    polar.D_total       # total near-field drag, in N
    polar.M_total       # total pitching moment, in N·m
```

!!! tip
    You can also slice at fixed alpha to get a sideslip sweep:
    ```julia
    lateral = VLMLoadSlice(setup.loads; alpha=4.0)
    lateral.beta        # sorted sideslip vector
    lateral.N_total     # yawing moment, in N·m
    ```

### Non-Dimensional Coefficients — [`VLMCoefficientsSlice`](@ref)

```julia
polar_nd = VLMCoefficientsSlice(setup; beta=0.0)

polar_nd.alpha           # angle of attack vector
polar_nd.CL_total        # lift coefficient
polar_nd.CD_total        # near-field drag coefficient
polar_nd.CD_trefftz_total  # Trefftz induced drag coefficient
polar_nd.Cm_total        # pitching moment coefficient
```

---

## Per-Surface Breakdown

Both `VLMLoadSlice` and `VLMCoefficientsSlice` return per-surface data as
`(n_surfaces x n_points)` matrices, allowing you to isolate contributions
from individual lifting surfaces:

```julia
polar = VLMCoefficientsSlice(setup; beta=0.0)

# CL contribution from surface 1 (e.g. main wing)
polar.CL[1, :]

# CL contribution from surface 2 (e.g. horizontal tail)
polar.CL[2, :]
```

This is useful for understanding how each surface contributes to total lift, drag,
and moments — especially for trim analysis of conventional configurations.

---

## Spanwise Distributions

Individual solved points store sectional load distributions accessible via
[`VLMLoadPoint`](@ref):

```julia
point = setup.loads.points[(4.0, 0.0)]

# Spanwise lift distribution for surface 1
point.L_dist[1]    # Vector{Float64} across span stations

# Spanwise pitching moment distribution for surface 2
point.M_dist[2]
```

!!! note "Distribution indexing"
    `*_dist[i]` is a vector across the spanwise panel stations of surface `i`.
    Values are normalised per unit span segment, so you can integrate them to
    recover the per-surface total.

---

## Coefficient Normalisation

Non-dimensional coefficients are normalised using the **primary surface** (first entry
in `plane.surfaces`) as the reference:

| Coefficient | Formula | Reference quantities |
|:------------|:--------|:---------------------|
| ``C_L, C_D, C_Y`` | ``F / (q \cdot S_\text{ref})`` | ``q = \tfrac{1}{2}\rho V_\infty^2``, ``S_\text{ref}`` |
| ``C_M`` | ``M / (q \cdot S_\text{ref} \cdot c_\text{ref})`` | MAC of primary surface |
| ``C_{Ml}, C_N`` | ``M / (q \cdot S_\text{ref} \cdot b_\text{ref})`` | wingspan of primary surface |

!!! note "Available coefficient fields"
    Each [`VLMCoefficients`](@ref) instance provides:

    | Field | Description |
    |:------|:------------|
    | `CL`, `CL_total`, `CL_dist` | Lift coefficient |
    | `CD`, `CD_total`, `CD_dist` | Near-field drag coefficient |
    | `CD_trefftz`, `CD_trefftz_total`, `CD_trefftz_dist` | Trefftz induced drag coefficient |
    | `CY`, `CY_total`, `CY_dist` | Side-force coefficient |
    | `CM`, `CM_total`, `CM_dist` | Pitching moment coefficient |
    | `CMl`, `CMl_total`, `CMl_dist` | Rolling moment coefficient |
    | `CN`, `CN_total`, `CN_dist` | Yawing moment coefficient |

---

## Comparing Drag Methods

The solver computes drag in two independent ways:

```julia
polar = VLMCoefficientsSlice(setup; beta=0.0)

Plots.plot(polar.alpha, polar.CD_total,
    label = "Near-field (surface integration)", lw = 2,
)
Plots.plot!(polar.alpha, polar.CD_trefftz_total,
    label = "Trefftz plane (far-field)", lw = 2, ls = :dash,
)
```

!!! tip "Which drag to report?"
    **Always prefer Trefftz plane drag** for induced drag reporting. Near-field drag
    is sensitive to chordwise panel density, while Trefftz drag converges with
    spanwise resolution alone. See [VLM Theory](../technical/vlm_theory.md) for the mathematical basis.
