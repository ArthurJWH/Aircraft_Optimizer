# Stability Analysis

Stability derivatives quantify how aerodynamic forces and moments change with small
perturbations in angle of attack (``\alpha``) and sideslip (``\beta``). This page shows
how to compute and interpret them.

!!! tip "Prerequisites"
    You need a configured [`VLMSetup`](@ref) with at least one solved point. See
    [Your First Simulation](first_simulation.md) if you haven't set one up yet.

---

## Computing Stability Derivatives

### Exact Method — [`VLMStabilityDerivatives`](@ref)

Computes central finite differences with full AIC matrix reassembly and wake
realignment at each perturbed point:

```julia
derivs = VLMStabilityDerivatives(setup, 4.0, 0.0)
```

This solves 4 additional conditions:
``(\alpha \pm \Delta\alpha,\; \beta)`` and ``(\alpha,\; \beta \pm \Delta\beta)``.

!!! note
    The perturbed points are cached in `setup.loads`, so subsequent calls at the same
    angles are free.

### Fast Method — [`VLMStabilityDerivativesFast`](@ref)

Reuses a single LU factorisation across all 4 perturbations — roughly **4x faster**:

```julia
derivs_fast = VLMStabilityDerivativesFast(setup, 4.0, 0.0)
```

!!! warning "Linearisation approximation"
    The fast method freezes the AIC matrix at the nominal ``(\alpha, \beta)`` and only
    updates the RHS. Because this solver aligns trailing wake filaments with the local
    wind direction (unlike AVL, where they are fixed to the body axis), this is a real
    linearisation approximation. It is typically accurate within 1–2% near trim.

!!! note
    Perturbed points from the fast method are **not** inserted into `setup.loads`.

---

## Choosing a Method

| | Exact | Fast |
|:---|:---|:---|
| **Cost** | 4 x ``O(N^3)`` factorizations | 1 x ``O(N^3)`` + 4 x ``O(N^2)`` |
| **Wake alignment** | Realigned per perturbation | Frozen at nominal |
| **Accuracy** | Exact numerical derivative | ≈ 1–2% linearisation error |
| **Cache** | Points stored in `setup.loads` | Points discarded |
| **Best for** | Final aerodynamic database, validation | Optimisation loops, trim solvers |

---

## Interpreting the Results

### Longitudinal Derivatives (w.r.t. ``\alpha``)

```julia
println("CL_alpha = ", derivs.CL_alpha, " /deg")
println("CM_alpha = ", derivs.CM_alpha, " /deg")
```

| Derivative | Physical meaning | Stability criterion |
|:-----------|:-----------------|:--------------------|
| ``C_{L\alpha}`` | Lift-curve slope | — |
| ``C_{M\alpha}`` | Pitch stiffness | ``C_{M\alpha} < 0`` for static stability |
| ``C_{D\alpha}`` | Induced drag slope | — |

### Lateral-Directional Derivatives (w.r.t. ``\beta``)

```julia
println("CMl_beta = ", derivs.CMl_beta, " /deg")
println("CN_beta  = ", derivs.CN_beta,  " /deg")
```

| Derivative | Physical meaning | Stability criterion |
|:-----------|:-----------------|:--------------------|
| ``C_{l\beta}`` (`CMl_beta`) | Dihedral effect / roll stability | ``C_{l\beta} < 0`` |
| ``C_{n\beta}`` (`CN_beta`) | Weathercock / directional stability | ``C_{n\beta} > 0`` |
| ``C_{Y\beta}`` (`CY_beta`) | Side-force sensitivity | — |

!!! tip "Field naming"
    In the [`VLMStabilityDerivatives`](@ref) struct, rolling moment derivatives use the
    prefix `CMl` (matching the `CMl` coefficient convention) rather than the lowercase
    ``C_l`` notation common in textbooks.

---

## Customising Step Size

Both methods accept optional `dalpha` and `dbeta` keyword arguments:

```julia
# Finer step for more precise derivatives
derivs = VLMStabilityDerivatives(setup, 4.0, 0.0; dalpha=0.25, dbeta=0.25)
```

!!! note
    The default step sizes are ``\Delta\alpha = 0.5°`` and ``\Delta\beta = 0.5°``.
    Smaller steps improve derivative accuracy but may amplify numerical noise at
    very coarse mesh resolutions.

---

## Full Example

```julia
using MyPackage.VLM

# Assuming setup is already configured and solved at alpha = 4°, beta = 0°

# Exact derivatives
derivs = VLMStabilityDerivatives(setup, 4.0, 0.0)

# Quick check: is the aircraft longitudinally stable?
if derivs.CM_alpha < 0
    println("✓ Longitudinally stable (CM_alpha = ", derivs.CM_alpha, " /deg)")
else
    println("✗ Longitudinally unstable (CM_alpha = ", derivs.CM_alpha, " /deg)")
end

# Quick check: does it have positive directional stability?
if derivs.CN_beta > 0
    println("✓ Directionally stable (CN_beta = ", derivs.CN_beta, " /deg)")
else
    println("✗ Directionally unstable (CN_beta = ", derivs.CN_beta, " /deg)")
end
```
