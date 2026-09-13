# Solver Architecture

Overview of the internal solver pipeline, data flow, caching strategy, and
thread-safety model.

---

## Pipeline Overview

```
[Airfoil .dat files] + [Aerosurface parameters]
              │
              ▼
          [Plane]  ─── Aircraft assembly, CG computation
              │
              ▼
       [VLMGeometry]  ─── Meshing, vortex rings, panel normals
              │
              ▼
        [VLMSetup]  ─── Flow conditions (V∞, ρ, ground, ε²)
              │
              ▼
       [VLMSolver!]  ─── AIC assembly, LU factorisation, force integration
              │
              ▼
        [VLMLoad]  ─── Sparse store of VLMLoadPoint results
         ┌───┴───────────────────────┐
         ▼                           ▼
  [VLMLoadSlice]             [VLMCoefficients]
  Dimensional polars         Non-dimensional polars
         │                           │
         ▼                           ▼
    Plots & tables          [VLMCoefficientsSlice]
                                     │
                                     ▼
                        [VLMStabilityDerivatives]
                         Exact or Fast variant
```

---

## Design Principles

### Geometry Decoupling

[`VLMGeometry`](@ref) is **decoupled from flow conditions**. The mesh, vortex ring
topology, panel normals, and collocation points depend only on the aircraft shape —
not on velocity, density, or angle of attack.

This means a single `VLMGeometry` can be shared across multiple [`VLMSetup`](@ref)
instances at different speeds, altitudes, or ground heights without re-meshing.

### Sparse Deduplication

Solved conditions are stored in a `Dict{Tuple{Float64, Float64}, VLMLoadPoint}` inside
[`VLMLoad`](@ref), keyed by ``(\alpha, \beta)``. Looking up or skipping an
already-solved point is an ``O(1)`` dictionary operation.

### Separation of Loads and Coefficients

Raw dimensional forces and moments (Newtons, Newton-metres) are stored in
[`VLMLoadPoint`](@ref). Non-dimensional coefficients are computed **on demand** via
[`VLMCoefficients`](@ref) using dynamic pressure and reference geometry. This avoids
storing redundant normalised data and keeps the solve path clean.

---

## Caching Strategy

### Ring Geometry Cache (Free Flight)

For free-flight simulations (`ground == false`), the Biot-Savart velocity influence
of the closed vortex rings on each panel collocation point depends only on geometry
and ``\epsilon^2`` — not on flow angles.

This angle-independent velocity field (three ``N \times N`` matrices for
``V_x, V_y, V_z``) is computed **lazily on the first solve** and cached in
`setup.ring_geom`. Subsequent angle sweeps skip this ``O(N^2)`` computation entirely.

!!! note
    The ring geometry cache is invalidated when ground effect is enabled, because the
    image vortex positions depend on ``\alpha``.

### Trailing Wake Cache

The **trailing wake** contribution to the AIC matrix changes with every
``(\alpha, \beta)`` pair because wake filaments realign with the local wind direction.
This portion is always recomputed per solve point.

---

## Thread Safety

[`VLMSolver!`](@ref) is designed to be thread-safe:

- A `ReentrantLock` in `VLMSetup` guards all mutations to `setup.loads` and the ring
  geometry cache.
- AIC matrix assembly and Trefftz drag integration are parallelised via
  `Base.Threads.@threads` within each solve point.

!!! warning
    Thread safety applies to concurrent *reads and writes* of the loads cache within a
    single `VLMSolver!` call. Do not call `VLMSolver!` from multiple threads on the
    same `VLMSetup` simultaneously — use separate `VLMSetup` instances instead.

---

## Data Aggregation Levels

All load and coefficient data is available at three levels:

| Level | Suffix | Type | Description |
|:------|:-------|:-----|:------------|
| Per-surface total | *(none)* | `Vector{Float64}` | One value per surface |
| Spanwise distribution | `_dist` | `Vector{Vector{Float64}}` | Sectional values per span station |
| Whole-aircraft total | `_total` | `Float64` | Sum across all surfaces |

For example, `point.L` is per-surface lift, `point.L_dist` is spanwise lift
distribution, and `point.L_total` is total aircraft lift.

---

## Module Structure

```
MyPackage
├── Utils          — Splines, integration, interpolation
├── PlaneInfo      — Coeffs and Data containers
├── Geometry       — Airfoil, Aerosurface, Plane
├── VLM            — Solver core
│   ├── VLMStructs     — Vec3, VortexRing, VLMSurface, GroundTransform
│   ├── structVLMMesh  — Mesh generation, cosine spacing, section transforms
│   ├── structVLMGeometry — Geometry assembly, vortex ring generation
│   ├── vlm_assemble_aic — AIC matrix and RHS assembly
│   ├── vlm_calc_forces — Kutta-Joukowski force integration
│   ├── structVLMSetup — Flow configuration and ring cache
│   ├── vlm_solver     — VLMSolver! entry points
│   ├── structVLMCoefs — Non-dimensional coefficient computation
│   └── structVLMStab  — Stability derivative computation
├── MyPlots        — Makie mesh plotting, Plots.jl airfoil plotting
└── IO             — .dat file I/O, CAD JSON export
```
