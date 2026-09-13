# Your First Simulation

This guide walks you through a complete VLM simulation — from loading airfoil profiles to plotting aerodynamic results.

!!! tip "Before you begin"
    Make sure you have completed the [Getting Started](../onboarding/getting_started.md) setup. All code below
    assumes the Julia REPL is running in the project environment.

## Step 1 — Load Airfoil Profiles

Airfoils are loaded from standard `.dat` coordinate files (Selig or Lednicer format).

```julia
using MyPackage.IO
using MyPackage.Geometry

# Load your airfoil from a .dat file
airfoil = airfoil_from_dat("assets/airfoils/NACA4412/NACA4412.dat")
```

The returned [`Airfoil`](@ref MyPackage.Geometry.Airfoil) contains spline closures for the upper surface, lower surface,
and mean camber line, each defined over ``x \in [0, 1]``.

!!! note
    You can visualise any loaded airfoil with [`plot_airfoil`](@ref MyPackage.MyPlots.plot_airfoil):
    ```julia
    using MyPackage.MyPlots
    plot_airfoil(airfoil)
    ```

---

## Step 2 — Define Lifting Surfaces

Each lifting surface (wing, tail, fin, canard) is represented by an [`Aerosurface`](@ref).
Geometric distributions are defined as functions of the normalised semi-span fraction ``y \in [0, 1]``.

```julia
# --- Main Wing ---
wing = Aerosurface(
    name       = "MainWing",
    airfoils   = [airfoil, airfoil],       # root and tip airfoils
    b          = 10.0,                      # full wingspan [m]
    chord      = y -> 1.5 * (1 - 0.4y),    # linear taper
    twist      = y -> -2.0y,               # 2° washout at tip
    sweep      = y -> 5.0,                 # constant 5° sweep
    dihedral   = y -> 3.0,                 # constant 3° dihedral
    sw_center  = 0.25,                     # quarter-chord sweep reference
    tw_center  = 0.25,                     # quarter-chord twist axis
    mirror_xz  = true,                     # mirror for port side
)
```

!!! tip "Computed properties"
    On construction, `Aerosurface` automatically computes reference quantities via
    numerical quadrature ([`IntegrateGLQ`](@ref)):
    planform area ``S``, aspect ratio ``AR``, mean geometric chord ``MGC``,
    and mean aerodynamic chord ``MAC``.

!!! note "Adding more surfaces"
    A tail, canard, or vertical fin follows the same pattern — just create additional
    `Aerosurface` objects with their own geometry. Set `vertical = true` for vertical
    surfaces (e.g. fins) and `mirror_xz = false` for single-sided surfaces (and vertical surfaces).

    ```julia
    htail = Aerosurface(
        name     = "HTail",
        airfoils = [plain, plain],
        b        = 3.0,
        chord    = y -> 0.6,
        pos      = (4.0, 0.0, 0.3),    # offset aft and up from the wing
    )
    ```

---

## Step 3 — Assemble the Aircraft

Combine all surfaces into a [`Plane`](@ref):

```julia
plane = Plane([wing])                              # single surface
# plane = Plane([wing, htail])                     # multi-surface
# plane = Plane([wing, htail]; CG=(0.5, 0.0, 0.0))  # explicit CG
```

!!! tip "Default center of gravity"
    If `CG` is not specified, it defaults to the quarter-chord of the primary wing's
    mean aerodynamic chord: `(pos[1] + 0.25 x MAC, pos[2], pos[3])`.

---

## Step 4 — Discretise the Geometry

Convert the continuous geometry into a VLM panel mesh with [`VLMGeometry`](@ref):

```julia
using MyPackage.VLM

# (chordwise panels, spanwise panels per semi-span) for each surface
geom = VLMGeometry(plane, [(20, 30)])
```

Panels are distributed using **cosine clustering**, concentrating resolution at the
leading edge, trailing edge, root, and wingtips where circulation gradients are steepest.

!!! warning "Panel count selection"
    More panels improve accuracy but increase solve time as ``O(N^2)`` for assembly
    and ``O(N^3)`` for factorisation. A good starting point for a single wing is
    `(15–20, 30–50)`. Always verify convergence by comparing results at different
    resolutions.

---

## Step 5 — Visualise the Mesh

Before solving, inspect the discretisation visually:

```julia
using GLMakie
GLMakie.activate!()
using MyPackage.MyPlots

fig = plot_mesh(geom.meshes)
```

[`plot_mesh`](@ref MyPackage.MyPlots.plot_mesh) renders an interactive 3D view in a separate window of all panels (including mirrored
port-side surfaces) with wireframe outlines.

!!! tip "What to check"
    - Panel aspect ratios look reasonable (not excessively stretched)
    - Camber and twist are oriented correctly
    - Mirrored surfaces appear where expected
    - Surface offsets (`pos`) are correct for multi-surface configurations

---

## Step 6 — Set Up the Flow Environment

Create a [`VLMSetup`](@ref) binding the geometry to specific flight conditions:

```julia
setup = VLMSetup(geom, 25.0; rho=1.225)
```

| Parameter   | Description                          | Default   |
|:------------|:-------------------------------------|:----------|
| `V_inf`     | Freestream airspeed, in `m/s`             | —         |
| `rho`       | Atmospheric density, in `kg/m³`           | `1.225`   |
| `ground`    | Enable ground effect                 | `false`   |
| `h`         | Height above ground, in `m`               | `0.0`     |
| `epsilon2`  | Vortex core regularisation           | `1e-10`   |

!!! note "Geometry reuse"
    A single `VLMGeometry` can feed multiple `VLMSetup` instances at different speeds
    or altitudes — the mesh only needs to be generated once.

---

## Step 7 — Run the Solver

Solve across a range of angles of attack with [`VLMSolver!`](@ref):

```julia
VLMSolver!(setup, collect(-4.0:2.0:12.0), [0.0])
```

This evaluates every ``(\alpha, \beta)`` combination in the Cartesian product of the
supplied vectors.

!!! note "Mutating function"
    `VLMSolver!` (note the `!`) stores results in `setup.loads` in-place. Re-calling
    with different angles appends new points without re-solving existing ones.

---

## Step 8 — Plot Results

### Dimensional Loads

Extract a sorted 1D polar slice with [`VLMLoadSlice`](@ref):

```julia
using Plots

polar = VLMLoadSlice(setup.loads; beta=0.0)

Plots.plot(polar.alpha, polar.L_total,
    xlabel = "Angle of Attack α [deg]",
    ylabel = "Lift [N]",
    label  = "Total Lift",
    lw     = 2,
)
```

### Non-Dimensional Coefficients

For a single solved flight condition, use [`VLMCoefficients`](@ref) with the
setup and its `(alpha, beta)` pair:

```julia
coefs = VLMCoefficients(setup, 4.0, 0.0)

# Whole-aircraft coefficients
println("CL = ", coefs.CL_total)
println("CD = ", coefs.CD_total)
println("CM = ", coefs.CM_total)

# Per-surface coefficients
println("Lift coefficients by surface = ", coefs.CL)
println("Drag coefficients by surface = ", coefs.CD)
```

!!! danger "Only for calculated data points!"
    The operating point must already have been solved by [`VLMSolver!`](@ref).
    The `*_total` fields contain whole-aircraft values, while fields such as `CL`
    and `CD` contain one value for each surface.

For coefficient curves across an angle-of-attack or sideslip sweep, use
[`VLMCoefficientsSlice`](@ref):

```julia
coefs = VLMCoefficientsSlice(setup; beta=0.0)

# Lift curve
Plots.plot(coefs.alpha, coefs.CL_total,
    xlabel = "α [deg]", ylabel = "CL", label = "CL", lw = 2,
)

# Drag polar
Plots.plot(coefs.CD_total, coefs.CL_total,
    xlabel = "CD", ylabel = "CL", label = "Near-field", lw = 2,
)
Plots.plot!(coefs.CD_trefftz_total, coefs.CL_total,
    label = "Trefftz Induced Drag", ls = :dash,
)
```

!!! tip "Trefftz drag vs near-field drag"
    The Trefftz plane induced drag (`CD_trefftz`) is significantly more accurate and
    mesh-independent than the near-field surface-integration drag (`CD`). Prefer it
    for reporting induced drag values. See [VLM Theory](../technical/vlm_theory.md) for details.

---

## Next Steps

| Topic | Page |
|:------|:-----|
| Coefficient breakdowns and per-surface analysis | [Post-Processing](post_processing.md) |
| Longitudinal and lateral stability derivatives | [Stability Analysis](stability_analysis.md) |
| Ground proximity effects | [Ground Effect](ground_effect.md) |
| Exporting geometry for CAD/SolidWorks | [CAD Export](cad_export.md) |
| Mathematical and numerical foundations | [VLM Theory](../technical/vlm_theory.md) |
| Solver pipeline and caching internals | [Solver Architecture](../technical/solver_architecture.md) |
