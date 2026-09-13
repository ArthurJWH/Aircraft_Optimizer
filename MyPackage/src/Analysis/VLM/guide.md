# VLM Solver — quick guide

Three steps: mesh, configure, solve. Each is a separate object so you can
build/inspect the mesh before paying for a solve, and reuse it across as
many flow conditions as you want.

## 1. Build the geometry (mesh only, no solving)

```julia
geom = VLMGeometry(plane, [(10, 5)])   # n_chordxspan per surface
```

`geom.meshes` is the per-surface mesh — inspect or plot it now, before
anything gets solved. `VLMGeometry` never runs the VLM; it just builds the
mesh and the derived vortex-ring geometry once, so it can be shared cheaply
across every condition below.

## 2. Fix a flow configuration

```julia
setup = VLMSetup(geom, 50.0; rho=1.225, ground=false, h=0.0, epsilon2=1e-10)
```

`V_inf` is positional; everything else has defaults. A `VLMSetup` pins down
everything that affects the *values* a solve produces — `V_inf`, `rho`,
`ground`, `h`, `epsilon2` — so its cache is always self-consistent. Want a
different altitude or a ground-effect case? Build another (cheap)
`VLMSetup` from the same `geom` rather than reusing one for both:

```julia
setup_ge = VLMSetup(geom, 50.0; ground=true, h=2.0)
```

## 3. Solve

```julia
VLMSolver!(setup, [0.0, 2.0, 4.0, 6.0], [0.0])   # sweep
VLMSolver!(setup, (8.0, 0.0))                     # single point
```

Both mutate `setup.polar` in place — that's the `!`. No return value to
reassign; `setup.polar` already has everything. Calling again with values
already solved (exact match) is a cheap no-op, not a recompute.

## 4. Read results

```julia
pt = setup.polar.points[(2.0, 0.0)]   # a VLMPolarPoint
pt.L                                   # per-surface lift, Vector{Float64}
pt.L_dist                              # per-surface spanwise lift, Vector{Vector{Float64}}
pt.L_total                             # whole-aircraft lift, Float64 (sum over surfaces)

slice = VLMLoadSlice(setup.polar; beta=0.0)   # lift-vs-alpha polar at beta=0
slice.alpha                            # sorted alpha values actually solved
slice.L[1, :]                          # surface 1's lift across that slice
slice.L_total                          # whole-aircraft lift across that slice

slice2 = VLMLoadSlice(setup.polar; alpha=4.0) # same, fixing alpha instead
```

`VLMLoadSlice` only exposes the per-surface and total levels (that's what
plots as a matrix/vector cleanly) — for spanwise distributions, index into
`setup.polar.points` directly per point.

`VLMLoadSlice` only returns points you actually solved for — it errors if
nothing matches rather than padding gaps, so if a plot looks sparse it's
telling you to call `VLMSolver!` for more points, not silently faking data.

## End to end

```julia
geom  = VLMGeometry(plane, [(10, 5)])
setup = VLMSetup(geom, 50.0)
VLMSolver!(setup, collect(-4.0:2.0:10.0), [0.0])
polar = VLMLoadSlice(setup.polar; beta=0.0)
# polar.alpha, polar.L, polar.D, ... ready to plot
```

## Gotchas

- Alpha/beta matching is **exact float equality** — reuse the same literals
  if you want a later call to hit the cache instead of resolving.
- A `VLMSetup`'s cache can be shared/mutated from multiple call sites
  (it's locked internally), but two `VLMSetup`s never share results even if
  built from the same `geom` — each has its own `polar`.
- No per-condition AIC matrix or rings are kept — this is a polar-(plus-
  spanwise-per-point) API, not a full-setup API. If you need the raw AIC or
  rings for one specific condition, `_solve_point(setup, alpha, beta)` is
  the internal function that computes them before discarding everything but
  the coefficients.
- For `ground=false`, the first solve under a `VLMSetup` computes and caches
  the angle-independent part of the AIC assembly (~3x one AIC matrix in
  memory), reused by every solve after that under the same `setup`. This
  does *not* happen for `ground=true` — the ground-image geometry rotates
  with alpha, so there's no angle-independent part to cache without a
  second, alpha-keyed cache layer that isn't built yet.