# Deliberately untyped: the shared core both other methods funnel into,
# over anything iterable of (Float64, Float64) pairs -- a Vector, a lazy
# alpha x beta Generator, or a 1-tuple from the single-point method below.
"""
    VLMSolver!(setup::VLMSetup, alpha::AbstractVector{<:Float64}, beta::AbstractVector{<:Float64})
    VLMSolver!(setup::VLMSetup, points::AbstractVector{NTuple{2, Float64}})
    VLMSolver!(setup::VLMSetup, rot::NTuple{2, Float64})

Executes the Vortex Lattice Method solver for specified angle of attack (`alpha`) and sideslip (`beta`) conditions (in `deg`).

All methods mutate `setup.loads` in place and return `setup.loads`.

# Signatures

  - `alpha, beta`: Evaluates every condition in the Cartesian product of `alpha x beta`.
  - `points`: Evaluates an explicit collection of `(alpha, beta)` coordinate pairs.
  - `rot = (alpha, beta)`: Convenience wrapper evaluating a single operating condition.

# Caching & Thread-Safety

  - **Deduplication**: If an `(alpha, beta)` point already exists in `setup.loads.points`, it is skipped.
  - **Thread Safety**: Access to `setup.loads` and the ring geometry cache is synchronized via `setup.lock`.
  - **Internal Parallelism**: For each condition solved, panel AIC matrix assembly and Trefftz plane induced drag integration are multi-threaded via `Base.Threads.@threads`.

# Returns

  - `setup.loads::VLMLoad`: The updated loads repository containing all solved operating points.

# Example

```julia
# Solve an alpha sweep from -4 to 12 degrees at beta = 0.0
VLMSolver!(setup, collect(-4.0:2.0:12.0), [0.0])

# Access results
polar = VLMLoadSlice(setup.loads; beta=0.0)
```
"""
function VLMSolver!(setup::VLMSetup, points)::VLMLoad
    lock(setup.lock) do
        for (a, b) in points
            haskey(setup.loads.points, (a, b)) && continue
            setup.loads.points[(a, b)] = _solve_point(setup, a, b)
        end
    end
    return setup.loads
end

function VLMSolver!(
    setup::VLMSetup,
    alpha::AbstractVector{<:Float64},
    beta::AbstractVector{<:Float64},
)::VLMLoad
    return VLMSolver!(setup, ((a, b) for b in beta, a in alpha))
end

function VLMSolver!(setup::VLMSetup, rot::NTuple{2, Float64})::VLMLoad
    return VLMSolver!(setup, (rot,))
end

# Solve a single flow condition, returning its coefficients (does not touch
# setup.polar -- that's VLMSolver!'s job).
function _solve_point(
    setup::VLMSetup, alpha::Float64, beta::Float64
)::VLMLoadPoint
    geometry = setup.geometry
    rings, n_panels, n_surfaces, surfaces, wake_map = geometry.rings,
    geometry.n_panels, geometry.n_surfaces, geometry.surfaces,
    geometry.wake_map
    CG = geometry.CG
    V_inf, rho, epsilon2 = setup.V_inf, setup.rho, setup.epsilon2
    rot = (alpha, beta)

    if setup.ground
        AIC, AIC_vx, AIC_vy, AIC_vz, RHS, V_dir = _assemble_sys(
            rings, n_panels, rot, wake_map, V_inf, setup.h, CG, epsilon2
        )
    else
        ring_geom = _ring_geom!(setup)
        AIC, AIC_vx, AIC_vy, AIC_vz, RHS, V_dir = _assemble_sys(
            rings, n_panels, rot, wake_map, V_inf, ring_geom, epsilon2
        )
    end

    gamma = lu(AIC) \ RHS

    forces, span_segments = _calc_forces(
        gamma,
        n_panels,
        n_surfaces,
        surfaces,
        rings,
        rho,
        V_dir,
        V_inf,
        AIC_vx,
        AIC_vy,
        AIC_vz,
    )

    (FX, FX_dist, FY, FY_dist, FZ, FZ_dist, L, L_dist, D, D_dist, M, M_dist, Ml, Ml_dist, N, N_dist) = _calc_loads(
        forces, span_segments, rings, CG, n_surfaces, rot
    )
    D_trefftz, D_trefftz_dist = _calc_trefftz_drag(
        gamma, rings, surfaces, rho, V_dir, epsilon2
    )

    return VLMLoadPoint(
        FX,
        FX_dist,
        sum(FX),
        FY,
        FY_dist,
        sum(FY),
        FZ,
        FZ_dist,
        sum(FZ),
        L,
        L_dist,
        sum(L),
        D,
        D_dist,
        sum(D),
        M,
        M_dist,
        sum(M),
        Ml,
        Ml_dist,
        sum(Ml),
        N,
        N_dist,
        sum(N),
        D_trefftz,
        D_trefftz_dist,
        sum(D_trefftz),
    )
end

# Populates (once) and returns setup's cached angle-independent ring-induced
# velocity components. Only valid to call while setup.lock is already held
# by the current task -- both call sites (_solve_point, via VLMSolver!) hold
# it for the whole operation, and ReentrantLock permits the same task to
# re-enter, so this doesn't re-lock itself.
function _ring_geom!(setup::VLMSetup)::NTuple{3, Matrix{Float64}}
    cached = setup.ring_geom[]
    cached !== nothing && return cached
    geometry = setup.geometry
    result = _assemble_ring_geom(
        geometry.rings, geometry.n_panels, setup.epsilon2
    )
    setup.ring_geom[] = result
    return result
end
