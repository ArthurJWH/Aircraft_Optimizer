# VLMSetup: one fixed flow configuration against a VLMGeometry, plus the
# per-condition solve (_solve_point) and the public VLMSolver! entry
# points. Depends on VLMGeometry/VLMPolar (structs, must be defined
# first) and calls into vlm_aic.jl/vlm_forces.jl (functions, resolved at
# call time -- no ordering constraint from that).

"""
    VLMSetup

Encapsulates a fixed flow physics configuration and solver state for a [`VLMGeometry`](@ref).

A `VLMSetup` couples a discretized aircraft geometry to specific atmospheric and boundary conditions
(freestream velocity, fluid density, ground effect settings, and numerical vortex core regularization).
It maintains a thread-safe storage container [`VLMLoad`](@ref) for solved operating points.

# Fields

  - `geometry::VLMGeometry`: Reference to the discretized mesh and vortex ring geometry.
  - `V_inf::Float64`: Freestream airspeed, in `m/s`.
  - `rho::Float64`: Atmospheric density, in `kg/m³` (default: `1.225`).
  - `ground::Bool`: Whether ground effect modeling is enabled (default: `false`).
  - `h::Float64`: Height of the aircraft center of gravity `CG` above the ground plane, in `m` (default: `0.0`).
  - `epsilon2::Float64`: Dimensionless squared core-radius fraction for vortex singularity regularization (default: `1e-10`).
  - `loads::VLMLoad`: Thread-safe dictionary container storing solved [`VLMLoadPoint`](@ref) results.
  - `ring_geom::Base.RefValue{Union{Nothing, NTuple{3, Matrix{Float64}}}}`: Lazy cache of angle-independent ring-induced velocities.
  - `lock::ReentrantLock`: Concurrency lock ensuring safe multi-threaded insertions into `loads` and cache population.

# Constructors

```julia
VLMSetup(
    geometry::VLMGeometry,
    V_inf::Float64;
    rho::Float64=1.225,
    ground::Bool=false,
    h::Float64=0.0,
    epsilon2::Float64=1e-10,
)
```

# Numerical Caching

For free-flight simulations (`ground == false`), the Biot-Savart velocity influence induced by the closed vortex rings
on each panel collocation point depends exclusively on geometry and `epsilon2`. This angle-independent velocity field
(three ``N \\times N`` matrices for ``V_x, V_y, V_z``) is computed lazily on the first solve and cached in `ring_geom`,
saving significant ``O(N^2)`` work across subsequent angle-of-attack sweeps.
"""
struct VLMSetup
    geometry::VLMGeometry
    V_inf::Float64
    rho::Float64
    ground::Bool
    h::Float64
    epsilon2::Float64
    loads::VLMLoad
    ring_geom::Base.RefValue{Union{Nothing, NTuple{3, Matrix{Float64}}}}
    lock::ReentrantLock
end

function VLMSetup(
    geometry::VLMGeometry,
    V_inf::Float64;
    rho::Float64=1.225,
    ground::Bool=false,
    h::Float64=0.0,
    epsilon2::Float64=1e-10,
)
    return VLMSetup(
        geometry,
        V_inf,
        rho,
        ground,
        h,
        epsilon2,
        VLMLoad(),
        Ref{Union{Nothing, NTuple{3, Matrix{Float64}}}}(nothing),
        ReentrantLock(),
    )
end
