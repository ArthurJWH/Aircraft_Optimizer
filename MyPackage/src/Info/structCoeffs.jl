"""
    Coeffs

Mutable container storing aircraft and per-surface aerodynamic coefficients.

# Fields
- `CX::Float64`: Total body-axis force coefficient in X (aft positive).
- `CY::Float64`: Total body-axis force coefficient in Y (starboard positive).
- `CZ::Float64`: Total body-axis force coefficient in Z (upward positive).
- `CL::Float64`: Total wind-axis lift coefficient.
- `CD::Float64`: Total wind-axis drag coefficient.
- `CM::Float64`: Total pitching moment coefficient.
- `CMl::Float64`: Total rolling moment coefficient.
- `CN::Float64`: Total yawing moment coefficient.
- `CX_surf::Vector{Float64}`: Per-surface X force coefficient.
- `CY_surf::Vector{Float64}`: Per-surface Y side-force coefficient.
- `CZ_surf::Vector{Float64}`: Per-surface Z force coefficient.
- `CL_surf::Vector{Float64}`: Per-surface lift coefficient.
- `CD_surf::Vector{Float64}`: Per-surface drag coefficient.
- `CM_surf::Vector{Float64}`: Per-surface pitching moment coefficient.
- `CMl_surf::Vector{Float64}`: Per-surface rolling moment coefficient.
- `CN_surf::Vector{Float64}`: Per-surface yawing moment coefficient.
"""
mutable struct Coeffs
    CX::Float64 # Force coefficient in x direction
    CY::Float64 # Force coefficient in y direction
    CZ::Float64 # Force coefficient in z direction
    CL::Float64 # Lift coefficient
    CD::Float64 # Drag coefficient
    CM::Float64 # Pitching moment coefficient
    CMl::Float64 # Roll moment coefficient
    CN::Float64 # Yaw moment coefficient

    CX_surf::Vector{Float64}
    CY_surf::Vector{Float64}
    CZ_surf::Vector{Float64}
    CL_surf::Vector{Float64}
    CD_surf::Vector{Float64}
    CM_surf::Vector{Float64}
    CMl_surf::Vector{Float64}
    CN_surf::Vector{Float64}
end

function Coeffs()
    return Coeffs(
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
    )
end
