"""
    Data

Mutable container holding general geometric and flight performance properties of a plane.

# Fields

  - `CG::Tuple{Float64, Float64, Float64}`: Center of gravity coordinates `(x, y, z)`, in `m`, in aircraft body axes.
  - `SM::Float64`: Static margin.
  - `MTOW::Float64`: Maximum takeoff weight, in `kg`.
  - `alpha_stall::Float64`: Stall angle of attack, in `deg`.
  - `beta_stall::Float64`: Stall sideslip angle, in `deg`.
  - `alpha_trim::Float64`: Trim angle of attack, in `deg`.
"""
mutable struct Data
    CG::Tuple{Float64, Float64, Float64}
    SM::Float64
    MTOW::Float64
    alpha_stall::Float64
    beta_stall::Float64
    alpha_trim::Float64
end

function Data(;
    CG::Tuple{Float64, Float64, Float64}=(0.0, 0.0, 0.0),
    SM::Float64=0.0,
    MTOW::Float64=0.0,
    alpha_stall::Float64=0.0,
    beta_stall::Float64=0.0,
    alpha_trim::Float64=0.0,
)
    return Data(CG, SM, MTOW, alpha_stall, beta_stall, alpha_trim)
end
