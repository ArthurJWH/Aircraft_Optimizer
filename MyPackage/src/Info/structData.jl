"""
    Data

    A mutable struct to hold general aerodynamic and geometric data of a plane.

    Fields
    ------
    CG : Tuple{Float64, Float64, Float64}
        The center of gravity of the plane.
    SM : Float64
        The static margin of the plane.
    MTOW : Float64
        The maximum takeoff weight of the plane.
    alpha_stall : Float64
        The angle of attack at which the plane stalls.
    beta_stall : Float64
        The sideslip angle at which the plane stalls.
    alpha_trim : Float64
        The angle of attack at which the plane is trimmed.
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
