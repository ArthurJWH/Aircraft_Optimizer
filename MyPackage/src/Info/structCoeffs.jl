"""
    Coeffs

    A mutable struct to store aerodynamic coefficients for a plane.

    Fields
    ------
    CX : Float64
        Force coefficient in x direction.
    CY : Float64
        Force coefficient in y direction.
    CZ : Float64
        Force coefficient in z direction.
    CL : Float64
        Lift coefficient.
    CD : Float64
        Drag coefficient.
    CM : Float64
        Pitching moment coefficient.
    CMl : Float64
        Roll moment coefficient.
    CN : Float64
        Yaw moment coefficient.
    CX_surf : Vector{Float64}
        Force coefficient in x direction per surface.
    CY_surf : Vector{Float64}
        Force coefficient in y direction per surface.
    CZ_surf : Vector{Float64}
        Force coefficient in z direction per surface.
    CL_surf : Vector{Float64}
        Lift coefficient per surface.
    CD_surf : Vector{Float64}
        Drag coefficient per surface.
    CM_surf : Vector{Float64}
        Pitching moment coefficient per surface.
    CMl_surf : Vector{Float64}
        Roll moment coefficient per surface.
    CN_surf::Vector{Float64}
        Yaw moment coefficient per surface.
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
    return Coeffs(ntuple(_ -> 0.0, 16)...)
end
