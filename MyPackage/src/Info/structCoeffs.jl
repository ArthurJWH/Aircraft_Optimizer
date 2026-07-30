mutable struct Coeffs
    CX::Float64 # Force coefficient in x direction
    CY::Float64 # Force coefficient in y direction
    CZ::Float64 # Force coefficient in z direction
    CL::Float64 # Lift coefficient
    CD::Float64 # Drag coefficient
    CM::Float64 # Pitching moment coefficient
    CMl::Float64 # Roll moment coefficient
    CN::Float64 # Yaw moment coefficient

    CX_surf::Float64
    CY_surf::Float64
    CZ_surf::Float64
    CL_surf::Float64
    CD_surf::Float64
    CM_surf::Float64
    CMl_surf::Float64
    CN_surf::Float64
end

function Coeffs()
    return Coeffs(ntuple(_ -> 0.0, 16)...)
end
