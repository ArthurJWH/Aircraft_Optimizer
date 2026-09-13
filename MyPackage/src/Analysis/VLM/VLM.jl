module VLM

using LinearAlgebra
using StaticArrays
using Base.Threads

using ..Utils
using ..Geometry

include("structVLMMesh.jl")
include("VLMStructs.jl")
include("structVLMLoad.jl")
include("structVLMGeometry.jl")
include("vlm_assemble_aic.jl")
include("vlm_calc_forces.jl")
include("structVLMSetup.jl")
include("vlm_solver.jl")
include("structVLMCoefs.jl")
include("structVLMStab.jl")

export VLMMesh, VLMLoad, VLMLoadPoint, VLMLoadSlice, VLMGeometry, VLMSetup,
    VLMSolver!, VLMCoefficients, VLMCoefficientsSlice, VLMStabilityDerivatives,
    VLMStabilityDerivativesFast, VortexRing, VLMSurface, GroundTransform, Vec3

end
