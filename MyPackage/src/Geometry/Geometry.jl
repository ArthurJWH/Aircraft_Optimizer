module Geometry

include("structAirfoil.jl")
include("structAerosurface.jl")
include("structPlane.jl")

export Airfoil, calc_surfaces, calc_camber
export Aerosurface
export Plane

end
