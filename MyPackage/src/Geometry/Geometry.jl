module Geometry

include("structAirfoil.jl")
include("structAerosurface.jl")
include("structPlane.jl")
include("AbstractMeshType.jl")

export Airfoil, calc_surfaces, calc_camber
export Aerosurface
export Plane
export AbstractMesh

end
