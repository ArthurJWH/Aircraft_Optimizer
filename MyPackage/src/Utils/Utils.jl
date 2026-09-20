module Utils

include("structBC.jl")
include("structSpline.jl")
include("structLSR.jl")
include("structInterpolate.jl")
include("structIntegrate.jl")
include("bisection.jl")

export AbstractBC,
    FirstDerivativeBC, SecondDerivativeBC, ThirdDerivativeBC, NaturalBC
export AbstractSpline,
    LinearSpline,
    QuadraticSpline,
    CubicSpline,
    LinearSplineFromCoeffs,
    QuadraticSplineFromCoeffs,
    CubicSplineFromCoeffs
export LSR
export Interpolate
export IntegrateGLQ
export bisection

end
