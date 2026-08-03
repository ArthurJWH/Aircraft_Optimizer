# TODO: Implement a integration method that takes a list of x values and corresponding f(x) values. This method should be able to handle uneven intervals.
# function integrate(f::Function; a::Float64 = 0.0, b::Float64 = 1.0, n::Int = 100)
#     h = (b - a) / n
#     x = a:h:b
#     fs = f.(x)

#     # Integrate the cubic spline using Simpson's rule
#     integral = f(x[1]) + 4*f(x[2])
#     for i in 3:2:n-3
#         integral += 2*f(x[i]) + 4*f(x[i+1])
#     end
#     integral *= h / 3

#     if n % 2 == 0
#         integral += (2*f(x[n-1]) + 4*f(x[n]) + f(x[n+1])) * h / 3
#     else
#         integral += (3/8 + 1/3) * (f(x[n-2])) * h + (3/8) * (3 * f(x[n-1]) + 3 * f(x[n]) + f(x[n+1])) * h
#     end

#     return integral
# end

"""
    legendre

    A constant that stores the roots and weights for Gauss-Legendre quadrature for n=1 to 5.
    Each entry is a tuple containing two arrays:
    - The first array contains the roots of the Legendre polynomial of degree n.
    - The second array contains the corresponding weights for the quadrature.
"""
const legendre = (
    nothing,
    ((0.0,), (2.0,)),
    ((-0.5773502691896257, 0.5773502691896257), (1.0, 1.0)),
    (
        (-0.7745966692414834, 0.0, 0.7745966692414834),
        (0.5555555555555556, 0.8888888888888888, 0.5555555555555556),
    ),
    (
        (
            -0.8611363115940526,
            -0.3399810435848563,
            0.3399810435848563,
            0.8611363115940526,
        ),
        (
            0.34785484513745385,
            0.6521451548625461,
            0.6521451548625461,
            0.34785484513745385,
        ),
    ),
    (
        (
            -0.906179845938664,
            -0.5384693101056831,
            0.0,
            0.5384693101056831,
            0.906179845938664,
        ),
        (
            0.23692688505618908,
            0.47862867049936647,
            0.5688888888888889,
            0.47862867049936647,
            0.23692688505618908,
        ),
    ),
)

"""
    IntegrateGLQ

    A struct for performing numerical integration using Gauss-Legendre quadrature.

    Fields
    ------
    f : F
        The function to be integrated.
    n : Int
        The number of points to use in the quadrature (1 to 5).

    Example
    -------
    ```julia
    f(x) = x^2
    integrator = IntegrateGLQ(f, 3)
    result = integrator(0.0, 1.0)  # Integrates f(x) from 0 to 1
    ```
"""
struct IntegrateGLQ{F}
    f::F
    n::Int
end

function IntegrateGLQ(f; n::Int=3)
    @assert n >= 1 && n <= 5 "n must be between 1 and 5"
    return IntegrateGLQ(f, n)
end

function (integrateglq::IntegrateGLQ)(a::AbstractFloat, b::AbstractFloat)
    return evaluate(integrateglq, a, b)
end

@inline function evaluate(
    integrateglq::IntegrateGLQ, a::AbstractFloat, b::AbstractFloat
)
    # Get the roots and weights for the Legendre polynomial of degree n
    x, w = legendre[integrateglq.n]
    f = integrateglq.f

    # Transform coefficients roots->[a, b]
    t1 = 0.5 * (b - a)
    t2 = 0.5 * (b + a)

    # Compute the integral using the weights and transformed roots
    integral = 0.0
    @inbounds for i in eachindex(x)
        integral += w[i] * f(t1 * x[i] + t2)
    end

    return integral * 0.5 * (b - a)
end

#TODO: Implement Richardson extrapolation for numerical integration to improve accuracy.
