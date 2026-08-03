raw"""
    Interpolate

    Create an interpolating function from given x-coordinates and function values.
    The interpolant is a polynomial of degree n-1, where n is the number of data points.

    Fields
    ------
    xs : Vector{Float64}
        The x-coordinates of the data points.
    fs : Vector{Float64}
        The function values at the data points.
    coeffs : Vector{Float64}
        The calculated coefficients of the interpolating polynomial.

    Example
    -------
    ```julia
    xs = [0.0, 1.0, 2.0]
    fs = [1.0, 2.0, 0.0]
    interp = Interpolate(xs, fs)
    f = interp(1.5)  # Evaluate the interpolating polynomial at x=1.5
    ```
"""
mutable struct Interpolate
    xs::Vector{Float64}
    fs::Vector{Float64}
    coeffs::Vector{Float64}
end

function Interpolate(coeffs::AbstractVector{<:AbstractFloat})
    return Interpolate([0.0], [0.0], coeffs)
end

function Interpolate(
    xs::AbstractVector{<:AbstractFloat}, fs::AbstractVector{<:AbstractFloat}
)
    @assert length(xs) == length(fs) "xs and fs must have same length"
    @assert issorted(xs) "xs must be sorted"
    n = length(xs)
    return LSR(xs, fs, n - 1)
end

(interp::Interpolate)(x::AbstractFloat) = evaluate(interp, x)

@inline function evaluate(interp::Interpolate, x::AbstractFloat)
    c = interp.coeffs
    result = 0.0
    for i in length(c):-1:1
        result = result * x + c[i]
    end
    return result
end
