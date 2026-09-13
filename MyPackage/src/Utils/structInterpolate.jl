"""
    Interpolate(xs, fs)
    Interpolate(coeffs)

Polynomial interpolation passing exactly through all `n` given data points `(xs, fs)` (degree ``n - 1``).

# Fields
- `xs::Vector{Float64}`: Monotonically increasing knot coordinates.
- `fs::Vector{Float64}`: Function values at knot coordinates.
- `coeffs::Vector{Float64}`: Calculated polynomial coefficients.

# Example
```julia
xs = [0.0, 1.0, 2.0]
fs = [1.0, 2.0, 0.0]
interp = Interpolate(xs, fs)
val = interp(1.5)  # Evaluates polynomial at x = 1.5
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
