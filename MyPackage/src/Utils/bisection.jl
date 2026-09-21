"""
    bisection(f, a, b; tol=1e-8, max_iter=100)

Finds a root of a continuous scalar function `f(x) = 0` on the interval `[a, b]` using the bisection method.

Requires `f(a)` and `f(b)` to have opposite signs (`f(a) * f(b) <= 0`).

# Arguments

  - `f::Function`: Univariate continuous function.
  - `a::Real`: Lower bound of the initial search interval.
  - `b::Real`: Upper bound of the initial search interval.
  - `tol::Float64`: Relative tolerance stopping criterion `|b - a| <= tol * |c|` (default: `1e-8`).
  - `max_iter::Int`: Maximum allowable bisection iterations (default: `100`).

# Returns

  - `c::Float64`: Approximated root such that `f(c) ≈ 0`.

# Example

```julia
using MyPackage.Utils
root = bisection(x -> x^2 - 4.0, 0.0, 5.0)  # Returns ≈ 2.0
```
"""
function bisection(f, a, b; ftol=1e-8, tol=1e-8, max_iter=100)
    fa = f(a)
    fb = f(b)

    if abs(fa) < ftol
        return a
    elseif abs(fb) < ftol
        return b
    elseif fa * fb > 0
        error("f(a) and f(b) must have opposite signs")
    end

    for i in 1:max_iter
        c = (a + b) / 2
        fc = f(c)

        if abs(b - a) <= tol * abs(c)
            return c
        end

        if fa * fc < 0
            b = c
            fb = fc
        else
            a = c
            fa = fc
        end
    end

    error("Maximum iterations reached without convergence")
end
