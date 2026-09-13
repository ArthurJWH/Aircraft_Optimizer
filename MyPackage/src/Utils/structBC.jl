"""
    AbstractBC

Abstract supertype for boundary conditions applied to spline interpolation in `MyPackage.Utils`.
"""
abstract type AbstractBC end

"""
    NopBC

Null boundary condition indicating no explicit derivative constraint applied.

# Example

```julia
bc = NopBC()  # No boundary condition
```
"""
struct NopBC <: AbstractBC end

"""
    FirstDerivativeBC(value, index)

Boundary condition constraining the first derivative ``f'(x)`` to `value` at knot index `index` (or `:left` / `:right`).

# Fields

  - `value::Float64`: Target first derivative value.
  - `index::T`: Index (or symbol `:left` / `:right`) where the boundary condition is enforced.

# Example

```julia
bc = FirstDerivativeBC(0.0, :left)  # Zero first derivative at left boundary
```
"""
struct FirstDerivativeBC{T} <: AbstractBC
    value::Float64
    index::T
end

"""
    SecondDerivativeBC(value, index)

Boundary condition constraining the second derivative ``f''(x)`` to `value` at knot index `index` (or `:left` / `:right`).

# Fields

  - `value::Float64`: Target second derivative value.
  - `index::T`: Index (or symbol `:left` / `:right`) where the boundary condition is enforced.

# Example

```julia
bc = SecondDerivativeBC(0.0, :left)  # Natural / zero-curvature boundary condition at left boundary
```
"""
struct SecondDerivativeBC{T} <: AbstractBC
    value::Float64
    index::T
end

"""
    ThirdDerivativeBC(value, index)

Boundary condition constraining the third derivative ``f'''(x)`` to `value` at knot index `index` (or `:left` / `:right`).

# Fields

  - `value::Float64`: Target third derivative value.
  - `index::T`: Index (or symbol `:left` / `:right`) where the boundary condition is enforced.

# Example

```julia
bc = ThirdDerivativeBC(0.0, :right)  # Third derivative constraint at right boundary
```
"""
struct ThirdDerivativeBC{T} <: AbstractBC
    value::Float64
    index::T
end
