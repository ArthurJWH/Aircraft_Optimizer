using SparseArrays

raw"""
    AbstractSpline

    Abstract type for spline interpolation.
    Used in function declarations to support different spline types.
"""
abstract type AbstractSpline end

raw"""
    LinearSpline

    Struct for linear spline interpolation.
    A relative-coordinate spline is applied, where the spline evaluated at a point `x` in the interval `[x_i, x_{i+1}]` is given by:

    ```math
        spl_i(x) = coeffs[2\cdot i - 1] + coeffs[2\cdot i] (x - x_i)
    ```

    Fields
    ------
    xs : Vector{Float64}
        The x-coordinates of the data points.
    fs : Vector{Float64}
        The function values at the data points.
    coeffs : Vector{Float64}
        The computed coefficients of the linear spline segments.

    Example
    -------
    ```julia
    xs = [0.0, 1.0, 2.0]
    fs = [1.0, 2.0, 0.0]
    spl = LinearSpline(xs, fs)
    f = spl(1.5) # Evaluate the spline at x = 1.5
    ```
"""
mutable struct LinearSpline <: AbstractSpline
    xs::Vector{Float64}
    fs::Vector{Float64}
    coeffs::Vector{Float64}
end

raw"""
    QuadraticSpline

    Struct for quadratic spline interpolation.
    A relative-coordinate spline is applied, where the spline evaluated at a point `x` in the interval `[x_i, x_{i+1}]` is given by:

    ```math
        spl_i(x) = coeffs[3\cdot i - 2] + coeffs[3\cdot i - 1] (x - x_i) + coeffs[3\cdot i] (x - x_i)^2
    ```

    Fields
    ------
    xs : Vector{Float64}
        The x-coordinates of the data points.
    fs : Vector{Float64}
        The function values at the data points.
    coeffs : Vector{Float64}
        The computed coefficients of the quadratic spline segments.
    bc : AbstractBC
        The boundary condition applied to the spline.

    Example
    -------
    ```julia
    xs = [0.0, 1.0, 2.0]
    fs = [1.0, 2.0, 0.0]
    bc = SecondDerivativeBC(0.0, :left)
    spl = QuadraticSpline(xs, fs; bc=bc)
    f = spl(1.5) # Evaluate the spline at x = 1.5
    ```
"""
mutable struct QuadraticSpline <: AbstractSpline
    xs::Vector{Float64}
    fs::Vector{Float64}
    coeffs::Vector{Float64}
    bc::AbstractBC
end

raw"""
    CubicSpline

    Struct for cubic spline interpolation.
    A relative-coordinate spline is applied, where the spline evaluated at a point `x` in the interval `[x_i, x_{i+1}]` is given by:

    ```math
        spl_i(x) = coeffs[4\cdot i - 3] + coeffs[4\cdot i - 2] (x - x_i) + coeffs[4\cdot i - 1] (x - x_i)^2 + coeffs[4\cdot i] (x - x_i)^3
    ```

    Fields
    ------
    xs : Vector{Float64}
        The x-coordinates of the data points.
    fs : Vector{Float64}
        The function values at the data points.
    coeffs : Vector{Float64}
        The computed coefficients of the cubic spline segments.
    bc1 : AbstractBC
        The first boundary condition applied to the spline.
    bc2 : AbstractBC
        The second boundary condition applied to the spline.

    Example
    -------
    ```julia
    xs = [0.0, 1.0, 2.0]
    fs = [1.0, 2.0, 0.0]
    bc1 = SecondDerivativeBC(0.0, :left)
    bc2 = SecondDerivativeBC(0.0, :right)
    spl = CubicSpline(xs, fs; bc1=bc1, bc2=bc2)
    f = spl(1.5) # Evaluate the spline at x = 1.5
    ```
"""
mutable struct CubicSpline <: AbstractSpline
    xs::Vector{Float64}
    fs::Vector{Float64}
    coeffs::Vector{Float64}
    bc1::AbstractBC
    bc2::AbstractBC
end

function LinearSpline(
    xs::AbstractVector{<:AbstractFloat}, fs::AbstractVector{<:AbstractFloat}
)
    @assert length(xs) == length(fs) "xs and fs must have same length"
    @assert issorted(xs) "xs must be sorted"
    n = length(xs)
    spl = LinearSpline(xs, fs, zeros(2 * (n - 1)))
    calc_spline!(spl)
    return spl
end

function QuadraticSpline(
    xs::AbstractVector{<:AbstractFloat},
    fs::AbstractVector{<:AbstractFloat};
    bc::AbstractBC=SecondDerivativeBC(0.0, :left),
)
    @assert length(xs) == length(fs) "xs and fs must have same length"
    @assert issorted(xs) "xs must be sorted"
    n = length(xs)
    spl = QuadraticSpline(xs, fs, zeros(3 * (n - 1)), bc)
    calc_spline!(spl)
    return spl
end

function CubicSpline(
    xs::AbstractVector{<:AbstractFloat},
    fs::AbstractVector{<:AbstractFloat};
    bc1::AbstractBC = SecondDerivativeBC(0.0, :left),
    bc2::AbstractBC = SecondDerivativeBC(0.0, :right),
)
    @assert length(xs) == length(fs) "xs and fs must have same length"
    @assert issorted(xs) "xs must be sorted"
    n = length(xs)
    spl = CubicSpline(xs, fs, zeros(4 * (n - 1)), bc1, bc2)
    calc_spline!(spl)
    return spl
end

@inline function calc_spline!(spl::LinearSpline)
    xs = spl.xs
    fs = spl.fs
    n = length(xs)
    deltax = diff(xs)
    deltaf = diff(fs)
    @assert n >= 2 "need at least two points"
    for i in 1:(n - 1)
        spl.coeffs[2 * i - 1] = fs[i]
        spl.coeffs[2 * i] = deltaf[i] / deltax[i]
    end
    return nothing
end

@inline function calc_spline!(spl::QuadraticSpline)
    xs = spl.xs
    fs = spl.fs
    n = length(xs)
    deltax = diff(xs)
    @assert n >= 2 "need at least two points"

    M = spzeros(Float64, 3 * (n - 1), 3 * (n - 1))
    rhs = zeros(Float64, 3 * (n - 1))
    si, xi = bc_index(spl.bc, n)

    @views begin

        # spl[1](x[1]) = f[1]
        M[1, 1:3] .= [1, 0, 0]
        rhs[1] = fs[1]

        # spl[1](x[2]) = f[2]
        M[2, 1:3] .= [1, deltax[1], deltax[1]^2]
        rhs[2] = fs[2]

        row = 3

        if si == 1
            # Apply BC
            if isa(spl.bc, FirstDerivativeBC)
                M[row, (3 * si - 2):(3 * si)] .= [0, 1, 2 * (xs[xi] - xs[si])]
                rhs[row] = spl.bc.value
                row += 1
            elseif isa(spl.bc, SecondDerivativeBC)
                M[row, (3 * si - 2):(3 * si)] .= [0, 0, 2]
                rhs[row] = spl.bc.value
                row += 1
            end
        end

        for i in 2:(n - 1)
            if si == i
                # Apply BC
                if isa(spl.bc, FirstDerivativeBC)
                    M[row, (3 * si - 2):(3 * si)] .= [
                        0, 1, 2 * (xs[xi] - xs[si])
                    ]
                    rhs[row] = spl.bc.value
                    row += 1
                elseif isa(spl.bc, SecondDerivativeBC)
                    M[row, (3 * si - 2):(3 * si)] .= [0, 0, 2]
                    rhs[row] = spl.bc.value
                    row += 1
                end
            end

            # spl[i]'(x[i]) = spl[i-1]'(x[i])
            M[row, (3 * i - 5):(3 * i)] .= [
                0, -1, -2 * (deltax[i - 1]), 0, 1, 0
            ]
            rhs[row] = 0
            row += 1

            # spl[i](x[i]) = f[i]
            M[row, (3 * i - 2):(3 * i)] .= [1, 0, 0]
            rhs[row] = fs[i]
            row += 1

            # spl[i](x[i+1]) = f[i+1]
            M[row, (3 * i - 2):(3 * i)] .= [1, deltax[i], deltax[i]^2]
            rhs[row] = fs[i + 1]
            row += 1
        end

        spl.coeffs .= M \ rhs
    end
    return nothing
end

@inline function calc_spline!(spl::CubicSpline)
    xs = spl.xs
    fs = spl.fs
    n = length(xs)
    deltax = diff(xs)
    @assert n >= 2 "need at least two points"

    M = spzeros(Float64, 4 * (n - 1), 4 * (n - 1))
    rhs = zeros(Float64, 4 * (n - 1))
    si1, xi1 = bc_index(spl.bc1, n)
    si2, xi2 = bc_index(spl.bc2, n)

    @views begin

        # spl[1](x[1]) = f[1]
        M[1, 1:4] .= [1, 0, 0, 0]
        rhs[1] = fs[1]

        # spl[1](x[2]) = f[2]
        M[2, 1:4] .= [1, deltax[1], deltax[1]^2, deltax[1]^3]
        rhs[2] = fs[2]

        row = 3

        if si1 == 1
            # Apply BC1
            if isa(spl.bc1, FirstDerivativeBC)
                M[row, (4 * si1 - 3):(4 * si1)] .= [
                    0, 1, 2 * (xs[xi1] - xs[si1]), 3 * (xs[xi1] - xs[si1])^2
                ]
                rhs[row] = spl.bc1.value
                row += 1
            elseif isa(spl.bc1, SecondDerivativeBC)
                M[row, (4 * si1 - 3):(4 * si1)] .= [
                    0, 0, 2, 6 * (xs[xi1] - xs[si1])
                ]
                rhs[row] = spl.bc1.value
                row += 1
            elseif isa(spl.bc1, ThirdDerivativeBC)
                M[row, (4 * si1 - 3):(4 * si1)] .= [0, 0, 0, 6]
                rhs[row] = spl.bc1.value
                row += 1
            end
        end

        if si2 == 1
            # Apply BC2
            if isa(spl.bc2, FirstDerivativeBC)
                M[row, (4 * si2 - 3):(4 * si2)] .= [
                    0, 1, 2 * (xs[xi2] - xs[si2]), 3 * (xs[xi2] - xs[si2])^2
                ]
                rhs[row] = spl.bc2.value
                row += 1
            elseif isa(spl.bc2, SecondDerivativeBC)
                M[row, (4 * si2 - 3):(4 * si2)] .= [
                    0, 0, 2, 6 * (xs[xi2] - xs[si2])
                ]
                rhs[row] = spl.bc2.value
                row += 1
            elseif isa(spl.bc2, ThirdDerivativeBC)
                M[row, (4 * si2 - 3):(4 * si2)] .= [0, 0, 0, 6]
                rhs[row] = spl.bc2.value
                row += 1
            end
        end

        for i in 2:(n - 1)
            if si1 == i
                # Apply BC1
                if isa(spl.bc1, FirstDerivativeBC)
                    M[row, (4 * si1 - 3):(4 * si1)] .= [
                        0, 1, 2 * (xs[xi1] - xs[si1]), 3 * (xs[xi1] - xs[si1])^2
                    ]
                    rhs[row] = spl.bc1.value
                    row += 1
                elseif isa(spl.bc1, SecondDerivativeBC)
                    M[row, (4 * si1 - 3):(4 * si1)] .= [
                        0, 0, 2, 6 * (xs[xi1] - xs[si1])
                    ]
                    rhs[row] = spl.bc1.value
                    row += 1
                elseif isa(spl.bc1, ThirdDerivativeBC)
                    M[row, (4 * si1 - 3):(4 * si1)] .= [0, 0, 0, 6]
                    rhs[row] = spl.bc1.value
                    row += 1
                end
            end

            if si2 == i
                # Apply BC2
                if isa(spl.bc2, FirstDerivativeBC)
                    M[row, (4 * si2 - 3):(4 * si2)] .= [
                        0, 1, 2 * (xs[xi2] - xs[si2]), 3 * (xs[xi2] - xs[si2])^2
                    ]
                    rhs[row] = spl.bc2.value
                    row += 1
                elseif isa(spl.bc2, SecondDerivativeBC)
                    M[row, (4 * si2 - 3):(4 * si2)] .= [
                        0, 0, 2, 6 * (xs[xi2] - xs[si2])
                    ]
                    rhs[row] = spl.bc2.value
                    row += 1
                elseif isa(spl.bc2, ThirdDerivativeBC)
                    M[row, (4 * si2 - 3):(4 * si2)] .= [0, 0, 0, 6]
                    rhs[row] = spl.bc2.value
                    row += 1
                end
            end

            # spl[i]'(x[i]) = spl[i-1]'(x[i])
            M[row, (4 * i - 7):(4 * i)] .= [
                0, -1, -2 * (deltax[i - 1]), -3 * (deltax[i - 1])^2, 0, 1, 0, 0
            ]
            rhs[row] = 0
            row += 1

            # spl[i]''(x[i]) = spl[i-1]''(x[i])
            M[row, (4 * i - 7):(4 * i)] .= [
                0, 0, -2, -6 * (deltax[i - 1]), 0, 0, 2, 0
            ]
            rhs[row] = 0
            row += 1

            # spl[i](x[i]) = f[i]
            M[row, (4 * i - 3):(4 * i)] .= [1, 0, 0, 0]
            rhs[row] = fs[i]
            row += 1

            # spl[i](x[i+1]) = f[i+1]
            M[row, (4 * i - 3):(4 * i)] .= [
                1, deltax[i], deltax[i]^2, deltax[i]^3
            ]
            rhs[row] = fs[i + 1]
            row += 1
        end

        spl.coeffs .= M \ rhs
    end
    return nothing
end

(spl::AbstractSpline)(x::AbstractFloat) = evaluate(spl, x)

@inline function evaluate(spl::LinearSpline, x::AbstractFloat)
    i = spl_index(spl.xs, x)
    c = spl.coeffs
    t = x - spl.xs[i]
    # Horner's polynomial evaluation
    @inbounds return c[2 * i - 1] + c[2 * i] * t
end

@inline function evaluate(spl::QuadraticSpline, x::AbstractFloat)
    i = spl_index(spl.xs, x)
    c = spl.coeffs
    t = x - spl.xs[i]
    # Horner's polynomial evaluation
    @inbounds return c[3 * i - 2] + (c[3 * i - 1] + c[3 * i] * t) * t
end

@inline function evaluate(spl::CubicSpline, x::AbstractFloat)
    i = spl_index(spl.xs, x)
    c = spl.coeffs
    t = x - spl.xs[i]
    # Horner's polynomial evaluation
    @inbounds return c[4 * i - 3] +
                     (c[4 * i - 2] + (c[4 * i - 1] + c[4 * i] * t) * t) * t
end

@inline function bc_index(bc::AbstractBC, n::Int)
    if bc.index == :left
        return 1, 1
    elseif bc.index == :right || bc.index == n
        return n - 1, n
    elseif isa(bc.index, Int)
        si = bc.index
        xi = si
        return si, xi
    end
end

@inline function spl_index(xs::Vector{<:AbstractFloat}, x::AbstractFloat)
    i = searchsortedlast(xs, x)
    if i == 0
        return 1
    end
    return min(i, length(xs) - 1)
end

raw"""
    LinearSplineFromCoeffs

    Create a LinearSpline from given x-coordinates and coefficients.
    The coefficients should be provided in the order of the linear segments, as follows:

    ```math
        f(x) = a_1 + a_2 (x - x_i)
    ```

    Example
    -------
    ```julia
    xs = [0.0, 1.0, 2.0]
    coeffs = [1.0, 1.0, 2.0, -2.0]  # Coefficients for two segments
    spl = LinearSplineFromCoeffs(xs, coeffs)
    ```
"""
function LinearSplineFromCoeffs(
    xs::AbstractVector{<:AbstractFloat}, coeffs::AbstractVector{<:AbstractFloat}
)
    return LinearSpline(xs, [0.0], coeffs)
end

raw"""
    QuadraticSplineFromCoeffs

    Create a QuadraticSpline from given x-coordinates and coefficients.
    The coefficients should be provided in the order of the quadratic segments, as follows:

    ```math
        f(x) = a_1 + a_2 (x - x_i) + a_3 (x - x_i)^2
    ```

    Example
    -------
    ```julia
    xs = [0.0, 1.0, 2.0]
    coeffs = [1.0, 1.0, 0.5, 2.0, -2.0, 1.0]  # Coefficients for two segments
    spl = QuadraticSplineFromCoeffs(xs, coeffs)
    ```
"""
function QuadraticSplineFromCoeffs(
    xs::AbstractVector{<:AbstractFloat}, coeffs::AbstractVector{<:AbstractFloat}
)
    return QuadraticSpline(xs, [0.0], coeffs, NopBC())
end

raw"""
    CubicSplineFromCoeffs

    Create a CubicSpline from given x-coordinates and coefficients.
    The coefficients should be provided in the order of the cubic segments, as follows:

    ```math
        f(x) = a_1 + a_2 (x - x_i) + a_3 (x - x_i)^2 + a_4 (x - x_i)^3
    ```

    Example
    -------
    ```julia
    xs = [0.0, 1.0, 2.0]
    coeffs = [1.0, 1.0, 0.5, 2.0, -2.0, 1.0, 0.5, -1.0]  # Coefficients for two segments
    spl = CubicSplineFromCoeffs(xs, coeffs)
    ```
"""
function CubicSplineFromCoeffs(
    xs::AbstractVector{<:AbstractFloat}, coeffs::AbstractVector{<:AbstractFloat}
)
    return CubicSpline(xs, [0.0], coeffs, NopBC(), NopBC())
end
