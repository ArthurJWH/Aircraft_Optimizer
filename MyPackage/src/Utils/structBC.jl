"""
    AbstractBC

    Abstract type for boundary conditions applied in Spline.
    Used in function declarations to support different BC types.
"""
abstract type AbstractBC end

"""
    NopBC

    Struct for no boundary condition applied. Used as a filler.

    Example
    -------
    ```julia
    bc = NopBC()  # No boundary condition
    ```
"""
struct NopBC <: AbstractBC end

"""
    FirstDerivativeBC

    Struct for boundary condition for the first derivative of a function.
    Defines the value of the first derivative at a given index.

    Fields
    ------
    value : Float64
        The value of the first derivative at the given index.
    index : T
        The index at which the boundary condition is applied.

    Example
    -------
    ```julia
    bc = FirstDerivativeBC(0.0, 1)  # First derivative is zero at index 1
    ```
"""
struct FirstDerivativeBC{T} <: AbstractBC
    value::Float64
    index::T
end

"""
    SecondDerivativeBC

    Boundary condition for the second derivative of a function.
    Defines the value of the second derivative at a given index.

    value : Float64
        The value of the second derivative at the given index.
    index : T
        The index at which the boundary condition is applied.

    Example
    -------
    ```julia
    bc = SecondDerivativeBC(0.0, 1)  # Second derivative is zero at index 1
    ```
"""
struct SecondDerivativeBC{T} <: AbstractBC
    value::Float64
    index::T
end

"""
    ThirdDerivativeBC

    Boundary condition for the third derivative of a function.
    Defines the value of the third derivative at a given index.

    value : Float64
        The value of the third derivative at the given index.
    index : T
        The index at which the boundary condition is applied.

    Example
    -------
    ```julia
    bc = ThirdDerivativeBC(0.0, 1)  # Third derivative is zero at index 1
    ```
"""
struct ThirdDerivativeBC{T} <: AbstractBC
    value::Float64
    index::T
end
