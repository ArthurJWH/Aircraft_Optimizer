"""
    AbstractMesh

Abstract supertype for mesh-like objects used across the package.
Subtyping `AbstractMesh` (for example `struct VLMMesh <: AbstractMesh`)
enables consistent, type-stable dispatch for plotting and processing functions.
"""
abstract type AbstractMesh end
