# Code Reference

This page provides the API and type reference for `MyPackage`, automatically generated from docstrings across the codebase.

The package is partitioned into focused, composable submodules:

| Submodule | Purpose |
|:----------|:--------|
| [`MyPackage.Geometry`](#Geometry) | Aircraft representation, parametric lifting surfaces, and 2D airfoil profiles. |
| [`MyPackage.VLM`](#VLM-Solver-&-Aerodynamics) | Vortex Lattice Method mesher, AIC solver, forces, moments, polars, and stability derivatives. |
| [`MyPackage.IO`](#Input-/-Output) | File parsers (`.dat` airfoils) and CAD outer-mold-line (OML) JSON exporters. |
| [`MyPackage.Utils`](#Numerical-Utilities) | Gauss-Legendre quadrature, spline interpolation, and least-squares regression. |
| [`MyPackage.PlaneInfo`](#Aircraft-Information-&-State) | Lightweight data and performance containers. |
| [`MyPackage.MyPlots`](#Visualization) | 3D mesh rendering (`Makie.jl`) and 2D airfoil plotting (`Plots.jl`). |

---

## Geometry

The `Geometry` module handles continuous parametric representations of 3D aircraft configurations, aerodynamic surfaces, and 2D airfoil sections.

```@autodocs
Modules = [MyPackage.Geometry]
```

---

## VLM (Solver & Aerodynamics)

The `VLM` module implements the 3D Vortex Lattice Method solver with realigning trailing wake filaments, thread-safe polar storage, far-field Trefftz drag integration, ground effect image systems, and linear stability derivative evaluation.

```@autodocs
Modules = [MyPackage.VLM]
```

---

## Input / Output

The `IO` module provides utilities for reading `.dat` airfoil coordinate files and exporting 3D surface geometry into JSON for downstream CAD and SolidWorks lofting.

```@autodocs
Modules = [MyPackage.IO]
```

---

## Numerical Utilities

The `Utils` module implements high-performance numerical routines tailored for aerodynamic geometry discretization and analysis, including Gauss-Legendre numerical quadrature, piecewise polynomial splines, and least-squares curve fitting.

```@autodocs
Modules = [MyPackage.Utils]
```

---

## Aircraft Information & State

The `PlaneInfo` module provides containers storing performance metadata, center-of-gravity locations, and whole-aircraft coefficient states.

```@autodocs
Modules = [MyPackage.PlaneInfo]
```

---

## Visualization

The `MyPlots` module provides visualization routines for inspectable 3D wireframe panel meshes and 2D airfoil profiles.

```@autodocs
Modules = [MyPackage.MyPlots]
```