# Understanding the Repository

This guide explains the directory structure of `Aircraft_Optimizer` and provides clear rules on **where to add new files** to keep the codebase clean, modular, and scalable.

---

## Directory Tree

```text
aircraft_optimizer/
├── MyPackage/                  # Core reusable package source code
│   ├── Project.toml            # Package dependencies and UUID
│   └── src/                    # Submodules (Geometry, VLM, IO, Utils, Info, MyPlots)
├── scripts/                    # Single-run scripts & quick experiments
├── projects/                   # Aircraft-specific design studies & optimization pipelines
├── assets/                     # Shared data & resources (airfoils, coordinates, media)
├── docs/                       # Documentation site built with Documenter.jl
│   ├── make.jl                 # Documentation build entry point
│   └── src/                    # Markdown pages and documentation assets
├── Setup.jl                    # One-time environment installation script
├── Startup.jl                  # Interactive REPL session startup helper
├── Format.jl                   # Code formatting runner
├── .JuliaFormatter.toml        # Code formatting rules
└── .githooks/                  # Git hooks (e.g. pre-commit formatting checks)
```

---

## Quick Repository Guide: Where Should I Put My Code?

| What are you creating? | Target Location | Why? |
|:-----------------------|:----------------|:-----|
| A new reusable struct, algorithm, solver function, or file parser | `MyPackage/src/<Submodule>/` | Core code intended to be imported and shared across multiple analyses. |
| A one-off script, quick plot, scratch test, or exploratory calculation | `scripts/` | Standalone scripts that execute once and are not imported by other files. |
| A complete aircraft design study, mission optimization, or thesis analysis | `projects/<project_name>/` | Multi-file, project-specific workflows (e.g. specific objective functions, design iterations). |
| A new airfoil coordinate file (`.dat`), material database, or CAD resource | `assets/airfoils/` or `assets/` | Static data files accessible across all scripts and projects. |
| A new how-to tutorial, technical derivation, or guide page | `docs/src/` | Documentation compiled into the Documenter.jl website. |

---

## Detailed Directory Breakdown

### 1. `MyPackage/` — Core Reusable Library

`MyPackage` contains all **reusable, generic code**. This is where structs, math routines, meshing algorithms, and physics solvers live. Code written here should be generic and not hard-coded for a single specific aircraft.

The package is organized into focused submodules under `MyPackage/src/`:

- **`Geometry/`** (`MyPackage.Geometry`):
  - Abstract representations of 2D airfoil profiles ([`Airfoil`](@ref MyPackage.Geometry.Airfoil)), 3D lifting surfaces ([`Aerosurface`](@ref MyPackage.Geometry.Aerosurface)), and whole-aircraft assemblies ([`Plane`](@ref MyPackage.Geometry.Plane)).
- **`Analysis/VLM/`** (`MyPackage.VLM`):
  - Vortex Lattice Method meshing ([`VLMMesh`](@ref MyPackage.VLM.VLMMesh), [`VLMGeometry`](@ref MyPackage.VLM.VLMGeometry)), flow setup ([`VLMSetup`](@ref MyPackage.VLM.VLMSetup)), matrix solver ([`VLMSolver!`](@ref MyPackage.VLM.VLMSolver!)), load/coefficient post-processing ([`VLMLoad`](@ref MyPackage.VLM.VLMLoad), [`VLMCoefficients`](@ref MyPackage.VLM.VLMCoefficients)), and stability derivatives ([`VLMStabilityDerivatives`](@ref MyPackage.VLM.VLMStabilityDerivatives)).
- **`IO/`** (`MyPackage.IO`):
  - File reading ([`read_dat`](@ref MyPackage.IO.read_dat), [`airfoil_from_dat`](@ref MyPackage.IO.airfoil_from_dat)) and 3D outer-mold-line CAD JSON export ([`export_plane_json`](@ref MyPackage.IO.export_plane_json)).
- **`Utils/`** (`MyPackage.Utils`):
  - Numerical quadrature ([`IntegrateGLQ`](@ref MyPackage.Utils.IntegrateGLQ)), root finding ([`bisection`](@ref MyPackage.Utils.bisection)), piecewise splines ([`LinearSpline`](@ref MyPackage.Utils.LinearSpline), [`QuadraticSpline`](@ref MyPackage.Utils.QuadraticSpline), [`CubicSpline`](@ref MyPackage.Utils.CubicSpline)), and regression ([`LSR`](@ref MyPackage.Utils.LSR)).
- **`Info/`** (`MyPackage.PlaneInfo`):
  - Aircraft performance metadata and center-of-gravity containers ([`Data`](@ref MyPackage.PlaneInfo.Data), [`Coeffs`](@ref MyPackage.PlaneInfo.Coeffs)).
- **`MyPlots/`** (`MyPackage.MyPlots`):
  - Visualization utilities for 3D panel meshes ([`plot_mesh`](@ref MyPackage.MyPlots.plot_mesh)) and 2D airfoil profiles ([`plot_airfoil`](@ref MyPackage.MyPlots.plot_airfoil)).

!!! tip "Adding New Code to `MyPackage`"
    1. Place your new file in the appropriate submodule folder (e.g. `MyPackage/src/Geometry/my_new_feature.jl`).
    2. Include it inside that submodule's entry file (e.g. `MyPackage/src/Geometry/Geometry.jl` using `include("my_new_feature.jl")`).
    3. Export any public structs or functions you want accessible to users (e.g. `export MyNewStruct`).

---

### 2. `scripts/` — Single-Run & Scratch Scripts

The `scripts/` directory is for **standalone, single-execution files**. In contrast to `MyPackage`, code here is **not meant to be imported** by other files.

Typical contents of `scripts/`:
- Quick sanity checks and exploratory scripts (e.g. running a test polar).
- Scratchpad calculations and temporary plotting scripts.
- Single-run batch processing jobs.

```julia
# Example script: scripts/test_wing_polar.jl
using MyPackage
using MyPackage.Geometry
using MyPackage.VLM
using MyPackage.IO

# Script-specific execution logic...
```

---

### 3. `projects/` — Engineering Studies & Aircraft Designs

The `projects/` directory holds **project-specific workflows, aircraft designs, and optimization studies**.

When you are designing a specific airplane or running an optimization study with multiple interrelated files (e.g. an objective function, constraints, geometry generators, mission simulators, and result plotters), create a dedicated folder under `projects/`:

```text
projects/
├── sae_aero_2026/              # Example design study
│   ├── geometry_builder.jl     # Assembles this specific airplane
│   ├── mission_analysis.jl     # Mission payload / takeoff scoring
│   ├── optimize_wing.jl        # Optimization driver script
│   └── outputs/                # Generated polars, figures, logs
└── uav_glider/
    └── run_study.jl
```

!!! note "Projects vs. MyPackage"
    - **`MyPackage`**: Provides generic tools (e.g. *how to solve any VLM mesh*).
    - **`projects/`**: Uses those tools for a concrete application (e.g. *optimizing the 2026 competition airplane*).

---

### 4. `assets/` — Data Files & Static Resources

The `assets/` directory stores **non-code data and resource files** shared across projects:

```text
assets/
└── airfoils/
    ├── Plain/
    │   └── Plain.dat           # Flat-plate / zero-camber reference airfoil
    ├── NACA0012/
    │   └── NACA0012.dat
    └── NACA4412/
        └── NACA4412.dat
```

!!! tip "Adding New Airfoils"
    Add new airfoil coordinate files under `assets/airfoils/<AirfoilName>/<AirfoilName>.dat`. Use standard Selig or Lednicer format (coordinates start at trailing edge, proceed over upper surface to leading edge, and return along lower surface to trailing edge).

---

### 5. `docs/` — Documentation Site

Contains all source files and configuration for the Documenter.jl documentation site:
- `docs/make.jl`: Build script that parses docstrings and renders HTML.
- `docs/src/`: Markdown source files organized by topic (`onboarding/`, `how_to/`, `technical/`, `code_reference/`).

To build the documentation locally:
```bash
julia --project=docs docs/make.jl
```

---

## Repository Maintenance & Developer Tools

### Environment Configuration
- **`Project.toml`**: Defines project dependencies (`LinearAlgebra`, `StaticArrays`, `Plots`, `Makie`, `JSON3`, etc.).
- **`Manifest.toml`**: Exact snapshot of all installed package versions, ensuring 100% reproducibility across different machines.
- **`Setup.jl`**: Run once when cloning the repo to instantiate the package environment.
- **`Startup.jl`**: Run at the beginning of an interactive Julia session to activate the environment.

### Code Formatting
- **`.JuliaFormatter.toml`**: Formatting configuration (4-space indentation, 80-character margin, trailing commas).
- **`Format.jl`**: Formats all Julia files in the repository automatically:
  ```julia
  include("Format.jl")
  ```
- **`.githooks/pre-commit`**: Optional Git hook that automatically runs formatting checks before committing.