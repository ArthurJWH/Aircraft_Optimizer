# Repository Structure

The repository tree consists of:

- [`MyPackage`](#MyPackage): main package with source code
- [`scripts`](#Scripts): single run codes
- [`projects`](#Projects): project specific codes
- [`assets`](#Assets): stores data and resource files
- [`docs`](#Docs): repository documentation
- [`miscellaneous`](#Miscellaneous): additional codes

## MyPackage

This is the main package with all the source code. These files consists of reusable codes, such as Struct and function definition.

The exposed Structs and functions can be accessed through:

```julia
using MyPackage
using MyPackage.Geometry
using MyPackage.Utils
...
```

## Scripts

In contrast to the source code, scripts are meant for single run. In other words, these are the codes that are not reused in other files.

## Projects

This folder contains project specific codes (e.g. objective function of an aircraft mission, analysis of a specific aircraft, outputs of an optimization).

## Assets

Assets are data or resource files that can be used in different projects. For example it can store:
- Airfoils data
- Materials data
- Static media

## Docs

It contains all project documentation, user guides, and technical references.

It also builds the website documentation through Documenter.jl.

## Miscellaneous

### Project.toml & Manifest.toml

Project.toml defines the project metadata, mainly its dependencies.

Manifest.toml tracks the full dependency tree alongside their versions, guaranteeing reproducibility across different machines.

### Setup.jl & Startup.jl

Setup.jl sets the work environment, guaranteeing that all packages are installed correctly

Startup.jl starts the Julia REPL at the beginning of a session, making sure that the correct environment is used.

### Format.jl and .JuliaFormatter.toml

Define the formatting guidelines for consistency between contributors and prevent `diff` from individual formatting.

### .githooks/pre-commit

Prevents code commit prior to formatting check and ensuring the format consistency.