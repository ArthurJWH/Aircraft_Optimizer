dir = abspath(@__DIR__)

pkgs = ["Plots", "Makie"]

reset = true

if reset
    for file in ("Project.toml", "Manifest.toml")
        path = joinpath(dir, file)

        isfile(path) && rm(path)
    end
end

using Pkg
Pkg.activate(dir)
Pkg.add("Revise")
using Revise

Pkg.add(; name="JuliaFormatter", version="1.0.39")

Pkg.add(pkgs)

Pkg.develop(; path=joinpath(dir, "MyPackage"))

Pkg.instantiate()

println("Setup complete for: $dir")
