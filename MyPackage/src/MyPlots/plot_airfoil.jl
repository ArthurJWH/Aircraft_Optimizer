using Plots

using ..Geometry: Airfoil

"""
    plot_airfoil(airfoil::Airfoil; save=false)

Plots the upper and lower surfaces of the given [`Airfoil`](@ref MyPackage.Geometry.Airfoil) alongside its mean camber line using `Plots.jl`.

# Arguments
- `airfoil::Airfoil`: The airfoil object to plot.
- `save::Bool`: Whether to save the plot as a PNG image adjacent to the source `.dat` file (default: `false`).

# Returns
- `p::Plots.Plot`: The generated plot object.

# Example
```julia
using Plots
airfoil = airfoil_from_dat("assets/airfoils/Plain/Plain.dat")
p = plot_airfoil(airfoil; save=false)
display(p)
```
"""
function plot_airfoil(airfoil::Airfoil; save=false)
    x = range(0.0, 1.0; length=100)
    y_top = airfoil.top_surface.(x)
    y_bottom = airfoil.bottom_surface.(x)
    y_camber = airfoil.camber.(x)

    title = airfoil.name

    p = Plots.plot(
        x,
        y_top;
        label="Top Surface",
        title=title,
        xlabel="x",
        ylabel="y",
        aspect_ratio=:equal,
    )
    Plots.plot!(p, x, y_bottom; label="Bottom Surface")
    Plots.plot!(p, x, y_camber; label="Camber Line", linestyle=:dash)

    if save
        save_path = replace(airfoil.datfile, ".dat" => ".png")
        Plots.savefig(p, save_path)
    end

    return p
end
