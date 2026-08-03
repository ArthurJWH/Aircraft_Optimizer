using Plots

using ..Geometry: Airfoil

"""
    plot_airfoil(airfoil::Airfoil; save=false)

    Plots the top and bottom surfaces of the given `Airfoil` along with its camber line.

    Arguments
    ---------
    airfoil : Airfoil
        The airfoil to plot.
    save : bool, optional
        Whether to save the plot as a PNG file (default is False).

    Returns
    -------
    p : Plot
        The generated plot object.

    Example
    -------
    ```julia
    airfoil = Airfoil("NACA2412", "path/to/NACA2412.dat")
    p = plot_airfoil(airfoil; save=true)
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
