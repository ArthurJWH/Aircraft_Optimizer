using Plots

using ..Geometry: Airfoil

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
