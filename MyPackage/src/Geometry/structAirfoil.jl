using ..IO
using ..Utils

"""
    Airfoil

    A struct that represents an airfoil.
    It includes the name, .dat file path, top and bottom surfaces, and camber line.

    Fields
    ------
    name : SubString{String}
        The name of the airfoil, extracted from the dat file name.
    datfile : String
        The path to the .dat file containing the airfoil coordinates.
    top_surface : T
        A function representing the top surface of the airfoil.
    bottom_surface : B
        A function representing the bottom surface of the airfoil.
    camber : C
        A function representing the camber line of the airfoil.
"""
struct Airfoil{T, B, C}
    name::SubString{String}
    datfile::String
    top_surface::T
    bottom_surface::B
    camber::C
end

function Airfoil(datfile::String)
    name = split(basename(datfile), ".")[1]
    println("Loading airfoil: $name from $datfile")
    airfoil_data = read_dat(datfile)
    top_surface, bottom_surface = calc_surfaces(airfoil_data)
    camber = calc_camber(top_surface, bottom_surface)
    return Airfoil(name, datfile, top_surface, bottom_surface, camber)
end

function calc_surfaces(airfoil_data)
    # find index of the leading-edge (minimum x)
    le_index = argmin(airfoil_data[:, 1])
    le_x = airfoil_data[le_index, 1]

    if le_x > 0.01
        # leading-edge not at x≈0 — prepend a (0,0) point
        top_data = vcat([(0.0, 0.0)], @views airfoil_data[le_index:-1:1, :])
        bottom_data = vcat([(0.0, 0.0)], @views airfoil_data[le_index:end, :])
    else
        top_data = @views airfoil_data[le_index:-1:1, :]
        bottom_data = @views airfoil_data[le_index:end, :]
    end

    # TODO: B-Spline top and bottom surfaces
    top_surface = LinearSpline(top_data[:, 1], top_data[:, 2])
    bottom_surface = LinearSpline(bottom_data[:, 1], bottom_data[:, 2])
    return top_surface, bottom_surface
end

function calc_camber(top_surface, bottom_surface)
    camber_func = x -> (top_surface(x) + bottom_surface(x)) / 2
    return camber_func
end

function calc_camber(datfile::String)
    airfoil_data = read_dat(datfile)
    top_surface, bottom_surface = calc_surfaces(airfoil_data)
    camber_func = x -> (top_surface(x) + bottom_surface(x)) / 2
    return camber_func
end
