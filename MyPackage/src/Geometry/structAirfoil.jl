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
