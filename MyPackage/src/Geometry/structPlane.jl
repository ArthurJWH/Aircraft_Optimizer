using ..PlaneInfo

struct Plane
    surfaces::Vector{<:Aerosurface}
    coeffs::Coeffs
    data::Data
end

function Plane(
    surfaces::Vector{<:Aerosurface}; CG::NTuple{3, Float64}=(0.0, 0.0, 0.0)
)
    coeffs = Coeffs()

    if CG === (0.0, 0.0, 0.0)
        s = surfaces[1]
        pos = s.pos
        CG = (pos[1] + 0.25 * s.MAC, pos[2], pos[3])
    end

    data = Data(; CG=CG)

    return Plane(surfaces, coeffs, data)
end
