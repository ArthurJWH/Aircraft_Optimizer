using ..Utils
using ..Geometry

struct VLMMesh
    vertices::Array{Float64, 3}
    mirror_xz::Bool
end

function VLMMesh(plane::Plane, n_chordxspan::Vector{NTuple{2, Int}})
    surfaces = plane.surfaces
    n = length(surfaces)
    ncs = length(n_chordxspan)
    meshes = Vector{VLMMesh}(undef, n)
    for i in 1:n
        n_chord, n_span = n_chordxspan[(i - 1) % ncs + 1]
        meshes[i] = VLMMesh(surfaces[i], n_chord, n_span)
    end
    return meshes
end

function VLMMesh(plane::Plane, n_chord::Int, n_span::Int)
    return [VLMMesh(surface, n_chord, n_span) for surface in plane.surfaces]
end

function VLMMesh(surface::Aerosurface, n_chord::Int, n_span::Int)
    vertices = _generate_vertices(
        surface, n_chord, n_span, Val(surface.vertical)
    )
    return VLMMesh(vertices, surface.mirror_xz)
end

@inline function _generate_vertices(
    surface::Aerosurface, n_chord::Int, n_span::Int, ::Val{false}
)
    geom = _generate_geom(surface, n_chord, n_span)

    vertices = Array{Float64, 3}(undef, 3, n_span + 1, n_chord + 1)

    @inbounds for i in eachindex(geom.x)
        xi = geom.x[i]
        spl_zi = geom.splines_z[i]
        for j in eachindex(geom.y)
            y = geom.y[j]
            chord_j = geom.chords[j]
            vertices[1, j, i] = xi * chord_j
            vertices[2, j, i] = y * geom.b / 2
            vertices[3, j, i] = spl_zi(y) * chord_j
        end
    end

    @inbounds for j in eachindex(geom.y)
        twist_j = geom.twists[j]
        sweep_length_j = geom.sweep_lengths[j]
        dihedral_length_j = geom.dihedral_lengths[j]

        _transform_section!(
            vertices,
            j,
            twist_j,
            geom.tw_center * geom.chords[j],
            sweep_length_j,
            geom.sw_center * (geom.root_chord - geom.chords[j]),
            dihedral_length_j,
            surface.pos,
        )
    end

    return vertices
end

@inline function _generate_vertices(
    surface::Aerosurface, n_chord::Int, n_span::Int, ::Val{true}
)
    geom = _generate_geom(surface, n_chord, n_span)

    vertices = Array{Float64, 3}(undef, 3, n_span + 1, n_chord + 1)

    @inbounds for i in eachindex(geom.x)
        xi = geom.x[i]
        spl_zi = geom.splines_z[i]
        for j in eachindex(geom.y)
            y = geom.y[j]
            chord_j = geom.chords[j]
            vertices[1, j, i] = xi * chord_j
            vertices[3, j, i] = y * geom.b / 2
            vertices[2, j, i] = spl_zi(y) * chord_j
        end
    end

    @inbounds for j in eachindex(geom.y)
        twist_j = geom.twists[j]
        sweep_length_j = geom.sweep_lengths[j]
        dihedral_length_j = geom.dihedral_lengths[j]

        _transform_section_v!(
            vertices,
            j,
            twist_j,
            geom.tw_center * geom.chords[j],
            sweep_length_j,
            geom.sw_center * (geom.root_chord - geom.chords[j]),
            dihedral_length_j,
            surface.pos,
        )
    end

    return vertices
end

@inline function _generate_geom(surface::Aerosurface, n_chord::Int, n_span::Int)
    b = surface.b
    ys = surface.ys
    cambers = getfield.(surface.airfoils, :camber)
    chord = surface.chord
    twist = surface.twist
    tw_center = surface.tw_center
    sweep = surface.sweep
    sw_center = surface.sw_center
    dihedral = surface.dihedral

    x = (1 .- cos.(range(0; stop=pi, length=n_chord + 1))) ./ 2
    y = (1 .- cos.(range(0; stop=pi, length=n_span + 1))) ./ 2

    chords = chord.(y)

    splines_z = Vector{LinearSpline}(undef, length(x))

    tmp = Vector{Float64}(undef, length(ys))

    @inbounds for i in eachindex(x)
        xi = x[i]

        for j in eachindex(ys)
            tmp[j] = cambers[j](xi)
        end

        splines_z[i] = LinearSpline(ys, copy(tmp))
    end

    sweep_integral = IntegrateGLQ(t -> tand(sweep(t)); n=5)
    dihedral_integral = IntegrateGLQ(t -> tand(dihedral(t)); n=5)
    sweep_lengths = Array{Float64}(undef, length(y))
    dihedral_lengths = Array{Float64}(undef, length(y))
    twists = Array{Float64}(undef, length(y))
    @inbounds for j in eachindex(y)
        yj = y[j]
        sweep_lengths[j] = b * sweep_integral(0.0, yj) / 2
        dihedral_lengths[j] = b * dihedral_integral(0.0, yj) / 2
        twists[j] = twist(yj)
    end

    return (
        b=b,
        ys=ys,
        tw_center=tw_center,
        sw_center=sw_center,
        x=x,
        y=y,
        chords=chords,
        root_chord=chords[1],
        splines_z=splines_z,
        sweep_lengths=sweep_lengths,
        dihedral_lengths=dihedral_lengths,
        twists=twists,
    )
end

@inline function _transform_section!(
    vertices, j, twist, tw_center, sweep_length, offset, dihedral_length, pos
)
    c = cosd(twist)
    s = sind(twist)

    @inbounds @simd for i in axes(vertices, 3)
        x0 = vertices[1, j, i]
        z0 = vertices[3, j, i]

        x = x0 - tw_center

        vertices[1, j, i] =
            x * c + z0 * s + tw_center + offset + sweep_length + pos[1]

        vertices[2, j, i] += pos[2]

        vertices[3, j, i] = -x * s + z0 * c + dihedral_length + pos[3]
    end

    return nothing
end

@inline function _transform_section_v!(
    vertices, j, twist, tw_center, sweep_length, offset, dihedral_length, pos
)
    c = cosd(twist)
    s = sind(twist)

    @inbounds @simd for i in axes(vertices, 3)
        x0 = vertices[1, j, i]
        z0 = vertices[2, j, i]

        x = x0 - tw_center

        vertices[1, j, i] =
            x * c + z0 * s + tw_center + offset + sweep_length + pos[1]

        vertices[2, j, i] = -x * s + z0 * c + dihedral_length + pos[2]

        vertices[3, j, i] += pos[3]
    end

    return nothing
end
