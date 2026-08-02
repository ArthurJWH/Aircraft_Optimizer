using GLMakie

using ..VLM: VLMMesh

"""
    plot_mesh(mesh::VLMMesh; min_extent = 5.0)

Interactive 3D visualization of a VLM mesh.

Arguments

  - `mesh`: the VLM mesh to plot
  - `min_extent`: minimum axis span in each direction
"""
function plot_mesh(mesh::VLMMesh; min_extent=5.0)
    return plot_mesh([mesh]; min_extent=min_extent)
end

function plot_mesh(meshes::AbstractVector{<:VLMMesh}; min_extent=5.0)
    fig, ax = _initialize_ax()

    xlims = (Inf, -Inf)
    ylims = (Inf, -Inf)
    zlims = (Inf, -Inf)

    for mesh in meshes
        _plot!(ax, mesh)
        xlims = _merge_bounds(xlims, _mesh_bounds(mesh, min_extent).x)
        ylims = _merge_bounds(ylims, _mesh_bounds(mesh, min_extent).y)
        zlims = _merge_bounds(zlims, _mesh_bounds(mesh, min_extent).z)
    end

    _apply_bounds!(ax, xlims, ylims, zlims)
    GLMakie.display(fig)

    return fig
end

function _initialize_ax()
    fig = GLMakie.Figure(; size=(1000, 800))

    ax = GLMakie.Axis3(
        fig[1, 1];
        xlabel="X",
        ylabel="Y",
        zlabel="Z",
        title="VLM Mesh",
        aspect=:data,
        perspectiveness=0.75,
    )

    return fig, ax
end

function _plot!(ax, mesh::VLMMesh)
    x = mesh.vertices[1, :, :]
    y = mesh.vertices[2, :, :]
    z = mesh.vertices[3, :, :]

    GLMakie.surface!(ax, x, y, z; shading=true, colormap=:viridis)

    if mesh.mirror_xz
        GLMakie.surface!(ax, x, -y, z; shading=true, colormap=:viridis)
        GLMakie.wireframe!(ax, x, -y, z; color=(:black, 0.4), linewidth=1)
    end

    GLMakie.wireframe!(ax, x, y, z; color=(:black, 0.4), linewidth=1)

    return nothing
end

function _mesh_bounds(mesh::VLMMesh, min_extent)
    x = mesh.vertices[1, :, :]
    y = mesh.vertices[2, :, :]
    z = mesh.vertices[3, :, :]

    xlims = _enforce_range(minimum(x), maximum(x), min_extent)
    ylims = _enforce_range(minimum(y), maximum(y), min_extent)
    zlims = _enforce_range(minimum(z), maximum(z), min_extent)

    if mesh.mirror_xz
        ymins = min(ylims[1], -ylims[2])
        ymaxs = max(ylims[2], -ylims[1])
        ylims = _enforce_range(ymins, ymaxs, min_extent)
    end

    return (x=xlims, y=ylims, z=zlims)
end

function _enforce_range(minv, maxv, minsize)
    minv = min(minv, 0.0)
    maxv = max(maxv, 0.0)
    span = maxv - minv
    if span < minsize
        center = (minv + maxv) / 2
        half = minsize / 2
        return (center - half, center + half)
    end
    return (minv, maxv)
end

function _merge_bounds(bounds1, bounds2)
    return (min(bounds1[1], bounds2[1]), max(bounds1[2], bounds2[2]))
end

function _apply_bounds!(ax, xlims, ylims, zlims)
    GLMakie.xlims!(ax, xlims)
    GLMakie.ylims!(ax, ylims)
    return GLMakie.zlims!(ax, zlims)
end
