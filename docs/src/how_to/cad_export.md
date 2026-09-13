# CAD Export

After aerodynamic analysis, you can export the aircraft outer mold line (OML) as a JSON
file for reconstruction in CAD software such as SolidWorks.

---

## Exporting the Geometry

Use [`export_plane_json`](@ref MyPackage.IO.export_plane_json) from the `IO` module:

```julia
using MyPackage.IO

export_plane_json(plane, "my_aircraft_oml.json"; n_chord=50, n_span=40)
```

| Parameter | Description | Default |
|:----------|:------------|:--------|
| `n_chord` | Number of chordwise points per loft station | `50` |
| `n_span` | Number of spanwise loft stations per semi-span | `40` |
| `min_thickness_rel` | Minimum relative thickness fraction | `0.002` |

!!! tip "Resolution independence"
    The loft station count (`n_span`) is independent of the VLM panel mesh
    and of how many airfoils you defined on the surface. You can use a coarse
    VLM mesh for fast analysis and a fine loft for smooth CAD surfaces.

---

## JSON Structure

The exported file contains:

```
{
  "units": "meters",
  "surfaces": [
    {
      "name": "MainWing",
      "vertical": false,
      "mirror_xz": true,
      "stations": [ ... ],         // starboard loft sections
      "guides": { ... },           // 3D guide curves
      "mirrored_stations": [ ... ] // port loft sections (if mirror_xz)
    }
  ]
}
```

### Loft Stations

Each station is a spanwise cross-section containing:
- `y_frac` — normalised span fraction
- `chord` — local chord length, in `m`
- `twist_deg` — local twist angle, in `deg`
- `points` — continuous contour from lower TE → LE → upper TE
- `top_points` / `bottom_points` — separate upper and lower surfaces

### Guide Curves

Three guide curves are generated per surface to prevent loft twisting in CAD:
- **Leading edge** — ``x = 0`` line
- **Trailing edge** — TE midpoint line
- **Reference line** — the `sw_center` sweep reference line

---

## Zero-Thickness Airfoil Protection

Flat-plate or symmetric airfoils with negligible thickness (e.g. a "Plain" airfoil)
would produce zero-volume solids that fail to loft in CAD.

The export pipeline automatically applies a minimum aerodynamic thickness envelope
controlled by `min_thickness_rel`:

```math
t_\text{eff}(x) = \max\bigl(t_\text{actual}(x),\; t_\text{min} \cdot 2\sqrt{x}(1 - 0.9x)\bigr)
```

This produces a thin but physically realisable airfoil shape with zero gap at the LE
and a finite trailing edge.

!!! note
    This protection only activates when the maximum section thickness is below
    `min_thickness_rel`. Standard cambered or thick airfoils are exported as-is.

---

## Workflow Example

```julia
using MyPackage.IO
using MyPackage.Geometry
using MyPackage.VLM

# 1. Define and analyse the aircraft (see First Simulation)
airfoil = airfoil_from_dat("assets/airfoils/NACA4412/NACA4412.dat")
wing = Aerosurface(name="Wing", airfoils=[airfoil, airfoil], b=10.0, chord=y -> 1.5)
plane = Plane([wing])

# 2. Export the OML
export_plane_json(plane, "wing_oml.json"; n_chord=60, n_span=50)
```

The resulting JSON can be read by a SolidWorks macro or any CAD tool that supports
loft-from-sections with guide curves.
