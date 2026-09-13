# Ground Effect

When an aircraft flies close to a surface (takeoff, landing, or ground-effect vehicles),
downwash is suppressed and effective lift increases. This page shows how to model ground
proximity in the VLM solver.

---

## Enabling Ground Effect

Pass `ground = true` and specify the height `h` (CG above ground) when creating
[`VLMSetup`](@ref):

```julia
setup_ge = VLMSetup(geom, 25.0; rho=1.225, ground=true, h=1.5)

VLMSolver!(setup_ge, collect(0.0:2.0:10.0), [0.0])
```

!!! note "Height definition"
    `h` is the vertical distance from the aircraft centre of gravity to the ground
    plane, **not** the wing clearance or landing gear height.

---

## Comparing with Free Flight

Run the same geometry and speed without ground effect for a direct comparison:

```julia
setup_free = VLMSetup(geom, 25.0; rho=1.225, ground=false)
VLMSolver!(setup_free, collect(0.0:2.0:10.0), [0.0])

polar_ge   = VLMCoefficientsSlice(setup_ge;   beta=0.0)
polar_free = VLMCoefficientsSlice(setup_free; beta=0.0)

using Plots
Plots.plot(polar_free.alpha, polar_free.CL_total,
    label = "Free Flight", lw = 2,
)
Plots.plot!(polar_ge.alpha, polar_ge.CL_total,
    label = "Ground Effect (h = 1.5 m)", lw = 2, ls = :dash,
)
Plots.xlabel!("α [deg]")
Plots.ylabel!("CL")
```

You should observe higher ``C_L`` and lower induced drag in ground effect at the same
angle of attack.

---

## How It Works

The solver models the ground plane at ``z = -h`` by creating an **image vortex system**
— a mirror copy of all vortex rings and trailing wake filaments reflected across the
ground plane. This enforces the non-penetration boundary condition (zero normal velocity
at the ground surface) via the method of images.

The image system is encapsulated in [`GroundTransform`](@ref), which precomputes the
necessary trigonometric coefficients for the reflection.

!!! warning "Performance note"
    Ground effect disables the ring geometry cache because the image vortex system
    depends on the angle of attack. Each solved point requires a full AIC assembly,
    making ground-effect sweeps slower than free-flight sweeps.

!!! warning "Limitations"
    - The ground plane is assumed to be flat and infinite.
    - The model does not account for ground boundary layer or runway roughness effects.
    - Very low heights (``h \to 0``) may produce numerical instabilities.
