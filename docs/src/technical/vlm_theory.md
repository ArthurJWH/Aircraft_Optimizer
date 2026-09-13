# VLM Theory

Mathematical and numerical foundations of the Vortex Lattice Method solver.

---

## Governing Assumptions

The VLM formulation operates under standard potential flow theory:

1. **Inviscid & incompressible** — zero viscosity, constant density.
2. **Thin surface approximation** — lifting surfaces are represented by their mean
   camber surfaces (thickness is not modelled in the aerodynamic solve).
3. **Attached flow** — no stall, separation, or boundary layer effects.
4. **Semi-infinite trailing wake** — trailing vortex filaments extend downstream aligned
   with the local freestream direction.

!!! warning
    VLM cannot predict viscous drag, stall behaviour, or transonic/supersonic effects.
    It solves for inviscid lift, induced drag, and aerodynamic moments only.

---

## Panel Layout — The Pistolesi Theorem

Each panel on the camber surface is modelled as a closed quadrilateral **vortex ring**
([`VortexRing`](@ref)) of unknown circulation ``\Gamma_j``:

- **Bound vortex** at the panel's **quarter-chord** line (25% of panel chord).
- **Collocation point** at the panel's **three-quarter-chord** line (75%), centred
  spanwise.
- **Surface normal** ``\vec{n}_i`` from the cross product of panel diagonals:

```math
\vec{n}_i = \frac{(D - B) \times (C - A)}{\|(D - B) \times (C - A)\|}
```

The 1/4–3/4 placement (Pistolesi's theorem) reproduces the exact 2D thin-airfoil
lift-curve slope ``C_{l\alpha} = 2\pi`` per radian with a single chordwise panel.

---

## Cosine Clustering

Panel edges are distributed using cosine spacing in both chordwise and spanwise
directions:

```math
x_i = \frac{1 - \cos(\pi i / N_c)}{2}, \qquad y_j = \frac{1 - \cos(\pi j / N_s)}{2}
```

This concentrates panels at the leading edge, trailing edge, wing root, and wingtips —
regions where circulation gradients are steepest and accuracy is most sensitive to
resolution.

---

## Biot-Savart Law with Core Regularisation

The velocity induced at a field point ``P`` by a straight vortex filament of unit
circulation running from ``A`` to ``B`` is:

```math
\vec{v}(P) = \frac{\Gamma}{4\pi}
\frac{\vec{r}_1 \times \vec{r}_2}
     {\|\vec{r}_1 \times \vec{r}_2\|^2 + \epsilon^2 r_0^2}
\left[
  \frac{\vec{r}_0 \cdot \vec{r}_1}{\sqrt{r_1^2 + \epsilon^2 r_0^2}}
  - \frac{\vec{r}_0 \cdot \vec{r}_2}{\sqrt{r_2^2 + \epsilon^2 r_0^2}}
\right]
```

where ``\vec{r}_0 = B - A``, ``\vec{r}_1 = P - A``, ``\vec{r}_2 = P - B``.

The dimensionless parameter ``\epsilon^2`` (default `1e-10` in [`VLMSetup`](@ref))
prevents singularities when collocation points lie close to vortex filaments.

---

## Wind Direction and Trailing Wake

The freestream wind direction is defined by angle of attack ``\alpha`` and sideslip
``\beta``:

```math
\vec{V}_\text{dir}(\alpha, \beta) = \begin{bmatrix}
  \cos\alpha \cos\beta \\
  -\sin\beta \\
  \sin\alpha \cos\beta
\end{bmatrix}
```

Trailing-edge vortex filaments extend semi-infinitely along ``\vec{V}_\text{dir}`` to
model the downstream wake. See [Sign Conventions](sign_conventions.md) for the coordinate system.

!!! note "Wake alignment vs AVL"
    Unlike AVL (where trailing vortices are fixed parallel to the body X-axis),
    this solver aligns trailing wake filaments with the true local wind direction
    at each ``(\alpha, \beta)``. This is more physically accurate for swept and
    high-dihedral configurations, but it means the AIC matrix depends on the flow
    angles.

---

## Kinematic Boundary Condition

The total flow must be tangent to the camber surface at every collocation point:

```math
(\vec{V}_\infty + \vec{v}_{\text{ind}, i}) \cdot \vec{n}_i = 0
```

Substituting the induced velocity summation
``\vec{v}_{\text{ind}, i} = \sum_{j=1}^{N} \vec{AIC}_{i,j} \Gamma_j``:

```math
\sum_{j=1}^{N} (\vec{AIC}_{i,j} \cdot \vec{n}_i)\, \Gamma_j
  = -(\vec{V}_\text{dir}\, V_\infty) \cdot \vec{n}_i
```

This yields the dense ``N \times N`` linear system:

```math
[\text{AIC}]\, \{\Gamma\} = \{\text{RHS}\}
```

The AIC matrix is factored via LU decomposition:
``\{\Gamma\} = \text{lu}(\text{AIC}) \backslash \text{RHS}``.

---

## Force Evaluation

Once panel circulations are known, local bound circulation differences across panel
boundaries are computed:

```math
\Delta\Gamma_c = \Gamma_j - \Gamma_\text{chordwise-upstream}, \qquad
\Delta\Gamma_s = \Gamma_j - \Gamma_\text{spanwise-inboard}
```

The local aerodynamic force on each panel is calculated via the **Kutta-Joukowski
theorem**:

```math
\vec{F}_j = \rho\, (\vec{V}_{\text{total},j} \times \vec{\ell}_j)
```

where ``\vec{V}_{\text{total},j} = \vec{V}_\infty + \vec{v}_{\text{ind},j}``.

---

## Near-Field Drag vs Trefftz Plane Drag

The solver computes drag in two independent ways:

### Near-field drag (``D``)

Obtained by projecting local panel forces onto the freestream direction. Sensitive to
**chordwise panel density** because it integrates small force components over discretely
angled panels.

### Trefftz plane drag (``D_\text{trefftz}``)

Evaluated by integrating the kinetic energy of the trailing vortex sheet on a 2D plane
infinitely far downstream, perpendicular to ``\vec{V}_\text{dir}``:

```math
D_\text{trefftz} = \rho \sum_{t=1}^{N_\text{TE}}
  \Gamma_t\, (u_\text{ind}\, \Delta v_t - v_\text{ind}\, \Delta u_t)
```

This depends only on trailing-edge circulation differences and spanwise wake geometry.
It is **immune to chordwise discretisation artifacts** and provides the true,
mesh-converged vortex-induced drag.

!!! tip "Best practice"
    Always prefer Trefftz drag (`CD_trefftz`) for reporting induced drag values.
    Use near-field drag primarily for convergence studies or debugging.

---

## Ground Effect — Method of Images

When ground effect is enabled ([`VLMSetup`](@ref) with `ground = true`), the flat ground
plane at ``z = -h`` is modelled via the **method of images**: a mirror copy of all
vortex rings and trailing wake filaments is created below the ground plane.

The image system cancels normal velocity at the ground surface, enforcing the
non-penetration boundary condition. The reflection transformation is encapsulated in
[`GroundTransform`](@ref).

!!! note
    Ground effect disables the angle-independent ring geometry cache because the image
    positions depend on ``\alpha``.
