# Pressure and traction in adaptive 3-D studies

Run the physical-size pressure/shear example through the normal study command:

```bash
cargo run -p fs-cli --features sdf3-study -- --json study \
  examples/marquee/bracket-3d-pressure.fsim /tmp/pressure-3d.db
```

This is a numerical design example with explicit 70 GPa modulus, 100 x 50 mm
footprint and reference loading. It is not a validated material/part model or
a stress-qualified design. The height field and design boundary stay fixed.

Use `:loads` instead of `:body-loads` to declare mixed independent cases:

```lisp
(scenario
  :fixed-boundary left
  :loads (
    (load :body-n-m3 (0.0 0.0 -26000.0)
      :surface (pressure :pa 1000.0 :x-fraction (0.5 1.0)) :weight 0.7)
    (load :body-n-m3 (0.0 0.0 -26000.0)
      :surface (traction :pa (0.0 1000.0 0.0) :x-fraction (0.5 1.0)) :weight 0.3)))
```

Every field is required and ordered. An explicitly absent surface is `:surface
none`; an absent body contribution is the zero vector. Both absent refuses.
Each surface law must be nonzero. Legacy `:body-loads` remains unchanged and
incurs no surface quadrature. The new load family works on both the original
unit cube and the opt-in physical-size domain.

Pressure is in pascals, positive INTO the solid (`traction = -pressure * normal`).
Negative pressure is outward tension. Traction is a global Cartesian vector in
N/m² per actual reference surface area, not a total force. Body density is N/m³
per reference volume. It is NOT automatically scaled by the design density:
this is prescribed dead loading, not density-dependent self-weight.

The loaded surface is the retained ZERO LEVEL SET of the implicit graph. Box
faces that clip the domain are not part of this surface. The `x-fraction`
interval selects a patch along x across the full y extent; `(0.5 1.0)` selects
the right half of the footprint. Patch endpoints must be initial-octree planes:
multiples of 1/2 at initial level 1, or 1/4 at level 2. This keeps discontinuous
patch boundaries aligned through every subsequent refinement. There is no
silent snapping, sampled load renormalization, or conversion into a body force.

Within a load case the body and surface contributions add before equilibrium.
Different cases remain separate solves; their compliances are weighted without
normalizing weights. Opposing loads in DIFFERENT cases cannot cancel each
other. Reports retain each case's compliance and separate physical displacement
field using the existing report and design artifacts.

The initial, enriched and proposed operators all retain oriented surface rules
from the same implicit field used for stiffness. All bulk/surface points share
the original cumulative quadrature budget, including failed preparations.
The identical reference law feeds both equilibrium and the enriched weak goal
residual. No nodal-load interpolation is used. Existing volume restoration,
gradient checks, accepted-state rollback, durable stage checkpoints and
same-executable resume apply unchanged. `--budget 1` caps new stages, not loads.

Current limits remain four independent cases, one pressure OR traction patch
plus an optional body contribution per case, fixed reference geometry and
homogeneous box clamps. These are linear-elastic numerical comparisons, not
follower pressure, deformed-surface loading, continuum stress bounds or
manufacturing certification. Surface quadrature is numerical, not a certified
surface-integral enclosure.

Focused native checks: `cargo test -p fs-cli --features sdf3-study --lib
study::elasticity::sdf3`. New tests cover physical force/moment totals, real
curved-surface loading on two grids, mixed RHS assembly, pressure/traction
agreement, shared quadrature stops, and exact two-stage checkpoint recovery.
