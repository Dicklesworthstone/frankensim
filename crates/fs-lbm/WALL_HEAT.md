# Equilibrium turbulent wall heat transfer

`conjugate::turbulence::thermal_wall` supplies a constant-property,
Prandtl-aware thermal wall model and couples it to the existing conservative
finite-volume energy equation. It is opt-in: existing LVEL and SGS runs do not
change automatically.

## Complete native cooling example

```sh
cargo run --release -p fs-lbm --example lvel_thermal_wall
```

The example solves native LVEL flow over a conducting slab, deposits 3 W in
the slab, solves fluid and solid temperatures together, and prints wall heat
transfer, maximum slab temperature and the original global heat balance.
Its small mesh illustrates the application path, not grid-converged turbulence.

## Model and ownership

For the existing Spalding profile, let `s = dy+/du+ = 1 + nu_t/nu`. The
constant-flux thermal wall equation is

```text
T+/Pr = integral_0^u+ s / [1 + (Pr/Pr_t)(s - 1)] du+
R_fluid = (distance/k) * (T+/Pr)/y+
```

The actual tangential wall-relative speed and wall-normal distance determine
`u+` and `y+`. The molecular Prandtl number comes from the supplied fluid;
turbulent Prandtl, quadrature tolerance and work limit are explicit controls.
At rest `R_fluid = distance/k` exactly. At `Pr = Pr_t` the exact Reynolds-analogy
limit is `T+ = Pr*u+`. Positive exponential remainders avoid sublayer
cancellation. Adaptive composite Simpson integration retains a successive-rule
convergence estimate; that estimate is not an interval error certificate.

`prepare_wall_heat` takes an accepted `FvFlow`, the original `ThermalSetup`,
and explicit `WallHeatFace` selections. Each selection is a fluid cell and its
face toward an internal stationary solid. The native fluid conductivity uses
`rho cp nu_t / Pr_t` from the tangent eddy field. The selected fluid half-cell
resistance is replaced by the integrated wall resistance by adding only the
missing nonnegative sublayer resistance to the original face coefficient.
Solid conduction and any physical contact stay in series. Temperatures on both
sides remain unknowns of the same solve. No external heat sink or second heat
source is introduced for an internal wall. An empty internal selection prepares
bulk tangent transport only, allowing composition with exterior walls.

`solve_energy_with_wall_heat` performs this preparation and calls the original
energy solver. `EnergySolution::solid_to_fluid_heat_w` and global balance consume
the same retained face resistance. Each wall report distinguishes the physical
contact, native fluid half-cell and wall-model correction.

A prepared setup can also be passed to the existing transient energy driver
while that velocity/eddy field is held fixed. When flow changes, prepare again
from the **original** setup; applying preparation to its own output refuses.
No new flow integrator, clock or temperature storage is introduced.

## Exterior walls and wall temperatures

`prepare_wall_heat_boundaries` and `solve_energy_with_wall_heat_boundaries`
compose the internal selection with explicit `BoundaryWallHeat` selections.
Each exterior selection gives its fluid cell, outward face, a Temperature,
Convective or inward HeatFlux condition, and independent physical contact
resistance. Mechanical face rules must identify a no-slip Wall with finite
purely tangential motion. Friction uses fluid velocity minus wall velocity;
normal velocity never contributes to the wall-law friction calculation.

The original thermal face must be Adiabatic, explicitly delegating transport
on the selected cell faces. No other cells on that side are reconfigured.
Temperature and film conditions contribute the exact native linear exchange
`G = area / (R_fluid + R_contact + 1/h)` (omit `1/h` for Temperature).
There is no second fluid half-cell resistance. Prescribed heat flux contributes
its exact watts, independently of the wall resistance. The solved cell field
then gives the fluid-side wall temperature, including its modeled temperature
drop. Coating resistance and external films have no storage in this model.

As with the existing `boundary_exchange` adapter, the native energy balance
counts delegated temperature/film heat in `sink_outflow_w` and imposed inward
wall flux in `source_w`. `boundary_exchanges` disaggregates those exact terms
as signed outward watts and surface temperatures. They must NOT be added again
to the native energy balance. Internal interfaces remain equal/opposite face
fluxes. Caller-authored sources, sinks and contacts are preserved.

A per-preparation cache shares identical exterior speed/distance evaluations
using exact floating-point keys, not rounding or interpolation. The returned
quadrature count reports actual work, while each face retains its own condition,
area, contact and coefficient. The existing equilibrium momentum model is not
modified; a thermal wall-relative calculation is not a claim that arbitrary
supplied momentum fields have been validated in a moving frame.

## Scope and refusals

This is a smooth-wall equilibrium closure suitable for explicitly selected
LVEL cooling studies. It does not establish the origin or validity of an
arbitrary supplied flow field. It does not model separation, roughness,
variable fluid properties, viscous heating or resolved-wall LES. Open faces,
symmetry boundaries and compact-component links are not thermal walls.
An internal bulk eddy field requiring a negative correction is incompatible
and refuses rather than inventing negative thermal resistance. Duplicate
walls/contacts, competing thermal owners, through-flow, invalid inputs,
exhausted quadrature and cancellation refuse without mutating inputs.

Equilibrium applicability is also stated in Idaho National Laboratory's
`INSFVTurbulentTemperatureWallFunction` documentation: detachment and
recirculation require a different model, not an uncalibrated correction to
this law. The same limitation applies during an evolving-flow calculation.

## Focused checks

```sh
cargo test -p fs-lbm --lib conjugate::turbulence::thermal_wall
```

Five kernel groups cover laminar and exact analogy limits, independently
integrated values, scaling and refusals. Five internal-coupling groups cover
six orientations on a graded grid, free solid temperature and series heat,
physical contacts, tangent/interior ownership, zero-speed equivalence and
atomic refusals. Five exterior groups add graded temperature/film boundaries,
imposed flux and wall temperature, relative motion and exact cache reuse,
joint internal/exterior transport, and boundary ownership/cancellation.

Native compilation/execution was unavailable in the authoring environment
(no Rust toolchain). Earlier independent Python checks covered 30 kernel and
18 conservative matrix cases. The exterior extension additionally passed 54
independent quadrature/series-resistance/matrix cases with maximum steady
energy imbalance 2.7e-16 W. These mathematical checks do not replace executing
the Rust tests or empirical turbulent heat-transfer benchmarks.
