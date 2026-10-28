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
source is introduced.

`solve_energy_with_wall_heat` performs this preparation and calls the original
energy solver. `EnergySolution::solid_to_fluid_heat_w` and global balance consume
the same retained face resistance. Each wall report distinguishes the physical
contact, native fluid half-cell and wall-model correction.

A prepared setup can also be passed to the existing transient energy driver
while that velocity/eddy field is held fixed. When flow changes, prepare again
from the **original** setup; applying preparation to its own output refuses.
No new flow integrator, clock or temperature storage is introduced.

## Scope and refusals

This is a smooth-wall equilibrium closure suitable for explicitly selected
LVEL cooling studies. It does not establish the origin or validity of an
arbitrary supplied flow field. It does not model separation, roughness,
variable fluid properties, viscous heating or resolved-wall LES. Open faces,
domain-boundary walls and compact-component links are not replaced. A bulk
eddy field requiring a negative correction is incompatible and refuses rather
than inventing negative thermal resistance. Duplicate walls/contacts,
competing thermal eddy fields, through-flow, invalid inputs, exhausted
quadrature and cancellation also refuse without mutating inputs.

## Focused checks

```sh
cargo test -p fs-lbm --lib conjugate::turbulence::thermal_wall
```

The five kernel groups cover laminar and exact analogy limits, independently
integrated values, scaling and refusals. Five coupled groups cover all six
wall orientations on a graded grid, free solid temperature and series heat,
physical contacts, tangent/interior ownership, zero-speed equivalence and
atomic refusals. These Rust tests and the native example were added but not
compiled or executed in the authoring session: no Rust toolchain was available.
Independent Python numerical comparisons passed 30 kernel cases and 18
conservative matrix cases; they are not a substitute for the native test run.
