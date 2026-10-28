# Temperature-driven cohesive separation

`fs_solid::cohesive::thermal_joint` solves a bounded, uniaxial series assembly:
one or more thermoviscoelastic members plus a single normal cohesive joint.
Temperatures and total end-to-end extension are supplied; joint opening,
member strains and the common reaction force are outputs. Holding total
extension at zero models restrained cooling, not a prescribed joint gap.

Each `Member` starts stress free with an explicit length, area, immutable
`ThermoMaxwell` law and birth temperature. Different members retain different
birth references and temperature histories. `Joint` retains the material
branch strains and irreversible cohesive point through every interval; a
new temperature-dependent tangent does not create a virgin material state.

## Physical solve

The existing material owner gives each implicit endpoint stress as
`sigma_i = C_i epsilon_i + b_i`. Axial series force balance and compatibility
then give

```text
S = sum_i L_i / (A_i C_i)
free = -sum_i L_i b_i / C_i
F = A_b traction(g; accepted damage)
g / S + A_b traction(g; accepted damage) = (D - free) / S
```

The existing cohesive Newton/FGMRES owner solves that last equation with the
same accepted damage on every nonlinear probe. The resulting force determines
each member's strain; the existing material owner then proposes all internal
state updates. Actual member forces, total extension and independent energy
checks must pass before the complete candidate is returned. Thermal/material
range errors, work limits and cancellation leave the original joint unchanged.
Memory is linear in members and retained Maxwell branches, with explicit
member, branch-state and Krylov-basis limits. Operator and iteration allowances
are supplied through the existing cohesive configuration.

Energy accounts retain bulk thermal-eigenstrain and modulus work, bulk stored
energy, physical and backward-Euler numerical dissipation, and cohesive path
work/storage/fracture. Endpoint actuator quadrature `F_new delta D` is NOT
called exact work. Its signed cohesive correction
`F_new delta g - integrated cohesive work` is reported separately; it can be
negative during softening and is not physical dissipation. The combined
check includes this correction, and independently checks force-work conjugacy
and total stored-energy change. A small defect is not physical validation.

## Coupled example

```sh
cargo test -p fs-solid --lib thermal_joint
cargo run -p fs-couple --example thermally_debonded_joint
```

The example cools a synthetic initially hot restrained member against a finite
thermal body. On each endpoint it advances heat, solves the temperature-driven
joint and synchronizes the accepted cohesive point with `FractureCooling`.
All three candidates are checked before publication. Separation removes the
thermal link and changes subsequent cooling. The internal gap is never a
prescribed input. Fracture energy is reported, not injected into enthalpy.

For explicit endpoint heat conversion, callers can synchronize the same staged
point through `FractureCooling::synchronize_with_heat`, supplying conversion and
two-cell partition fractions. That adapter preserves latent storage and reports
the unconverted fracture energy; it does not update toughness within the step.

This uses endpoint-staggered coupling; refine both the thermal timestep and
the mechanical/coupling interval. Unlike an explicitly prescribed linear gap
history, sparse equilibrium endpoints do not certify intervening contact events.
The example must not be passed to the prescribed-motion API as exact motion.

## Limits and checks

This is an axial coupon/restraint, not a spatial 3D peel, bridge, or growing
mesh. Member temperatures are uniform, areas and frames are fixed, and there
is no automatic contact discovery, bond creation, thermal-material calibration,
dynamic snap-back, bulk failure, simultaneous healing/fracture, fracture heat,
or full entropy model. The nonlinear solver establishes a numerical equilibrium,
not uniqueness or stability; subdivide paths at actual reversals.

Six native regression functions cover analytical restrained cooling, thermal
failure and compression after reheating, Maxwell/series-spring relaxation,
member partition, nonuniform member temperatures, restart and atomic refusals.
They are authored but not run in this environment: DSR/RCH/Cargo are absent.
Independent Python/SciPy references compared 1,280 eliminated and monolithic
member/joint equilibria and their discrete energy identities. Those references
do not execute the Rust implementation or qualify a real material.
