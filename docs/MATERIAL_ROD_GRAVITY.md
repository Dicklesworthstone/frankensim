# Material weight, thermal softening and rod equilibrium

`SpatialMaterialRod::with_gravity` connects the existing source-resolved
concentric sections to actual body loading in the existing Cosserat solver.
Gravity is an explicit world-frame acceleration vector in m/s2. There is no
implicit Earth gravity and no material-name branch.

## One reference mass, two consumers

`reference_mass()` retains each original section and its density receipts,
reference segment length, segment masses, nodal shape-function masses and total
mass. Each uniform-reference segment contributes half its mass to each endpoint,
including the clamped root. This integrates the linear centreline shape functions
and preserves the weight and its first moment; it is not a dynamic mass matrix.

`RodMaterialSource::resolve_at` shares this immutable reference inventory with
every temperature snapshot. Current rho(T) remains in the current source bundle,
but it does not redefine mass on frozen geometry. The heat-driven example uses
those same segment masses for enthalpy storage. Heating can soften E/G and change
expansion, deformation and support forces without silently changing weight.

```rust
let mut weighted = source.resolve_at(&layer_temperatures_k, &gate)?
    .with_gravity([0.0, 0.0, -1.0], &gate)?;
weighted.solve_static(&tip_load, 1, 1e-8, &gate)?;
let stored_elastic_j = weighted.energy_j()?;
let gravity_potential_j = weighted.gravity_potential_j()?;
let root = weighted.root_support_reaction(&gate)?;
let other_supports = weighted.support_reactions(&tip_load, &gate)?;
```

The gravity potential is referenced to the undeformed straight pose and enters
all existing Newton probes, line searches and reaction recovery. Elastic storage
and gravitational potential remain separate observations. `with_gravity` replaces
rather than adds another weight field; a zero vector removes it. It preserves
supports and the previous pose as a static initial guess. During `solve_static`,
gravity is permanent and the existing load stages ramp only the separate tip
load. This is not work-conserving dynamic activation of gravity or supports.

Root force and body-frame moment come from centred energy derivatives, including
root weight, not from forcing a global balance. Additional supports expose their
constrained reactions and the free-coordinate residual. Check that residual when
interpreting reactions as equilibrium. Probe/domain/cancellation failures retain
the original specimen.

## Heat-driven cantilever executable

The existing `heated_material_rod` example now supports `cantilever` (root clamp
only) and an optional explicit gravity vector. The original `free` and `fixed`
choices still mean straight axial guides with respectively free or restrained
axial motion at the far end. Omitting gravity still means no weight, including
for the new cantilever mode.

These inputs differ only in synthetic material heat capacity:

```sh
cargo run --release -p fs-couple --example heated_material_rod -- --demo 500 cantilever 1 4 20 --gravity 0 0 -1
cargo run --release -p fs-couple --example heated_material_rod -- --demo 1000 cantilever 1 4 20 --gravity 0 0 -1
```

Use `--pack MATERIAL.fsmat` instead of `--demo CP` to read the existing canonical
material-card format with the same source requirements as `HEAT_DRIVEN_ROD.md`.
The same heater, thermal network, source reevaluation and mechanical solve are
used in every support mode. Heat/time/pose publication remains transactional.

Cantilever or explicit-gravity invocations append transverse tip coordinates,
reference mass, gravitational potential, root reactions and free residual to
CSV. Original invocations retain their old CSV columns. The fixed thermal
geometry, film coefficient and source conditions are printed explicitly.

Thermal-rod translation derivative probes now scale with reference segment
length, while rotational probes remain angles. Previously the nested 0.0003 m
and 0.000003 m probes could exceed the example's unchanged 0.03 axial-strain
limit on 0.01 m segments before any Newton step. The same scaling now applies
in tangent, residual and root/additional reaction evaluations. The raw legacy
`Rod` and unheated elastic material paths retain their original probe scale
and exact uniform/spatial parity. This does not guarantee
admission at a material boundary or establish derivative accuracy for every
geometry; the physical strain limits are not relaxed.

## Verification and domain

Seven library and three example regressions are authored for mass/first moments,
axial and transverse equilibrium, density-only substitution, mass-preserving
thermal softening, reactions, zero-gravity behavior, short-element probes and
transactional failures. The original rod inline tests and four existing heating
example tests are unchanged.

```sh
rch exec -- cargo test -p fs-solid --test rod_gravity --test rod_material_source
rch exec -- cargo test -p fs-couple --example heated_material_rod
```

Native Rust/RCH/DSR executables were absent in the implementation environment;
compilation and these Rust tests remain unverified. Independent NumPy/SciPy
three-dimensional rod equations reproduce density-only sag and axial loading.
A separate four-cell heat/rod reference reproduces the old short-element refusal
and passes with scaled probes under the same strain limit. At two seconds,
Cp=500 versus 1000 J/(kg K) gives 0.340815 versus 0.331445 mm downward tip motion,
from 0.322481 mm initially, at unchanged 0.000502655 N weight. These are synthetic
reference-equation outputs, not native executable results or measured materials.

The model remains bounded quasistatic, linear-section Cosserat mechanics with
prescribed or one-way solved thermal states. Heat geometry is fixed; gravitational
and elastic work are not returned as heat. No dynamics, buoyancy, radial expansion,
yielding, creep, helical winding, buckling/stability qualification, moving thermal
mesh or total thermomechanical conservation is claimed. Refine thermal and coupling
steps as well as spatial resolution before interpreting a trajectory quantitatively.
