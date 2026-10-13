# Deposition temperature history to bond development

This implements the bounded R3 exposure/strength handoff of issue #1 on the
existing **explicit** growing-domain thermal kernel. It does not complete the
printing umbrella or turn illustrative inputs into calibrated filament data.

## Run

```sh
cargo test -p fs-material --lib healing
cargo test -p fs-conduction --lib bonding
cargo run -p fs-conduction --example deposited_bond_healing
```

Use the repository's required execution/toolchain setup for these commands.
The example deposits a hot second road onto an earlier road, compares two
cooling conditions, and freezes the resulting exposure through an explicit
synthetic strength/fracture table. A new displacement-controlled cohesive pull
then measures reaction and fracture work per actual contact area. It is not a
whole-structure failure calculation or a calibrated strength prediction.

## Physical data and APIs

`fs_material::healing::MobilityTable` accepts nonnegative mobility rate `r(T)`
in 1/s, absolute-temperature knots in K, and a law/version identity. The
history is dimensionless `D = integral r(T) dt`. The table is integrated
analytically on its piecewise-linear rate curve and a declared linear
within-step temperature path, splitting at every crossed knot. This is not
exact integration of an unknown thermal path. Zero-rate cold intervals must
be supplied explicitly; no default glass transition or Arrhenius constants
are inferred. Reheating adds to the same dose. An application can instead
implement the bounded, immutable `HealingLaw` trait.

Equivalent-time welding provides the motivation (Seppala et al., *Soft Matter*
13, 2017, 6761-6769, doi:10.1039/C7SM00950J). It does not justify a universal
strength law. Orientation, disentanglement, crystallization, degradation,
contact pressure and evolving intimate-contact geometry remain outside this
scalar model. The table and its applicable temperature range need independent
material data.

`fs_conduction::transient::enthalpy::deposition::bonding::BondedDeposition`
binds declared healing laws to individual, ordered thermal interface IDs.
Unbound thermal contacts remain thermal-only. Each contact cohort retains its
original area, actual contacting duration, dose and separation time. The
interface temperature is an explicitly chosen weighted mean of its two cells;
this assumption must be justified or refined by the application.

The existing thermal kernel supplies **every actual CFL-limited substep** and
each individual equal-time event in original order through `DepositionObserver`.
No final-temperature-only reconstruction, output resampling or second thermal
solver is used. A birth is not a contact. Separation stops exposure while
preserving the old history; recontact under an old bond ID refuses, including
separation/recontact at the same timestamp. New area needs a new cohort ID,
so newly contacting material cannot inherit the age of an older larger patch.

`StrengthMap` is a separate, explicit dose-to-traction/fracture-energy relation
bound to the mobility identity. It linearly interpolates supplied Pa and J/m2
values only inside its dose domain. Zero developed strength returns `None`,
and missing maps remain missing. This API provides no default saturation or
calibration authority. Positive outputs are admitted by the existing
`BilinearCohesive` constructor with explicitly supplied Pa/m stiffness,
compression stiffness and metre separation limit.

`freeze_bond` carries the actual area, cell IDs, dose and mapping into a fresh
mechanical stage. A separated or never-contacted interface refuses. Never use
this handoff to reset existing cohesive damage. Concurrent healing/fracture,
thermal conductance updates from damage, fracture-heat feedback and a coupled
thermodynamic balance are not implemented.

## Execution and numerical scope

Thermal state, event cursor and bond histories are cloned and published only
as one successful candidate after all thermal and bond checks. Domain, energy,
work-budget and cancellation failures publish neither partial heat nor dose.
Limits bound retained interfaces and per-advance history updates; the latter
reset per call but physical exposure does not. History loops poll cancellation
every 256 records. As with the original thermal kernel, bounded cloning and
schedule construction are not internally cancellable. Generic observers must
deep-copy mutable state and avoid external side effects.

Clone/restart on the same thermal step grid preserves numerical history.
Inserting observation stops may change that grid; compare time refinement,
spatial refinement and interface-temperature modelling separately. Mobility
integration exactness is not a temperature or physical-strength error bound.

Six material tests and eight coupling tests cover isothermal dose, crossing
knots, cooling/reheating, independent cooling-integral refinement, thermal
primal/audit parity, contact windows and area partition, equal-time recontact,
strength/fracture handoff, atomic refusals and same-grid replay. Native Rust
execution is unverified in this editing environment. Independent Python/SciPy
checks accompany the patch but do not substitute for compiling these tests.
