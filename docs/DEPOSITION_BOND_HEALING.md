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
synthetic opening/shear strength and fracture-energy table. Two independent
fresh-specimen tests then measure opening and shear reactions and fracture work
per actual contact area. It is not a
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

`StrengthMap::new_mixed_mode` accepts `MixedModeStrengthKnot` observations with
separate opening and shear peak traction, and separate mode-I and mode-II
fracture work at each dose. `MixedModeStrengthConfig` supplies the original
normal/shear stiffnesses, BK exponent, compression penalty and separation limit.
No ratio between opening and shear is guessed. The existing mixed-mode
`BilinearCohesive` owner applies its quadratic initiation criterion and freezes
the BK energy mix at first initiation; unloading, irreversible damage, shear
rotation and compression after failure retain that owner's behavior.

Each supplied column must be finite, nonnegative and nondecreasing with dose.
All-zero knots may explicitly represent an induction interval without bonding.
A developed bond needs all four mode observations and must satisfy the original
cohesive admissibility checks, both at supplied knots and after interpolation.
Selection is logarithmic in the bounded knot count; it requires no numerical
quadrature or new constitutive solve. The same `freeze_bond` method carries
either common-mode or mixed-mode laws into the existing mechanical consumers.

`freeze_bond` carries the actual area, cell IDs, dose and mapping into a fresh
mechanical stage. A separated or never-contacted interface refuses. Never use
this handoff to reset existing cohesive damage. Concurrent healing/fracture
and a full coupled thermodynamic balance are not implemented.
The separate `fracture::FractureCooling` stage now consumes real cohesive
endpoints and changes the following thermal interval's contact. It can import
existing damage; reclosure restores only the explicitly supplied thermal
conductance, never a mechanical bond. See the staged coupling below.

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

Eleven material tests and nine coupling tests cover isothermal dose, crossing
knots, cooling/reheating, independent cooling-integral refinement, thermal
primal/audit parity, contact windows and area partition, equal-time recontact,
strength/fracture handoff, atomic refusals and same-grid replay. These are
numerical regressions, not experimental calibration of a filament or adhesive.

## Cohesive failure changes subsequent cooling

```sh
cargo test -p fs-conduction --lib fracture
cargo run -p fs-couple --example debonded_cooling
```

`fracture::FractureCooling` binds existing cohesive points, original thermal
links and explicit `ContactPolicy` values. Intact/partially damaged patches
keep the bonded conductance. After complete failure, a normal opening threshold
removes the link and restores ambient-facing area on both cells. A separate
closure threshold permits frictionless thermal recontact, with an independently
supplied W/K conductance. Damage and dissipated fracture energy remain intact.
Shear failure at a closed normal gap changes conductance without fabricating
lost area. Excessive tangential travel refuses the fixed-area approximation.
Threshold hysteresis is a declared numerical closure model, not a contact theorem.

Advance heat with the accepted closure, stage a mechanical equilibrium using
those temperatures where supported, and call `synchronize` at the exact thermal
endpoint with one `(interface_id, CohesivePoint)` per binding. Commit the staged
mechanical result only after synchronization succeeds. The owner reproduces
constitutive increments against retained history and checks exact law identity;
it does not solve equilibrium again. Lost or inconsistent irreversible history,
wrong IDs/time, stale thermal links and cancellation refuse without
publishing a partial interface update. Scheduled changes of mechanically owned
links also refuse, including same-time remove/reinsert sequences. Unrelated
births and thermal boundary changes continue normally.

The `debonded_cooling` example actually solves the existing `fs-solid` cohesive
operator under a compliant actuator, opens and recloses the interface, and
reports force balance, damage, both cell temperatures and changed exposure.
All inputs are synthetic. The example stages heat, mechanics and contact before
publishing their common endpoint. Its mechanical coefficients stay fixed: no
thermoviscoelastic feedback or calibrated detachment time is claimed. The
contact classification governs the next thermal interval; a missed within-step
opening is not retroactively corrected. Refine physical coupling times to
resolve events. Full geometry motion, new contact search, simultaneous healing,
within-step toughness/temperature feedback and a closed thermodynamic balance
remain outside this stage. Explicit endpoint heat conversion is available below.

Eight additional Rust regressions cover analytical contact/no-contact cooling,
independent reclosure conductance without healing, shear/partial damage,
cohort/history/time/area refusals, scheduled ownership conflicts, cancellation,
area partition and restart.

## Optional heat from new fracture work

Use `FractureCooling::synchronize_with_heat` with one `FractureHeatPartition`
for every bound interface. `fraction_to_heat` declares the fraction of new
cohesive dissipation converted to immediate heat; `fraction_to_first` partitions
that heat between the original two cells. No universal conversion or partition
is inferred. For each original area, the same authoritative cohesive trial gives
`D = area * new_dissipation`; the method deposits `fraction_to_heat * D` and
reports the remainder as nonthermal fracture energy. Do not add either heat or
the upstream mechanical debit again to a combined balance.

The existing `DepositionSimulation::apply_heat_j` owner updates each cell through
its original enthalpy law, including latent plateaus. No thermal time advances,
source watts are not reintegrated, and no elastic storage, old damage, unloading,
reclosure or numerical dissipation becomes heat. Replaying an accepted point
produces zero new heat. The original `synchronize` remains contact-only.

Contact changes, irreversible history and thermal enthalpy form one candidate.
Domain, cancellation or energy refusal publishes none of them. An explicit
joule tolerance covers local partition defects, heat aggregation, represented
cell-storage errors and final work/nonthermal/thermal report rounding. The
returned `transfer_residual_j` equals reported thermal storage plus reported
nonthermal fracture energy minus reported fracture work, without subtracting
large absolute thermal reference energies.

Seven added regressions cover analytic work and unequal masses, latent heat,
replay/unloading/imported damage, zero conversion, shared-cell area partition,
atomic refusals and disparate-magnitude report rounding. This is an endpoint
staggered model with fixed areas and immutable toughness; it does not resolve
crack-tip temperatures or certify a material's heat-conversion fraction.

## Physical bond viscosity and timed heat transfer

A bond created with `BilinearCohesive::with_viscosity` needs
`FractureCooling::synchronize_over_time`, or
`synchronize_over_time_with_heat` when its physical losses enter the thermal
cells. The original `trial_over_time` law supplies both intrinsic fracture work
and integrated surviving-bond dashpot loss along the declared straight segment.
Viscous loss accumulated before failure remains present even when the final
dashpot traction is zero. Elastic unloading may generate viscous heat without
any new fracture work.

Compute the duration as `target_time_s - state.last_synchronized_time_s()` and
use that same value in the mechanical material trial and synchronization.
Heat-transport stops leave this mechanical clock unchanged. The target must
equal the current thermal time; the supplied duration must equal its represented
difference from the previous mechanical time. A repeated timed synchronization
refuses, while a later stationary hold contributes no repeated loss.

One `CohesiveHeatPartition` per interface supplies independent fracture and
viscous heat-conversion fractions, with separate receiving-side splits.
`ContactReport` reports both losses, and `FractureHeatReport` retains both
nonthermal remainders. The transfer residual closes their sum against the
native enthalpy change. Neither channel implicitly consumes the other, and
the caller must avoid double counting the upstream mechanical work.

The complete contact/history/clock/enthalpy candidate publishes atomically.
Final points cannot reveal every unreported trajectory: an undamaged closed
loop, for example, can dissipate viscous work while returning to the same
cohesive point. Supply each physical reversal as a separate segment. These APIs
retain the staggered endpoint model; they do not infer a within-step fracture
time, update toughness from the generated heat, or implement concurrent healing.
