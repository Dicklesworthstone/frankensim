# fs-nuclear contract

## Prescribed fission, finite inventory and delayed energy

`fs-nuclear` is an L3 source owner. The source kernel uses the shared deterministic
`fs-math` decay primitives and shares immutable nuclear model data with `Arc`.
It accepts a prescribed total
fission rate in expected nuclei per second, finite fuel inventory, explicit
initial delayed-energy reservoirs, SI energy coefficients and caller budgets.
An identified dataset and its domain of applicability remain the caller's
responsibility. There is no built-in isotope or silently assumed irradiation
history. The supplied data identity is a label, not scientific validation.

`EnergyChannels` separates energy deposited in the source, outward radiation
available to another transport owner, and nonrecoverable export. Each model
has a prompt partition and an ordered list of nonnegative exponential decay
groups. A group's coefficients are integrated joules per fission, with a
positive decay constant in inverse seconds. They describe effective energy
reservoirs, not delayed-neutron precursor populations. For impulse-power fit
coefficients `alpha * exp(-lambda * t)`, integrated energy is `alpha/lambda`.

`DecayHeatGroup::from_half_life_s` admits a finite positive effective half-life
and derives `lambda = ln(2)/half_life`; `half_life_s` returns the corresponding
finite half-life of a directly authored rate. These group lifetimes govern
stored energy release, not a claim about any individual isotope's population.

For a constant prescribed rate `F`, each group obeys

```text
dE_i/dt = q_i F - lambda_i E_i.
```

The accepted event count is limited by actual remaining fuel. If it is exhausted
within an interval, fission stops at `N/F` and all groups continue decaying to
the interval endpoint. For active duration `a`, cooling duration `c`, and
`phi(x) = (1-exp(-x))/x`, the endpoint is

```text
E_i,new = E_i,old exp(-lambda_i (a+c))
        + q_i (F a) phi(lambda_i a) exp(-lambda_i c).
```

Integrated release is computed with `expm1` and a short-interval series,
independently of subtracting large stored-energy totals. It therefore retains
small positive decay heat. No explicit decay substeps or endpoint-power
rectangle are needed. Group storage remains nonnegative. Turning the rate to
zero does not erase the accepted reservoirs or stop their heat release.

Scaled survival uses `fs_math::decay::attenuate`, retaining a representable
large-reservoir tail even when the unscaled exponential underflows. Overflow
of `lambda * duration` uses the fully decayed limit for existing energy;
uniformly generated surviving energy evaluates the equivalent divided product
without first rounding its retention fraction to zero. Unrepresentable positive
exposures and nonfinite physical outputs still refuse through the original
transaction and accounting budgets.

## Accounting, numerics and work

Every accepted interval reports fissions, active duration, generated energy,
three released channels, surviving delayed-storage change, a fuel residual,
an energy residual and accumulated absolute numerical energy defects. The
energy identity is

```text
Q * accepted_fissions = delta_delayed_storage
                     + deposited + outward_radiation + nonrecoverable_export.
```

Generated energy is independently anchored to the immutable total yield Q and
the same event count that debits fuel. Rounding in the sum of individual
component products is charged against that independent total. Local group
defects cannot cancel one another to satisfy the energy budget. Compensated accumulators retain small repeated fuel debits and
lifetime release credits. An unresolved reservoir debit consumes the explicit
energy allowance; it cannot silently create heat. Nonfinite inputs, negative
rates/inventories/channels, unrepresentable physical quantities, excess work
and excess accounting defects refuse the candidate. The tolerances bound
numerical account closure, not nuclear-data uncertainty.

Each trial uses O(groups) work and storage, polls the caller checkpoint before
work, at every group, and before publication. `trial` never modifies its input;
`advance` publishes only a fully accepted candidate. The group count is bounded
by the caller. Cloning retains the full accepted inventory and accumulated
accounts. Replay with identical coefficients and interval partition is bitwise
stable on the same ISA/toolchain profile. Different interval partitions
agree to floating-point tolerance; cross-ISA equality is not claimed.

## Physical scope

The fission source is a reduced prescribed-source model. It does not infer neutron flux,
cross sections, criticality, reaction rates, temperature reactivity feedback,
isotopic product yields, activation, poisoning or a nuclide depletion network.
It does not solve radiation transport or assign dose. The radiation channel is
an energy export, not proof of where that radiation subsequently deposits.
Group fits require independent calibration over the intended cooling times
and irradiation histories. No regulatory or reactor-design qualification is
claimed.

Fission cannot be represented by an ordinary chemical reaction preserving
each element. These source accounts therefore do not mutate `fs-thermochem`
species or invent a chemical stoichiometry. Likewise, fuel nuclei consumed by
fission become products; they are not lost bulk material mass. Thermal material
mass, composition evolution and mass-energy changes require their own physical
owners. A fixed material caloric law is a declared reduced-model assumption.

## Evidence and references

The initial 12 owning source tests cover G0 account closure/admission, G1 independent
constant-source and shutdown solutions, separated channels, multiple decay
timescales, short and stiff intervals, finite fuel exhaustion, G3 interval
partitioning, G4 cancellation/refusal rollback, and G5 identical checkpoint
replay. These passed using the pinned Rust 1.99.0-nightly toolchain through
`rustc --edition=2024 --test crates/fs-nuclear/src/lib.rs`. Coefficients are
explicitly synthetic: equation verification is not isotope validation or a
performance benchmark. Workspace Cargo/DSR verification remains separate.

Six additional source regressions cover half-life admission and conversion,
scaled normal/subnormal late tails, overflowing exposure, representable forced
source retention, and accounting refusal. All 19 current source tests passed
through the pinned native `rustc` test build linked to the real `fs-math` crate.
This is a focused source-kernel result, not a workspace or thermal-feature gate.

The exponential source-history formulation and distinction between prompt and
delayed fission energy follow the public [Argonne SAS4A/SASSYS-1 decay-heat
documentation](https://sas-doc.nse.anl.gov/5.7/Part02/Ch04/decay-heat.html),
especially sections 4.4.1.1 and 4.4.1.3. The implementation does not copy or
supply the ANS coefficient tables. [INL BISON's power and burnup
documentation](https://mooseframework.inl.gov/bison/theory/power_burnup.html)
provides context for prescribed fission sources and why energy production,
fuel history and material evolution are distinct quantities.

## Shared material thermal integration

The opt-in `thermal` feature provides `thermal::FissionHeating` on the existing
`fs_thermal::enthalpy::DepositionSimulation`. Bind each nuclear inventory to
an already active `CellId` and its original `MaterialId`. Original material
calorics, reference mass, geometry, storage basis, contact links and boundary
reservoirs keep their physical owners. Where the thermal law supplies a full
caloric identity it is also checked; an absent identity remains unavailable.
Burnup does not silently change a material chart, thermal mass or phase law.

A validated right-continuous `FissionRateSchedule` supplies total fissions per
second. Binary-search transition lookup clips the original thermal substep at
every authored rate change. Source trials run through `step_sources` only after
CFL, event, transition and target clipping have established the final duration.
Integrated deposited joules divided by that same duration enter the original
heat-source path once, in addition to independently authored boundary heat.
The emitted-radiation and nonrecoverable channels are not credited as heat.
`requires_step_interval` prevents an instantaneous balance from silently
omitting this duration-dependent source.

The original thermal step must succeed before its matching nuclear state is
accepted. The existing deep-cloned observer transaction discards both histories
on material-domain, budget, numerical or checkpoint refusal. The exact original
thermal checkpoint error is preserved through source cancellation. Per-advance
compensated accounts bound accumulated absolute core and thermal-credit defects;
ordinary thermal storage/transfer admission remains independently mandatory.
The direct `advance_to` wrapper also validates an idle call's clock and bindings.

The observer composes with other existing physical observers through
`fs_thermal::observers::CoupledObservers`. It supports the original enthalpy
and explicitly selected internal-energy storage bases; any sensible/latent
temperature response comes from the original `ThermalLaw`. Heat-transfer
accuracy remains the native explicit, fixed-reference-geometry, uniform-cell
approximation with independent maximum-step refinement. Source integration is
analytic; that does not make the spatial thermal evolution exact or implicit.

Six owning integration regressions cover analytic shutdown calorimetry with
additive boundary heat, finite-coolant contact and first-order refinement,
multiple material/source schedules, identical-partition restart, fuel
exhaustion, and joint rollback/admission failures. They were written and
reviewed against the native APIs, but could not be compiled or run in the
implementation session: the shared filesystems filled and the execution
environment then became unavailable. Run the real native gate with

```bash
cargo test -p fs-nuclear --features thermal --locked
```

before treating those coupled regressions as executed evidence. DSR and the
complete constellation/layout checks were not available as successful gates.

## Radioactive parent/daughter inventories

`decay::DecayModel` binds caller-owned unique nuclide identities, a dataset
identity, admitted `fs_thermochem::decay::DecayLawV1` laws, and ordered decay
branches. Each unstable parent selects one represented daughter, with explicit
stable endpoints. A branch retains its own probability and deposited, outward
radiation and escaped joules per event. Multiple branches can share a daughter
without merging their energy destinations. Fractions are never silently
normalized and no isotope catalog or equilibrium history is inferred.

`DecayInventory` owns nonnegative expected populations, compensated lifetime
transformation counts and released-energy channels. `activity_bq` evaluates
lambda times the current population; `decay_power_w` reports instantaneous
channel power. Neither diagnostic is used as an interval-energy rectangle.
For every accepted interval, branch events equal the independently integrated
parent transformations times the branch probability. Those events release the
branch's authored energy, including daughters produced and decayed within the
same interval. Stable species have exactly zero activity and transformation count.

Numerical propagation remains owned by `fs_la::decay::DecayPlan`: bounded
nonnegative uniformization, exact survival/loss diagonals at every composition
scale, equal-rate support, arbitrary original node order, work limits and
cooperative cancellation. `DecayModel::prepare` retains both endpoint and
integrated-count maps. Reusing a `PreparedDecay` avoids repeated matrix
preparation; application costs O(nodes squared + branches). The same immutable
model Arc and numerical limits must be retained. Source-bound maps cannot be
silently applied to inventories carrying different half-lives or branch energies.

`trial` and `trial_prepared` return independent candidates; `advance` and
`advance_prepared` publish only after propagation, physical accounts and the
final checkpoint all succeed. Signed total population residuals are supplemented
by a noncancelling sum of per-nuclide transfer, branch-product and lifetime-counter
defects. Tiny unresolved parent debits consume that explicit population budget.
Emission-product and accumulated-energy rounding consume a separate joule
budget. These numerical accounting limits are not nuclear-data uncertainties or
certified forward-error enclosures.

The model is closed, constant-rate, acyclic, and limited to the numerical owner's
64-node bound. It has no irradiation, injection, removal, activation, fission
product yields, temperature-dependent rates or transport of emitted particles.
The represented daughter's population is not a full nuclear mass/charge ledger.
Original chemical composition, bulk thermal mass and material caloric laws remain
with their existing owners. Large inventories that could amplify vanished or
subnormal propagator coefficients retain the numerical owner's explicit refusal;
the physical wrapper never bypasses its range checks or invents a zero tail.

Ten focused physical regressions cover analytic repeated-rate ingrowth and
integrated emissions, distinct energy channels for duplicate daughters, activity
versus interval energy, preparation reuse and model binding, partition/restart,
tiny heat with charged unresolved depletion, stable and zero-duration identities,
input admission, early/mid/final cancellation, and energy/range refusals. All ten
passed, alongside the 19 fission-source tests, in a native Rust source build
using the unchanged production decay modules from `fs-la` and `fs-thermochem`
and the real `fs-math` and `fs-qty` crates. This isolated the applicable owners
from unavailable constellation dependencies; it is not a workspace Cargo gate.

## Shared coupling facade and original-water example

Enable `fs-couple/nuclear-thermal` to use `fs_couple::fission`. This facade
re-exports the same `fs-nuclear` source and thermal types; it introduces no
second inventory, scheduler, heat solver or source implementation. Its dependency
is optional, and it composes on the original thermal observer protocol.

The `prescribed_fission_water` example uses the existing
`fs_material::phase::water::IsobaricWater` IF97 caloric law and native finite-cell
contact. Water is a thermal receiver; its mass is not interpreted as nuclear
fuel. Source coefficients and the prescribed history are explicitly synthetic.
The isobaric caloric chart retains sensible and latent heat, while reference
geometry and supplied contact conductance retain their original reduced-model
assumptions. The example runs a source/shutdown/cooling history and checks
continuing decay heat, fixed post-shutdown fission count, coolant warming and
total deposited-energy/thermal-storage closure.

```bash
cargo run -p fs-nuclear --features thermal --example prescribed_fission_water --locked
cargo check -p fs-couple --features nuclear-thermal --locked
```

The example, facade check, and a thirteenth core regression for independently
budgeting total-yield product roundoff were added after the execution environment
became unavailable. They have not been compiled or executed in this session;
the final snapshot must be rerun through the native checks above. The previously
executed 12-test result is not a claim that these later changes passed.

