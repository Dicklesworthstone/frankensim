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

The base `FissionState` is a prescribed-source model. The irradiation layer
below additionally computes reaction rates and parent removal from supplied
multigroup flux and cross sections. Neither fission lane infers neutron flux,
evaluates microscopic nuclear data, solves criticality, or evolves irradiation
product yields, activation or poisoning. Radioactive daughter histories have
the separate decay owner described below.
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
exhaustion, and joint rollback/admission failures. All six now pass in the
focused native source build described below. That build exposed and fixed a
method-resolution error in the original fission schedule hook: its iterator
accessor is now explicitly qualified as `FissionHeating::sources`.
The complete Cargo gate remains

```bash
cargo test -p fs-nuclear --features thermal --locked
```

DSR and the complete constellation/layout checks were not available as
successful gates.

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

## Radioactive decay heat on the shared thermal timeline

With the existing `thermal` feature, `decay::thermal::DecayHeating` binds each
explicit inventory to an already active thermal cell and its original material.
It uses the original `DepositionObserver` transaction and the final accepted
interval after thermal CFL, event and target clipping. Only integrated deposited
joules divided by that same duration enter the native heat-source path. Authored
boundary heat remains additive; outward radiation and escaped energy are never
credited as local heat.

The original thermal step must succeed before the corresponding populations,
transformation counts and energy history are accepted. Refused material-domain,
work, numerical-account or checkpoint candidates discard both histories. The
caller's exact thermal checkpoint error is preserved. Per-advance accounts sum
absolute population and energy defects across all cells and substeps, so an
unresolved parent debit cannot repeatedly spend a fresh tolerance each substep.
The original thermal energy budget remains independently mandatory.

Each binding caches its immutable prepared decay map by exact interval duration;
repeated durations use `trial_prepared`, and clones share maps through `Arc`.
Changing the interval prepares a matching map before any publication. Inventories
retain their own history; sharing a model or map never shares mutable populations.
Original thermal material identity, caloric identity when available, storage basis,
reference mass and geometry are checked before and after the common step.
Decaying nuclei do not silently change the material caloric chart or bulk mass.

`fs_couple::radioactive`, enabled by the existing `nuclear-thermal` feature,
re-exports these original owners. `CoupledObservers` can combine radioactive
inventories and the existing prescribed fission source in one thermal advance.
The facade adds no second source implementation or thermal solver. Source
integration remains analytic up to the numerical propagator's admission; thermal
contact remains the original explicit uniform-cell approximation and needs its
own timestep refinement.

Eight native coupling regressions pass: separated-channel half-life calorimetry,
daughter ingrowth heat within one interval, contact convergence, multiple-source
restart and internal-energy storage, common rollback and checkpoint preservation,
accumulated unresolved-debit refusal, binding/instantaneous-balance admission,
and composition with the existing fission observer. The passing run compiled the
actual `fs-nuclear` root with `thermal` enabled: 43 tests total (29 source, six
fission thermal, eight radioactive thermal). Its thermal dependency selected the
unchanged current production enthalpy, observer, coverage, boundary and composition
code while excluding unrelated module declarations; no replacement heat law or
integrator was used. This focused source result does not claim a full workspace
or `fs-couple` Cargo build, which still needs the missing constellation checkouts.

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

The original-water example and full facade Cargo check remain unexecuted; those
checks still need the complete constellation. The thirteenth fission regression
for independently budgeting total-yield product roundoff is included in the
passing source runs described above.

## Runnable parent/daughter example

`examples/decay_chain.rs` constructs a synthetic two-second parent, a half-second
daughter and an explicit stable endpoint. It starts with 1e12 expected parent
nuclei and reuses one prepared half-second map for 16 advances. CSV output reports
all populations, activities, separate lifetime energy channels and daughter heat.
Assertions check population conservation, branch counts, four parent half-lives
and lifetime-energy consistency. It requires no thermal feature:

```bash
cargo run -p fs-nuclear --example decay_chain --locked
```

The exact public-source compile and example execution passed. At eight seconds,
the parent population is 6.25e10 and cumulative deposited energy is
0.4625015258789 J, including 0.2750015258789 J from daughter transformations.
All coefficients are explicitly synthetic; this verifies integration and accounts,
not an empirical isotope model. This direct source run is separate from the
unavailable full Cargo workspace gate.


## Multigroup irradiation and competing parent removal

`irradiation` admits an explicit `NuclideId`, immutable identified energy-group
edges in joules, microscopic fission/capture cross sections in square metres,
an absolute group-integrated scalar flux in neutrons per square metre per second,
a fixed coefficient applicability interval in kelvin, the existing fission
energy model, and an explicit prompt capture-energy partition. The flux is not
a normalized probability spectrum or a surface current. Grid identity and all
edges must match, including for zero flux; no silent regrouping or unit conversion
occurs. Data labels and temperature limits are declarations, not evaluated-data
validation. This original fixed-coefficient constructor does not process data.
The shared temperature tables and explicitly flat-in-energy processed-ENDF
collapse described below provide opt-in data paths; self-shielding and spatial
transport remain with the supplying physics owner.

For each frozen incident spectrum the kernel computes `kf = sum(sigma_f,g phi_g)`
and `kc = sum(sigma_c,g phi_g)`, then solves `N' = -(kf+kc)N` analytically. Fission
and capture counts are the corresponding shares of the same integrated parent
removal. A capture is never counted as a fission or credited the fission Q value.
Its supplied immediate energy has its own accounting entry. Capture products and
fission fragments leave the tracked parent inventory but remain bulk material;
this lane neither invents a daughter chain nor subtracts thermal mass.

The original delayed-energy groups are driven by `F(t)=kf N0 exp(-(kf+kc)t)`.
A stable divided difference integrates that exponential source against each
energy group's decay. Coincident depletion/decay constants use the finite
Bateman limit; short intervals use separately computed heat release, and stiff
intervals preserve representable surviving tails. Neutron shutdown stops new
reaction events and retains the accepted delayed-energy reservoirs. The report
separates fission energy, capture energy, released channels and surviving storage,
with independent total-yield accounting and noncancelling numerical defects.

`IrradiationState::trial` shares the original immutable fission model and common
compensated accounts. Work and storage are linear in the admitted neutron and
decay-group counts. It polls the caller between groups and publishes no partial
state on domain, work, arithmetic, accuracy or cancellation refusal. Instantaneous
reaction rates and channel-separated power are available from the same owner.
This exact frozen-source integration does not bound the error of a changing flux,
a collapsed spectrum, a cross-section dataset or a subsequently coupled material.

The equation form follows [OpenMC's depletion formulation](https://docs.openmc.org/en/latest/methods/depletion.html);
source-history convolution follows the [Argonne decay-heat formulation](https://sas-doc.nse.anl.gov/5.8/Part02/Ch04/decay-heat.html).
The implemented lane retains only independent parent removal and effective
fission-energy reservoirs from that broader physics. No evaluated isotope
coefficients, criticality prediction, transport solution or daughter evolution
are bundled.

Twelve irradiation regressions passed in the follow-up native source build,
including independent reaction/capture balances, equal-rate convolution, spectral
identity and temperature refusals, shutdown, replay, sub-ULP depletion, surviving
stiff tails and resolved slow heat after extreme exposure. The build retained
the actual production `fs-la::decay` and `fs-thermochem::decay` modules and rebuilt
the real `fs-math`, `fs-qty` and `fs-blake3` crates. It is focused original-owner
evidence, not a full workspace Cargo, thermal-feature or DSR gate.

## Point-neutron kinetics and explicit temperature feedback

`kinetics` keeps effective delayed-neutron precursor coordinates separate from
both isotope populations and delayed-energy reservoirs. `PointKineticsModel`
requires an identified effective generation time, delayed fractions and precursor
decay constants. `PointKineticsState` requires the actual initial neutron amplitude
and every precursor coordinate in one caller-defined normalization. Calling the
explicit equilibrium constructor is an authored history assumption; state
construction never guesses equilibrium from temperature or power.

For constant reactivity `rho` and external source `S`, the model solves

```text
n'   = (rho - beta) n / Lambda + sum(lambda_i C_i) + S
C_i' = beta_i n / Lambda - lambda_i C_i.
```

A positive shifted Taylor matrix exponential with scaling and squaring advances
the state. Augmented coordinates integrate the constant external source and
`integral(n dt)` in the same solve; no endpoint-rate rectangle supplies the
returned population-time integral. The dimension is the delayed-group count
plus three. Storage is O(groups squared), and arithmetic is
O(groups cubed times (Taylor terms + squarings)). Prompt stiffness increases
scaling logarithmically rather than requiring explicit prompt-lifetime time
steps. This is an algorithmic work description, not a performance benchmark.

Explicit budgets bound group count, Taylor terms, squarings, estimated population
error, estimated population-time error and the population balance residual.
Nonnegative matrix arithmetic does not clip negative physical populations after
the fact. Finite supercritical growth is integrated; unrepresentable growth or
insufficient work/accuracy budget refuses the whole candidate. Checkpoints run
through matrix work and before publication. Clone retains complete precursor
history and permits deterministic same-partition replay on the same platform.
The reported truncation/roundoff estimates are not interval certificates or
nuclear-data uncertainty bounds.

`TemperatureFeedback` admits identified linear and logarithmic reactivity fits,
reference temperatures and explicit Kelvin domains. Its evaluator returns a
reactivity shift only in-domain. It does not determine a coefficient from an
isotope name, solve criticality, choose a flux normalization, or close a thermal
feedback iteration by itself. Neutron amplitude and its integral are not fission
counts without a separately admitted rate or flux mapping. Frozen coefficients
and the fixed spatial/spectral shape are declared reduced-model assumptions.

Eleven kinetics regressions passed with the merged production source: independent
prompt/source exponentials and one-group eigenmodes, critical equilibrium, stiff
prompt relaxation, source buildup and shutdown, positive tiny-time integrals,
partitioned integral agreement, replay, cancellation/work/accuracy/overflow
rollback, and explicit temperature-feedback behavior. Together with the
prescribed source, radioactive inventories and irradiation kernels, the merged
focused source run passed 52 tests. The dependency build retained the complete
original decay modules and real math/quantity/hash crates described above; it
is not the full Cargo workspace or thermal-feature gate.

Equation references: [Argonne reactor power](https://sas-doc.nse.anl.gov/5.8/Part02/Ch04/power.html)
and [delayed-neutron precursors](https://sas-doc.nse.anl.gov/5.8/Part02/Ch04/delayed-neutrons.html).
No evaluated neutron coefficients or reactor-design qualification are supplied.

## Shared thermal irradiation of isotope mixtures

`thermal::NuclearHeating<S>` now owns the common material binding, absolute clock,
authored-transition clipping, pending interval, additive heat credit and atomic
publication. `FissionHeating` remains its prescribed-source specialization;
`IrradiationHeating` specializes that same owner for neutron irradiation. The
sealed `NuclearHeatSource` seam does not duplicate a thermal integrator, nuclear
inventory or energy-account implementation. The optional `fs-couple` facade
exposes these types plus the original `irradiation` and `kinetics` namespaces.

An `IrradiationCellSource` binds a bounded ordered mixture of distinct parent
nuclides to one original active thermal cell. Its immutable
`NeutronFluxSchedule` uses one exact group identity/edge grid and right-continuous
absolute-time changes. Each isotope sees the same incident spectrum and the
original cell temperature. Each advances its own fission/capture inventory and
delayed-energy history; emitted channels and reaction counts are then aggregated.
Fission and capture products remain in the thermal material's mass coordinate.
No composition-dependent caloric chart or daughter yield is silently assigned.

Both initial and final thermal temperatures must remain within every bound
nuclear dataset's declared range. A final-temperature refusal occurs before
pending inventories are accepted and discards the whole thermal advance. Flux
transitions clip the actual native timestep. Only integrated deposited fission,
capture and delayed heat divided by that final duration enters the native source
path, in addition to authored boundary heat. Outward radiation and escaped energy
remain exports. The same per-advance numerical allowance includes absolute
source and heat-credit defects across cells/substeps, including final rounding
when multiple isotope channels become one cell heat rate. A one-joule isotope
contribution cannot vanish beside a 1e16-joule contribution outside that budget.

Five additional thermal regressions passed: independent two-group/two-isotope
reaction calorimetry with additive boundary heat and shutdown clipping; decay-heat
restart at flux shutdown; joint refusal when the thermal endpoint leaves the
nuclear-data domain; duplicate-isotope/group-identity/work admission; and mixed-
magnitude heat-reduction accounting. Together with all source kernels, the six
original prescribed thermal tests and eight preserved radioactive thermal tests,
the focused production-owner build passed **71 tests** with `thermal` enabled.
It used the unchanged native enthalpy, heat-rate, timestep, observer, coverage,
boundary and species implementations. No heat law, integrator, physical type or
nuclear algorithm was substituted. Unrelated thermal wrappers were excluded from
that focused dependency build; full workspace and `fs-couple` Cargo gates remain
separate. The bounded native `cargo check -p fs-thermal --locked --offline --lib`
attempt timed out after 180 seconds in real dependency compilation, without
producing thermal metadata or a source diagnostic.

## Material-temperature feedback through neutron kinetics

`thermal::KineticHeating` specializes the same `NuclearHeating` transaction for
one explicitly normalized fixed-shape neutron region per source cell. Its
`KineticCellSource` retains the original isotope mixture, point-neutron and
precursor history, immutable control schedule, and explicit temperature-feedback
fits. The irradiation schedule is now the reference absolute multigroup flux at
neutron amplitude one; actual exposure scales with the admitted amplitude.
Initial admission checks the actual reference-times-population spectrum as well
as the reference data, so an unrepresentable physical flux cannot enter through
an otherwise finite initial neutron amplitude.

On each final native thermal interval the source reads the original material
temperature, adds the admitted feedback shifts to scheduled control reactivity,
and advances point kinetics. Its independently integrated neutron-population time
divided by that same duration scales the reference spectrum. Existing irradiation
owners then advance actual fission/capture counts, finite parents, delayed-energy
reservoirs and disjoint release channels. The original thermal owner receives
only their integrated deposited energy. The neutron state and all isotope
histories are accepted only after the original material step succeeds.

`KineticsControlSchedule` supplies right-continuous reactivity and external-source
changes; its transitions and reference-flux changes both clip the native timeline.
Explicit feedback, neutron-group, precursor and isotope work limits remain in
force. Both endpoints must satisfy feedback-fit and nuclear-data temperature
domains. Nested checkpoint failures preserve the original caller error. The
report retains both physical interval reports, applied reactivity, mean neutron
amplitude and integral-to-mean conversion roundoff. The latter consumes the
kinetics integral allowance rather than silently becoming extra exposure.

For fixed microscopic data and spectral shape within an interval, total neutron
fluence determines independent parent removal and the competing reaction counts.
Using the mean flux approximates the timing of new delayed-energy production;
material-temperature feedback is also frozen at the beginning of each thermal
interval. Therefore this is an explicitly partitioned coupled model. Refine
`StepConfig::max_dt_s`; the source does not claim an exact or unconditionally
accurate coupled trajectory. Numerical population estimates, nuclear/thermal
account closure, data uncertainty and this partitioning error remain distinct.
The parent-only specialization does not infer burnup-dependent reactivity.
There is no intercell neutron transport, neutron exchange, automatic criticality
calculation or inferred nuclear dataset.

Five new coupled regressions passed: independent exponential population-time
calorimetry; feedback refinement against an independently derived neutron/fluence
invariant; control/reference-spectrum transitions with restart and retained
neutron precursors and decay heat; endpoint and nested-checkpoint rollback; and
control/work/idle/initial-flux admission. The final focused production-owner run
passed **76 tests** with `thermal` enabled: 19 prescribed-source, 10 radioactive
inventory, 12 irradiation, 11 point-kinetics, 11 prescribed/irradiation thermal,
eight preserved radioactive thermal and five feedback thermal tests. It compiled
the complete nuclear crate and rebuilt the real math/quantity/hash owners plus
the complete original decay modules and verified native thermal core/type slice.
No physical type, heat law, numerical algorithm or integrator was mocked. Full
workspace Cargo and `fs-couple` build validation remain separate; the existing
native reproduction gates above still apply when that dependency environment
is available.

## Prescribed production during radioactive decay

`decay::production` integrates an explicitly authored vector of nuclide
production rates while the existing inventory and every newly produced daughter
continue decaying. `DecayModel::prepare_production` binds the original half-life
data and topology to `fs_la::decay::source::DecaySourcePlan`; it adds no second
propagator. The prepared maps are independent of production rates and initial
populations, so one interval can be reused across both inventory and source changes.

`trial_with_production` and `trial_production_prepared` report a
`DecayProductionStep` containing original-order produced counts and the original
`DecayStep` emission ledger. The corresponding advance methods publish atomically.
A single shared physical accounting path serves closed and produced inventories:
per-node balance is final minus initial minus production plus transformations
minus incoming daughters. Source rate-times-duration roundoff and compensated
lifetime production-counter defects consume the same noncancelling nuclei
allowance. `cumulative_produced_nuclei` survives shutdown, restart and later
ordinary closed-decay advances. Initial inventories are distinct from subsequent
source production.

`DecayProductionSchedule` binds nonnegative rate vectors to the exact model Arc,
admits a caller-bounded strictly increasing change history, and performs binary
search for right-continuous rates and the next transition. Invalid rates, times,
shapes, model bindings, source products, work limits and cancellations refuse
without changing an accepted inventory. Positive but unrepresentable produced
counts are not silently credited as zero.

Production is prescribed within the tracked nuclear inventory. It does not
infer activation cross sections, fission yields, source-parent depletion, or
the bulk mass and enthalpy of a transported material stream. In particular, an
irradiated parent's exponentially falling reaction rate cannot be replaced by
its constant interval average and described as an exact source convolution.
The existing irradiation owner and any explicit product binding must preserve
that history and keep prompt capture/fission energy distinct from later
radioactive emissions.

All 60 default-feature source tests passed in the actual native `fs-nuclear`
crate root linked to the applicable unchanged production dependencies: 52 prior
tests plus eight production-owner regressions. The additions check independent
constant-source and repeated-rate ingrowth integrals, shutdown, stable and zero
intervals, partition/restart, transition admission, noncancelling per-node birth
roundoff, source-product defects, cancellation and prepared-map reuse. Full
workspace Cargo verification still requires the missing constellation checkouts.

`examples/decay_production.rs` demonstrates a source increase followed by
shutdown, using one prepared quarter-second map across changing source vectors.
It reports produced populations and separate lifetime emissions; the authored
schedule aligns with its explicit grid. Integrating consumers must clip
intervals at arbitrary schedule transitions.

The example also compiled and executed successfully through the same exact-source
native build, including its production/conservation/shutdown assertions.

## Scheduled production on the shared nuclear/thermal transaction

`thermal::DecayProductionHeating` is a specialization of the existing
`NuclearHeating<S>` owner. `DecayProductionCellSource` binds a radioactive
inventory and an immutable `DecayProductionSchedule` to the original thermal
cell/material. The common scheduler clips steps at authored source changes;
the original final accepted duration drives both the production convolution
and integrated deposited-heat credit. Source shutdown preserves radioactive
populations, daughter ingrowth and later heat release.

The source caches its prepared production map by duration and propagation
limits. Changing the rate vector alone reuses the same immutable map. A clone
shares model, schedule and prepared data while retaining an independent accepted
inventory. The original shared material bindings, caloric identity, storage basis,
geometry, temperature admission, checkpoint propagation and rollback remain in
one owner. Deposited energy remains distinct from radiation and escaped energy,
and production does not invent a bulk-material mass or enthalpy influx.

The common source seam now has default-zero produced/decayed counters and an
optional advance-wide nuclei allowance. Existing fission, irradiation and
kinetics sources retain their original behavior. Production enables the nuclei
gate, which sums per-node physical defects and scalar counter-reduction defects
across every cell and interval. Nonnegative defect totals preserve a positive
low part at an allowance boundary. The source's original energy budget and the
native thermal energy budget remain independently mandatory. Typed decay
refusals are retained by `FissionThermalError::Decay`.

`fs_couple::radioactive` re-exports the production types and shared thermal
specialization under the existing `nuclear-thermal` feature. No additional
thermal solver or independent source scheduler is introduced.

All 90 native tests passed with the actual nuclear crate root and `thermal`
enabled: 60 core, eight closed-decay thermal, six production thermal, and 16
existing fission/irradiation/kinetics thermal regressions. The production checks
cover schedule shutdown clipping and continuing heat, newborn daughter heat,
rate-change map reuse and restart, composition with prescribed fission, joint
rollback, whole-advance population allowances, and zero-source compatibility.
The dependency build used current production thermal runtime and exact boundary/
composition definitions at c9d1b20e, excluding unrelated module declarations.
This focused source run does not replace the unavailable full workspace gate.

The constant-forcing interpretation uses the standard matrix-exponential
convolution and related phi functions; see [Al-Mohy and Higham (2011)](https://eprints.maths.manchester.ac.uk/1591/).
The implementation retains this project's positive DAG uniformization and its
own explicit numerical admission rather than claiming the paper's general
matrix-exponential error analysis.

## Product-resolved fission, capture and radioactive inventories

`burnup` represents an explicitly supplied network of unique `NuclideId`s,
admitted spontaneous-decay laws, disjoint neutron reaction channels, product
yields and event energy. The same frozen linear reaction operator evolves fuel,
capture products, fission products and radioactive daughters together. Products
created during an interval can decay or react again inside that interval.
Integrated parent removals come from `fs_la::decay::reaction::ReactionPlan`;
each channel's relative hazard attributes those removals to actual events.

Fission product yields are expected multiplicities, never normalized branch
probabilities. The network may contain sinks, repeated transformations, cycles
and self-return reactions. Nodewise production/removal residuals are checked;
total tracked nucleus count is not falsely conserved across fission. All
unstable nuclides must declare their complete spontaneous-decay branching,
with probabilities summing to one. Empty product lists explicitly export the
products from this tracked inventory. Particle mass/charge accounting and
composition-dependent thermal material laws remain outside this owner.

Microscopic cross sections use square metres and group-integrated absolute flux
uses neutrons/m²/s. The full named energy grid and edges must match. The original
fixed-coefficient variants stay fixed inside the explicit Kelvin range; the
shared tabulated variants interpolate supplied temperature rows as described
below. No Doppler broadening, self-shielding or source normalization is inferred. Energy-dependent yields can
be supplied as disjoint partial channels; duplicated total cross sections would
double-count a physical reaction and are not an evaluated-data admission.

Each reaction's energy belongs to that event alone. Fission/capture energies
exclude future daughter-decay energy; the explicit product populations generate
that later release, including after zero-flux shutdown. An all-inclusive fission
Q value or an additional effective decay-heat fit would count the same energy
twice and must not be assigned to this network. Deposited, transportable radiation
and escaped channels remain disjoint. Cumulative event and energy accounts use
the existing compensated accumulator. Absolute noncancelling population and
energy defects consume separate user-supplied allowances.

`PreparedBurnup` binds immutable interval maps to the exact model, interval,
channel rates and work limits and can serve multiple inventories.
`BurnupInventory::trial` retains one shared immutable map for exact repeated
operator inputs; current flux, temperature and work are still admitted before
reuse. Preparation cost is paid only on a cache miss; application and retained
maps are quadratic in node count. Changed limits cannot bypass numerical
preparation. Candidate
state, cumulative accounts and the cache are published only after acceptance.

The formulation follows the production-minus-loss depletion equation described
in [OpenMC's depletion methodology](https://docs.openmc.org/en/latest/methods/depletion.html).
This is a bounded dense reduced network using the native positive propagator,
not OpenMC's CRAM or a spatial neutron transport solver. Accounting admission,
series truncation, coupled timestep error and nuclear-data uncertainty are
separate claims. No evaluated isotope dataset, criticality prediction or
reactor-design qualification is supplied.

## Product heating and composition feedback on the shared thermal timeline

`thermal::BurnupHeating` specializes the existing `NuclearHeating` transaction
for a complete product-resolved inventory per cell. The original temperature
and final accepted interval feed the reaction owner; only integrated deposited
energy enters the original thermal heat path. Capture, fission and product-decay
counts remain distinct. Flux changes clip the shared timeline, zero flux stops
neutron reactions while radioactive products keep evolving, and restart retains
all populations and cumulative accounts. Material identity, bulk mass, geometry
and the original caloric law remain bound to the existing thermal cell.

`FluxDrivenHeatSource` extends the sealed source seam so the original neutron
kinetics bridge can operate on either the parent-only irradiation model or a
full product network. `KineticCellSource<S>`, `KineticHeating<S>` and their
report/budget types preserve the original parent-only default types. The same
kinetics population-time integral scales the reference multigroup flux. No
duplicate neutron solver, source clock, heat integrator or thermal transaction
is introduced. The common source audit and advance-wide nuclei allowance also
retain and bound product-network population defects and radioactive counts.

`kinetics::CompositionFeedback` is an explicitly identified linear nuclide
worth fit, with absolute reference/minimum/maximum populations and a coefficient
per nucleus. Nuclide identity binds a fit to the inventory independently of node
order. Its coefficient must be calibrated for the same region and neutron
normalization; a microscopic capture cross section is not itself a reactivity
worth. Initial and prospective product inventories must stay in every fit's
domain. An unknown nuclide, domain exit, exhausted resource/account allowance,
or checkpoint failure refuses before inventories and heat are published.

Scheduled controls, temperature feedback and admitted composition worth sum
to the beginning-of-interval reactivity. This closes the reduced loop from
fuel reactions to products, product worth, neutron amplitude, subsequent
reactions and material heat. Temperature/composition feedback and the mean
within-step spectrum remain partitioned approximations; decay competes with
time-dependent neutron reactions, so matching total fluence alone is not an
exact product history. Timestep refinement is required. No automatic criticality,
spatial flux redistribution, burnup-dependent delayed-neutron coefficients,
evaluated worth calibration or empirical reactor validation is claimed.

Focused regressions exercise analytic competing reaction calorimetry with
same-interval product decay/capture, fission multiplicity, additive boundary
heat, shutdown/restart, unchanged bulk mass, and composition feedback against
an independent neutron/fluence invariant. Halving the coupled timestep gives
the expected first-order convergence in that manufactured fixture. Domain and
cancellation checks preserve the original source and thermal states.

The final merged native build passed **141 tests**: 35 numerical-owner tests
(13 reaction and 22 existing decay/source) and 106 complete nuclear tests with
`thermal` enabled, including the concurrent radioactive-production work. The
new work adds 11 burnup inventory and five coupled product/feedback regressions.
The production nuclear library and the exact `fs-couple` fission-facade imports
also compiled. The build used actual math, quantity, decay and native thermal
owners; unrelated thermal module declarations were excluded from that focused
dependency build. Full workspace Cargo, the complete `fs-couple` crate, examples,
performance benchmarks and cross-ISA validation were not run for this change.

## Radiative-capture activation on the original irradiation history

`decay::activation::CaptureActivationState` binds an admitted irradiation parent
to one explicitly identified capture-product node in a `DecayInventory`. The
binding retains its dataset label and typed `(Z, A, isomer)` identity; admission
requires the parent's same Z and A+1. Decay-model string identities are never
parsed or silently matched to an isotope table.

A frozen-flux interval uses the original irradiation removal rate `k_f+k_c` and
its actual declining parent population. Two private numerical coordinates add
that parent and a fission-removal sink to the existing positive acyclic decay
propagator. The capture edge feeds the selected product at `k_c*N_parent(t)`;
all original daughter half-lives, branches and integrated transformations remain
active during that interval, including coincident parent-removal/product-decay
rates. This retains the age distribution of produced nuclei. A uniform mean
production rate would give a different inventory and heat history.

The original irradiation owner remains authoritative for parent depletion,
fission and capture counts, prompt release and effective delayed fission-energy
reservoirs. Its capture energy coefficient must represent prompt capture energy;
subsequent activation-product decay release belongs to the decay inventory.
Numerical parent endpoints and capture/fission counts are checked against the
original irradiation candidate before its actual capture count is entered in
the product production ledger. Noncancelling discrepancies and existing owner
roundoff must fit the joint absolute nucleus and energy allowances. Export
roundoff is charged separately for each energy destination, both in the shared
decay ledger and when joining irradiation/product release; opposite channel
errors cannot cancel. Positive defect totals round conservatively at an
allowance boundary.

Both histories and the prepared-map cache publish together only after the
complete trial succeeds; any refusal or checkpoint error preserves the accepted
state. The immutable map is reused for the exact reaction-rate pair, duration
and propagation limits. The two private coordinates count against node, edge
and work limits, leaving at most 62 physical product nodes. Zero flux stops new
production while accepted product populations and fission-energy reservoirs
continue decaying. Flux, temperature-domain and group-grid admission remain with
the original irradiation data.

This adapter preserves existing irradiation semantics. It does not infer isotope
cross sections, fission yields, product recapture, spontaneous decay of the
irradiation parent, neutron balance or an emitted-particle mass/charge ledger.
Use an explicitly authored reaction-network owner for those different model
requirements; this binding does not silently replace either original history.

The separation of prompt capture and subsequent induced-radioactivity emission
follows the distinction described in the [IAEA neutron activation training
material](https://nucleus.iaea.org/sites/connect/RRIHpublic/CompendiumDB/Shared%20Documents/Malaysia/Protocols%20in%20PDF/MalaysiaProtocol_NeutronActivationAnalysis.pdf).
Coefficients and half-lives remain caller-supplied, identified data; synthetic
verification fixtures are not a calibrated isotope library. Regressions cover
analytic equal/distinct rates, competing capture/fission energy, shutdown,
interval composition, restart/cache reuse, channel export and atomic binding,
budget, range and cancellation refusals.

## Capture-activation heat and neutron kinetics on the shared timeline

`thermal::{CaptureActivationCellSource, CaptureActivationHeating}` specialize
the same `NuclearHeating` transaction used by fission, irradiation, kinetics,
burnup and prescribed radioactive production. Each source binds the original
irradiation parent and its explicit radioactive-product inventory to one active
material cell and the existing multigroup neutron-flux schedule. Flux
transitions, including shutdown, clip the final thermal interval before either
nuclear history advances.

Only the coupled step's integrated deposited release is credited to heat, once.
Prompt capture/fission release, retained fission-energy reservoirs and radioactive
capture-product emission retain their original physical owners. Radiation and
escaped channels remain in the source account. Existing authored boundary heat,
material calorics, fixed cell mass, geometry, storage basis and contact/boundary
transport continue through the original thermal runtime.

The adapter reuses the activation state's immutable propagation cache, carries
both nuclear histories through thermal cloning/restart, preserves original
checkpoint errors and reports typed activation refusals. The original
irradiation temperature domain is admitted at binding and accepted thermal
endpoints. The shared advance-wide audit includes captures, fissions, produced
nuclei, actual radioactive transformations and noncancelling population/energy
defects; all accepted intervals share the caller's combined allowances.
A thermal or nuclear refusal publishes neither candidate history.

The same source implements `FluxDrivenHeatSource`, so
`KineticHeating<CaptureActivationCellSource>` uses the existing normalized
neutron-population integral, controls and temperature feedback to drive its
incident spectrum. Scheduled and neutron-driven paths call the same activation
trial and accounting implementation. Composition reactivity defaults to zero;
no nuclide-worth fit is inferred. For changing neutron amplitude, the shared
kinetics bridge retains its mean-flux approximation: parent removal follows
integrated fluence, but within-interval product ages and decay heat require
coupled timestep refinement. Exact frozen-rate activation does not make that
changing-rate thermal partition exact.

The `fs_couple::radioactive` facade exposes the capture binding, state, budget,
report and heat adapter under the existing `nuclear-thermal` feature. Five
native thermal regressions exercise analytic capture/shutdown heat, competing
fission against the original irradiation owner, domain/grid/checkpoint refusal,
cached restart, and analytic constant-amplitude neutron-driven activation with
late-checkpoint rollback. These preserve the concurrently added full burnup and
generic kinetics owners; explicit burnup fission products and effective
fission-energy reservoirs remain distinct model choices.

## Temperature-dependent microscopic reaction tables

`cross_sections::CrossSectionTable` retains identified processed microscopic
coefficients on an immutable named energy grid and strictly ordered absolute
temperature knots. Nuclide and reaction identity are explicit. Coefficients
are nonnegative square metres; rates use absolute group-integrated scalar flux
in neutrons/m²/s. Interpolation is linear in Kelvin between supplied rows,
preserves exact knots and constant rows, and refuses temperature extrapolation.

One binary search brackets a temperature, followed by an allocation-free group
pass. Flattened rows and their temperature/grid data are shared on clone.
Near-endpoint interpolation computes the two temperature distances separately
and uses a fused update from the nearer endpoint. This preserves small
coefficients beside much larger endpoint values. Positive unrepresentable
interpolants or rate products refuse. Explicit temperature/group/value bounds
precede retained allocation; construction and rate evaluation poll cancellation.
Zero flux still requires the complete data-domain and grid checks.

`IrradiationModel::new_tabulated(TabulatedIrradiationData)` binds fission and
capture tables to the same parent and energy grid, with full coverage of the
declared model temperature domain. The existing fixed constructor remains.
All inventory, delayed-energy convolution and capture-activation paths evaluate
the current temperature through the original `rate_constants` owner. The legacy
`cross_sections()` observation returns fixed coefficients or a tabulated model's
lower-domain reference row; `cross_section_reference_temperature_k()` makes
that reference explicit. Physical evaluations must use `rate_constants`.
Original tables remain available through `temperature_tables()`.

Supplied rows must already embody the intended processing and within-group
spectrum assumptions. Linear temperature interpolation does not perform
resonance reconstruction, Doppler broadening, self-shielding or transport.
This deterministic multigroup law is distinct from the statistical temperature
selection described in [OpenMC's cross-section methodology](https://docs.openmc.org/en/latest/methods/cross_sections.html).
Dataset identity does not certify nuclear-data accuracy.

Nine table regressions and three irradiation/activation regressions passed in
the actual-owner nuclear build. They cover analytic temperature-dependent rates,
exact knots, positivity and endpoint scales, limits/identity/range/cancellation,
actual fuel competition, the original delayed-heat convolution and existing
capture-product evolution. Full workspace Cargo and empirical nuclear-data
validation remain separate.

## Temperature-dependent depletion on the original coupled timeline

`BurnupRate::TabulatedFission` and `TabulatedCapture` use the same immutable
`cross_sections::CrossSectionTable` as parent-only irradiation. Construction
binds the table's isotope, reaction, named neutron grid and complete temperature
coverage to its parent and model. The original fixed-coefficient variants remain
available. Rates are collapsed from the current material temperature and supplied
absolute multigroup flux; the reaction operator is frozen over the accepted step.

The prepared-map cache compares the entire channel-rate vector, exact interval,
model identity and work limits. A changed fission/capture split invalidates it
even if total removal is unchanged. Different temperatures or spectra that yield
identical channel rates share the same immutable map. Each call still admits
its current temperature, grid and work before reuse; a cache miss collapses the
rates only once. Preparation's recorded flux and temperature describe where
that map was first built, not every later state that can share its operator.

No second thermal or neutron solver is introduced. The existing irradiation,
capture-activation, burnup and kinetic heating adapters pass the material's
actual beginning-of-step temperature to their original reaction owner. Their
shared transaction checks the prospective endpoint against the data domain
before publishing either heat or populations. Temperature-dependent nuclear
data, material heat, product worth and neutron amplitude therefore participate
in the existing partitioned feedback loop. This remains a first-order coupled
approximation, not an implicit temperature/transport solve.

Four focused burnup regressions cover analytic temperature-dependent event and
product counts, reaction-channel cache invalidation, reuse at equivalent rates,
and table identity/grid/domain admission. Three shared-timeline regressions
compare both original reaction owners against an exact nonlinear cooling-feedback
solution, refine the neutron/temperature/product-worth loop against an independent
invariant, and verify joint rollback on a final-temperature domain exit. The
combined native nuclear lane passed 133 tests when this increment was prepared.
These are manufactured verification problems, not evaluated-data validation.

## Processed ENDF/PENDF fission and capture data

`endf::PointwiseCrossSection::parse` reads caller-supplied text and selects an
explicit MAT, isotope/isomer, total-fission (MT18) or radiative-capture (MT102)
reaction and exact positive processed temperature. MF1/MT451 supplies ENDF-6
incident-neutron identity, library/version/release/revision and temperature;
selected MF3 identity must agree. LRP=0 or 2 is required. LRP=1 MF3 resonance
backgrounds are refused rather than admitted as complete reaction curves.
Natural-element mixtures, primary zero-temperature evaluations and unsupported
TAB1 laws require an external processor. No network access or built-in isotope
data is used. The format is specified in the official
[ENDF-102 (2023) manual, §§0.5–0.6, 1.1 and 3.2–3.3](https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf).

The reader handles fixed-width records, implicit signed exponents, histogram
and linear-linear interpolation regions, and paired equal-energy points for
discontinuities. It checks selected-section control records and SEND payloads,
rejects ambiguous selections, and converts eV to joules and barns to square
metres. Source Q values remain observable metadata; neither QM nor QI is
assigned to the prompt deposited-energy channel. Event energy and later
radioactive-product release retain their existing explicit owners.

`collapse_flat_in_energy` computes each group mean by integrating every
intersected histogram or linear segment. Its declared weighting is flat flux
density per unit energy inside each group. Thresholds and discontinuities are
integrated rather than sampled at group midpoints, and source energy coverage
must contain the full group grid. `collapse_temperature_table` assembles curves
with the same isotope, reaction and evaluation revision at strictly increasing
temperatures into the original shared table. Each source curve retains its own
identity; the caller identifies the assembled table and its grouping assumption.
A single forward traversal visits the groups and source segments, while point
queries use binary search. This is an import-time transformation, outside the
reaction evaluation and prepared propagation paths.

Byte, record, point, interpolation-region, group, retained-value and segment-visit
limits bound work and allocation. All temperature rows share the collapse-work
allowance. Checkpoints are polled during scanning and integration. Conversion,
interpolation and integration refuse unrepresentable positive quantities instead
of silently turning a physical coefficient into zero. No energy or temperature
extrapolation, resonance reconstruction, Doppler broadening, self-shielding,
spectrum calculation or evaluated-data certification is supplied.

Nine importer regressions exercise identity/temperature selection, malformed
controls, numeric fields and SI conversion, histogram/linear integration with
thresholds and discontinuities, work/cancellation refusal, small endpoint tails,
and imported fission/capture tables driving actual competing burnup and explicit
prompt calorimetry. The final focused native run compiled the production nuclear
crate with thermal support and the exact public fission facade, including both
new module exports. All 147 nuclear tests passed, including the concurrent
capture-activation and thermal additions. These manufactured fixtures establish
implementation behavior; full workspace Cargo and empirical validation against
processed evaluated nuclear data remain separate gates.

## External nuclide production during full burnup

`burnup::production` adds `PreparedBurnupProduction`,
`BurnupProductionStep`, `BurnupProductionSchedule` and
`BurnupProductionChange` to the original physical burnup owner.
`BurnupInventory::trial_with_production` and its atomic advance/prepared
variants evolve an explicit expected-nuclei-per-second vector under the current
admitted flux and temperature. Supplied fuel and products can decay, capture,
fission and revisit cycles during the same interval. Half-lives, reaction
multiplicity and integrated event energies retain their original definitions.

The shared `ReactionSourcePlan` integrates actual source ages; supplied counts
are not inserted as one end- or beginning-of-step dose. One private channel
compiler and one physical event/product/energy ledger serve both homogeneous
and supplied evolution. Neutron-group collapse and temperature-table evaluation
occur once on a preparation miss. No dense homogeneous maps are duplicated to
obtain event attribution.

Prepared source maps compare exact interval, model identity, work limits and
the complete current channel-rate vector. Changing the feed vector alone reuses
the map; changing the fission/capture split invalidates it even when total
removal is unchanged. Equivalent spectra/temperatures may share the same
operator after their current data domains and rates are admitted. Accepted
mode changes retain one cache family, and failed trials preserve all accepted
physical counters and caches.

`cumulative_produced_nuclei()` retains compensated external-supply accounts
separately from reaction-created products. Node balances include supplied
counts explicitly; production products, source-counter errors and existing
reaction attribution consume the absolute nuclei allowance. Positive defect
totals retain low parts conservatively at representable tolerance boundaries.
No total-nucleus conservation is imposed on fission or other multiplying
channels. Production introduces no inferred incoming enthalpy, emitted-particle
account or bulk material-mass change; these remain separate physical inputs.

Schedules are bounded, exact-model-bound and right-continuous, with binary
search for rates and the next transition. Consumers must split at production
changes. All rates are constant within an admitted interval, including the
existing frozen flux/temperature operator; this extension does not remove the
coupled changing-rate timestep requirement.

The default-feature actual-root native gate passed all 112 nuclear tests,
including eight source regressions for saturation/shutdown, newborn fission
products and decay heat, cyclic/self-return counts, partition/restart/cache
reuse, temperature-dependent channel changes, model/schedule/work binding,
atomic budget/grid/cancellation refusal and finite integrated supplies despite
unrepresentable sums of instantaneous rates.

## Scheduled production coupled to heat and neutron kinetics

`thermal::BurnupProductionCellSource` binds the original `BurnupCellSource`
to an exact-model-bound `BurnupProductionSchedule`.
`BurnupProductionHeating` uses the existing `NuclearHeating` transaction,
material binding, inventory cache and event-energy accounts. Each accepted
interval integrates continuous supply and subsequent decay, capture and fission
together. The next source boundary is the earlier of the neutron-flux and
production transitions. Either history may stop independently; inventories
continue their admitted reactions and radioactive daughters continue releasing
their own decay energy.

The source also implements the existing `FluxDrivenHeatSource` interface.
Its backward-compatible `trial_flux_at` hook carries absolute start time from
the neutron-kinetics owner, including initial admission and each reaction trial.
Existing flux-only sources keep their original behavior through the default
implementation. Scheduled production selects rates at the actual start time;
an untimed trial with a changing production history, or a timed trial crossing
its next production change, is refused. Constant production remains usable
through the original untimed hook.

The wrapper retains the original caller-supplied composition-worth law and
checks its applicability on both the accepted and proposed inventories before
publishing heat or history. Source-count reduction error is charged in addition
to the original reaction ledger. Failed work, domain, budget and cancellation
checks preserve accepted state. Incoming nuclei do not imply incoming enthalpy
or a change of thermal cell mass; deposited reaction energy and explicit
boundary heat retain their existing separate accounts.

Neutron-driven use retains the original mean-flux approximation and frozen
temperature during an interval. Exact integration of continuous supply under
that admitted operator does not make the changing-rate coupled problem exact;
flux, composition and thermal feedback still require timestep refinement.

Five focused thermal regressions cover analytic competing-reaction heat,
independent production/flux changes and shutdown tails, scheduled production
under neutron kinetics with nonzero composition feedback, restart/cache reuse,
and atomic model/work/domain/callback refusals. The final actual-root native
nuclear run with thermal enabled passed all 160 tests, including the eight new
physical production tests. The separate numerical decay/reaction gate passed
44 tests. These manufactured checks do not substitute for the unavailable full
workspace Cargo gate or empirical evaluated-data validation.

## Spectrum-weighted processed-data collapse

`endf::spectrum::FluxSpectrum` admits an explicit identified incident flux
density in neutrons/m²/s/J on its own energy knots, with histogram or linear
interpolation and paired discontinuities. `collapse_with_spectrum` returns
both the original `GroupFlux` and microscopic coefficients satisfying
`sigma_g = integral_g sigma(E) phi(E) dE / integral_g phi(E) dE`.
This preserves the source reaction integral for the supplied within-group
spectrum, following the scalar-flux weighting in
[OpenMC's MGXS introduction](https://docs.openmc.org/en/v0.12.1/examples/mgxs-part-i.html).
It does not infer a spectrum or a scattering/fission-emission distribution.

Integration traverses both knot meshes, retains thresholds/discontinuities,
and integrates the product of each admitted pair of segment laws analytically.
The four linear-basis product weights are nonnegative. Exponents remain separate
through compensated accumulation, so a representable narrow-interval integral
need not overflow merely because its raw endpoint product does. Final fluxes,
integrals, coefficients and coefficient/flux products undergo representability
checks. Every group requires positive integrated flux and complete source-domain
coverage; no coefficient is invented for a group with zero incident flux.

The spectrum is integrated on its own mesh once. Temperature-table assembly
uses that same group flux and the same supplied spectrum at every processed
temperature, even when the cross-section meshes differ. All rows share explicit
work limits and cancellation. The original flat-in-energy operation remains
available for callers that explicitly choose that assumption.

Seven focused regressions passed in the actual-owner nuclear build: analytic
nonuniform weighting, mixed/mismatched meshes and jumps, flux-scale invariance,
representable extreme-scale products, domain/zero-flux/work/cancellation
refusals, and processed ENDF temperature rows driving original burnup counts,
products and prompt heat. These are manufactured verification fixtures; the
supplied spectrum and processed evaluation remain caller-owned physical data.

## Spatial multigroup neutron diffusion and delayed emission

`diffusion` retains extensive neutron counts per cell/energy group and explicit
stationary precursor counts per cell/family. Scalar group flux is `v_g N_cg/V_c`.
An authored macroscopic reaction removes incident neutrons at `Sigma_g v_g`;
capture, group scattering, prompt fission multiplicities and delayed-family
yields remain distinct channels. A precursor disappearance emits its declared
normalized group spectrum. Fission yields are expected multiplicities and are
never normalized as probabilities. No equilibrium precursor history is inferred.

Reciprocal interfaces carry `K_g(phi_a-phi_b)`, with K in square metres. Two
positive transfers give the same conservative current across unequal-volume
cells. A planar orthogonal face can supply `K=A/(d_a/D_a+d_b/D_b)` using the
physical diffusion coefficients D in metres and centre-to-face distances.
Exterior conductances explicitly remove neutrons; omitted faces are reflecting.
Boundary extrapolation lengths and incoming exterior sources are not invented.
The conventional multigroup/precursor equations and the limits of the diffusion
approximation are described in the original
[Moltres theory documentation](https://arfc.github.io/moltres/getting_started/theory.html).

The existing positive reaction owner propagates the frozen spatial/group operator
and its population integrals. Physical reaction, face-transfer, leakage and
precursor counts are attributed independently; fluence includes particles born
and transferred during the interval, including zero-removal receiving groups.
A single immutable preparation is shared across clones and reused for exact
model/duration/work inputs. Coefficient rebinding preserves volumes, ordered
face/boundary correspondence, grid, speeds and complete precursor identities.

Work, coefficients, allocation, cancellation and numerical particle-account
allowances are explicit. The exported noncancelling defect is rounded upward
when its retained low part exceeds the ordinary scalar export. Trials preserve
all accepted populations on refusal. The current dense numerical owner limits
the total neutron-plus-precursor population nodes to 64; this is a bounded
spatial diffusion model, not a scalable angular transport or criticality solver.
Cross sections, speeds, conductances, emission spectra and their physical
applicability remain explicitly supplied data. No precursor advection or
resonance/self-shielding processing is inferred.

Twelve focused regressions cover analytic absorption/leakage and multiplication,
down-scattering into a zero-loss group, unequal-volume conservative exchange,
delayed release and same-interval precursor ingrowth, self-scattering, interval
composition, cache/rebinding behavior, work/cancellation and representability,
and upward defect export. The final merged actual-owner nuclear lane passed
188 tests, including the separately integrated burnup/thermal paths and current
uniform-production additions. Full workspace Cargo and evaluated reactor
benchmarks remain separate validation gates.

## Spatial neutrons, finite burnup and original material heat

`diffusion::burnup::DiffusionBurnupState` binds every neutron-induced channel
of each original `BurnupInventory` to one diffusion reaction per incident group.
Grid and fission/capture identity must agree; missing or duplicate bindings
refuse. Unbound diffusion channels are explicitly fixed background interactions
whose target inventories are not tracked by this coupled owner.

At the actual cell temperature, each bound channel uses the original fixed or
tabulated microscopic data and `Sigma = mean_parent_population * sigma(T)/V`.
A bounded relaxed Picard iteration alternates the same initial neutron history
and the same initial product inventories: diffusion supplies integrated fluence,
original burnup consumes its mean group flux, and integrated parent removals
supply the next mean parent populations. Products born and transformed during
the interval remain in the original burnup propagator. A candidate is accepted
only when the sum of absolute neutron/burnup reaction-count mismatches fits the
separate caller allowance. Failed iterations never become physical substeps.

The mismatch comparison retains compensated low parts across incident groups,
and positive per-cell/global defect exports round upward where necessary. A
small group's reaction count therefore cannot disappear next to a much larger
channel count. Nonconvergence, numerical/work refusal or cancellation publishes
neither the neutron field nor any fuel/product history. This closes the discrete
channel accounts within the supplied allowance; frozen material temperature,
within-group shapes and mean within-step flux remain model approximations.
Timestep refinement and suitable nuclear data remain necessary.

`thermal::DiffusionCellSource::bind` maps the complete spatial cohort to distinct
original thermal cells and materials. The shared `NuclearHeatSource::trial_sources`
hook advances the neutron field once over the final interval selected by the
original thermal owner. Its default preserves independent source laws. Cell
handles share one immutable accepted cohort through Arc; a partial or mixed
cohort refuses. Number-density volumes must match the original cell geometries.
No second heat integrator, source clock, caloric law or bulk-mass update is used.

The existing `NuclearHeating` transaction credits each cell's integrated
locally deposited energy, keeps radiation and escaped channels separate, checks
the prospective material-domain endpoints, and publishes all cell histories
atomically. Source scheduling, current uniform-production/kinetic adapters and
native thermal boundary/contact work retain their original owners. The spatial
cohort currently tracks closed nuclear inventories; supplied-nuclide schedules
remain available through the existing production source types.

Six coupled-state regressions cover finite-fuel neutron multiplication and
calorimetry, an independent nonlinear analytical solution under refinement,
temperature-dependent capture and spatial current, full channel binding,
nonconvergence/domain/cancellation rollback, and retained small-group mismatch.
Three native thermal regressions cover one shared neutron history across cells
in different index orders, additive boundary heat, unchanged mass and calorics,
geometry/coverage refusal, and joint rollback. The final actual-owner run
compiled production fs-nuclear with thermal support and the exact fs-couple
fission facade, then passed all 188 nuclear tests, including concurrent
uniform-production additions. The merged numerical owner separately passed all
50 tests. Full workspace Cargo, large-scale transport performance and evaluated
reactor benchmarks were not established by these focused source builds.

## Evaluated radioactive-decay records

`endf::decay::DecayEvaluation::parse` imports one explicit MAT and
nuclide/isomer from caller-supplied ENDF-6 MF8/MT457 text, using the original
bounded ENDF scanner. NST distinguishes explicit stability from a finite
half-life. Radioactive half-lives enter the shared `DecayLawV1`, retaining
the supplied seconds. The format follows
[ENDF-102 (2023), §§8.4–8.4.2](https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf).

NC3/NC17 mean emitted energies, their uncertainties, mode RTYP/RFS/Q and
branching probabilities are retained in SI. The unknown-energy sentinel
becomes `None`. Detailed components in the seventeen-entry form overlap
the first three aggregate categories; they are not additional heat.
Parent-wide means do not supply branch-conditional spectra. Q values and
emitted energies imply no local deposition. Source uncertainties remain
data, not certified intervals or assumed distributions.

Raw branching probabilities and their compensated sum are never normalized;
physical assembly must still admit complete branch closure and daughter data.
The caller's source label and MF8 identity are retained declarations, with no
claim of independently verified MF1 library provenance or an isotope catalogue.

Discrete, continuum and covariance records are traversed by their LIST/TAB1
counts through SEND without allocating spectral arrays or interpreting their
radiation transport. Input bytes, all scanned records and retained modes have
explicit ceilings checked before data-dependent allocation. Energy storage is
bounded by seventeen pairs; checkpoints preserve the caller's original error.

Seven regressions cover the pinned published ENDF/B-VIII.0 tritium/helium-3
files, unknown NC17 components, spectral/covariance traversal, malformed
identity/count/padding/SEND, exact resource and cancellation limits, and
unrepresentable input. The default-feature actual-root native nuclear gate
passed all 126 tests. The fixture source and exact original Git blobs are
recorded in `tests/data/endf_b_viii_0/README.md`. These checks establish
import interoperability and numerical behavior, not empirical nuclear-data
validation.

## Evaluated decay data in inventories and burnup

`endf::decay::model::EvaluatedDecayAssembly` resolves original-order
nuclide/isomer identities into the existing `Arc<DecayModel>`. It retains the
evaluations, uncertainties, routing assumptions and identity-to-population
index map. `burnup_data` emits ordinary `BurnupNuclide` and source-labelled
decay `BurnupReaction` values for composition with separately supplied capture
and fission channels. The original inventory, numerical propagator, energy
ledger and thermal adapters remain the evolution owners.

Supported RTYP sequences derive the residual daughter for beta-minus,
electron-capture/positron emission, isomeric transition, alpha, neutron and
proton emission, including immediate combinations such as 1.5 and 1.44.
RFS selects its isomer. Every daughter must be explicitly present: absent
evaluations never become stable sinks. Duplicate identities, unsupported
origin codes and spontaneous fission without resolved product yields refuse.
Branch closure, topology and model work limits are admitted by the original
`DecayModel`, without branch normalization.

Callers explicitly route each of the ELP, EEM and EHP aggregate emitted-energy
categories among deposited, radiation and escaped destinations. Fractions are
admitted as supplied and applied once. Unknown aggregate energies refuse model
assembly. Detailed NC17 components are not added again; its known neutrino
mean may optionally enter only the escaped channel. Q never supplies heat or
an inferred energy remainder. Each branch receives the parent's routed mean,
preserving expected calorimetry without inventing daughter/emission
correlations. The heavy-particle aggregate includes neutrons and recoil, so its
shared routing fraction is an explicit transport approximation.

Source uncertainties remain observable data; they are not propagated into
solver tolerances or treated as validation evidence. No emitted-particle
inventory, material-mass adjustment or nuclear mass/charge ledger is inferred.
Assembly is bounded by the original node/edge limits and cancellable; burnup
data construction preserves deterministic ordering and source mode labels.

Six focused regressions cover published tritium/helium-3 populations and
calorimetry, equal-half-life chains and partitioned evolution, residual
daughter identities, aggregate/optional-neutrino routing, admission refusals,
and imported decay combined with an authored reciprocal capture channel in
the original cyclic burnup solver. After including concurrently published
spatial-neutron modules, the default-feature actual-root native gate passed
all 150 tests. The `evaluated_decay` example also compiled and ran: at the
imported half-life it retained half the original tritium, produced the other
half as helium-3, and deposited the explicitly routed 5690 eV mean electron
energy per transformation rather than the 18590 eV Q value.

## Continuous spatial sources and moving delayed-neutron precursors

The bounded spatial diffusion owner accepts an explicit extensive neutron and
precursor influx, in expected particles/s per original cell/group or
cell/precursor family. The original positive reaction propagator evolves each
supplied particle for its actual residence time within the frozen interval;
integrated fluence therefore includes source-created particles, including in
zero-removal receiver groups. Injection is a separate balance debit, not an
initial or endpoint inventory dose. Changing source magnitudes reuses an
unchanged operator preparation.

Optional directed precursor transfers use authored volumetric flow Q [m³/s]
and donor-cell volume V [m³], giving a first-order transfer rate Q/V for every
well-mixed delayed family. The receiving cell owns subsequent delayed-neutron
emission. Internal transfers cancel in the combined precursor balance;
explicit outflows export precursors without emitting or heating inside the
model. This describes prescribed precursor transport only. It does not solve
carrier-fluid momentum, advect bulk fuel or enthalpy, infer an incoming
composition, or model an exterior return history. External return precursors
can be supplied explicitly through the source vector.

Finite burnup accepts the same neutron/precursor sources plus per-cell,
original-nuclide production rates. Newly supplied nuclei participate in the
existing simultaneous decay/capture/fission network, with event energies and
future daughter decay retained by that owner. Coupling iterations restart the
original accepted histories; they match neutron reaction counts to burnup
counts under the existing frozen-temperature, mean-flux/mean-inventory
approximation. Changing external rates requires explicit interval splitting; spatial transport,
finite fuel and temperature feedback retain their timestep-refinement requirements.

Limits remain the numerical owner's bounded dense population graph, with at
most 64 combined neutron and delayed-precursor nodes. Prescribed transport is
well mixed per cell and has no angular neutron transport, fluid-flow solution,
resonance reconstruction or self-shielding. Manufactured analytic,
conservation and timestep-refinement checks establish the implemented discrete
equations; they do not qualify a particular reactor or evaluated data library.

## Scheduled spatial sources on the native thermal clock

`DiffusionCellSource::bind_with_sources` attaches one immutable
`DiffusionSourceSchedule` to the complete shared spatial cohort. Its absolute
right-continuous transitions clip native thermal intervals, including separate
source startup and shutdown. Each accepted step evolves the spatial cohort
once, then credits only its original per-cell deposited reaction energy.
External nuclide production has a separate thermal audit account; source
neutron/precursor injection and outflow remain in the spatial particle report.
Incoming particle energy, carrier mass and feed enthalpy are not silently
credited as heat or added to the original bulk thermal material.

## Adaptive time refinement of coupled nuclear heat

`NuclearHeating::advance_adaptive_to` compares one original native interval
with two native half intervals and publishes only the fine trajectory. The
original decay, reaction, burnup, neutron, material and thermal owners perform
every trial. This estimates first-order frozen-temperature, mean-flux and
explicit heat-transport time error; it neither extrapolates a new state nor
certifies a global error bound.

`AdaptiveNuclearHeatSource` covers the existing decay, production, burnup,
irradiation, capture-activation, point-kinetic and spatial-diffusion sources.
The original `decay::thermal::DecayCellSource` also implements
`NuclearHeatSource`, admitting plain radioactive decay through this shared
owner without requiring a production schedule. Its original `DecayHeating`
API remains available.

`NuclearTimeTolerances` scales every isotope/fuel population, delayed-energy
reservoir, released-energy channel, thermal stored-energy coordinate and
temperature independently. Point-neutron amplitudes and spatial actual
particle counts have separate scales. Interval comparisons retain individual
reaction and production counts rather than subtracting lifetime totals.
Spatial reports include group reaction, directional exchange, leakage,
delayed decay, injection, precursor transfer and outflow counts, compared once
per complete shared cohort. Each cell's supplied-nuclide vector is also
compared. Relative scales use changes from the accepted state, so an arbitrary
caloric datum or constant temperature on a latent-storage plateau cannot
conceal stored-energy error.

Native CFL, source and deposition-event cuts are discovered before nuclear
evolution; a clipped comparison is retried on the actual native span.
Identical hidden subdivisions therefore cannot masquerade as refinement.
Unequal representable halves use their corresponding first-order fine-error
factor. The deterministic controller halves rejected spans and may double
well-resolved accepted spans; the report exposes counts, accepted error and
the suggested next span.

The original continuing-observer transaction keeps one compensated nuclear
account over all accepted intervals. Thermal accounting uses original and
actual birth storage coordinates and retains both signed and noncancelling
defects. The caller's whole-call thermal numerical allowance is apportioned
over duration and the bounded accepted-step count; retries do not renew it.
Explicit trial/substep ceilings and the original thermal substep budget bound
work. A smaller final/event/CFL interval may be attempted below the requested
comparison floor, but its accuracy gate is never waived.

Only excessive estimated temporal error requests a smaller interval.
Physical, material-domain, numerical-budget and checkpoint failures preserve
their original refusal and publish neither a partial accepted prefix nor
trial histories. Nuclear and thermal numerical tolerances remain separate
from the local temporal estimate.

Seven new regressions cover analytic delayed-heat cooling and CFL cuts,
independently scheduled supply and neutron flux with restart, temperature
feedback against an independent invariant, latent storage and caloric datum
invariance, work/floor/numerical-budget/late-callback rollback, plain decay
through the shared owner, and one shared two-cell spatial cohort. The final
native integration gate passed all 228 tests, compiling the production nuclear
crate with thermal support against the actual numerical owners and the
selected original thermal runtime. Full workspace Cargo and empirical reactor validation remain
separate gates.

## Finite radioactive inventories across material compartments

`burnup::transport::BurnupTransportState` binds complete existing
`BurnupInventory` histories to explicit positive control volumes. All
compartments share the same admitted `BurnupModel` and original nuclide order.
Directed per-nuclide first-order rates transfer nuclei from finite donors;
separate exterior outflows remove tracked nuclei without subsequent local
decay or heat. An authored well-mixed volumetric flow Q uses Q/V_donor.
Selective extraction is expressed by different per-nuclide coefficients.
No carrier-flow solution, bulk material-mass change or transported sensible/
latent enthalpy is inferred.

Local decay, capture, fission, movement and nonnegative continuous external
production enter one block operator on (compartment, nuclide). The existing
`ReactionPlan`/`ReactionSourcePlan` supplies positive populations, removals
and integrated residence. Nuclei born or transferred during an interval can
react and move again during that same interval. The physical formulation
follows the coupled depletion/transfer equations described in the original
[OpenMC depletion methods](https://docs.openmc.org/en/stable/methods/depletion.html#transfer-rates);
this implementation uses FrankenSim's existing numerical owners.

Nuclear-event and movement counts are attributed directly from their rates
times integrated donor residence. This keeps a small decay channel separate
from a much larger transport count. The original shared inventory ledger
continues to own reaction products, per-event energy and compensated lifetime
accounts. `BurnupStep::parent_removals` remains nuclear-event removals;
transfers and exterior outflows have their own per-path/per-nuclide accounts.
Only true external supply enters production counters. Earlier emitted heat
and reaction histories stay in the compartment where those events occurred.

One immutable preparation is shared across cohort clones and reused for
matching complete channel/transport rates, duration and work limits. Initial
populations and supply magnitudes may change without rebuilding a compatible
map. Closed intervals use the existing homogeneous residence path; a retained
source-capable map also admits zero supply. Current temperature, flux, data
and work domains remain checked before cache reuse.

Per-cell balances include imports and exports. An independent whole-cohort
balance cancels internal transfers and retains exterior outflow, nuclear
products and supply separately. Numerical removal attribution, event products,
scalar exports, lifetime counters and collective residuals contribute to
noncancelling extensive allowances. Whole-cohort defect totals reduce the
per-cell reports, so thermal consumers charge each contribution once.
Every refusal or checkpoint failure publishes neither a partial compartment
state nor proposed movement rates.

Path counts, rate coefficients, joint graph shape and numerical work are
bounded. The original dense numerical owner still limits the flattened joint
population graph to 64 nodes. Rates are frozen within an admitted interval;
changing conditions require explicit splitting or the coupled adaptive driver.
Extreme cyclic stiffness or unrepresentable contributions retain the original
numerical refusal. These are local accounting/admission guarantees, not a
certified floating-point forward-error bound or a carrier-flow model.

Eight focused regressions cover finite donor depletion and retained old heat,
unequal-volume circulation and selective outflow, continuously supplied
particle ages, fission-product transport and competing capture, preparation
reuse and invalidation, partition/restart, slow decay amid much faster
circulation, and atomic resource/domain/callback refusal. The
`evaluated_decay_transport` example imports the original ENDF/B-VIII.0 H3/He3
records and checks all populations and local emitted heat against an
independent two-volume analytic solution over two physical half-lives.

The default-feature production-source integration gate passed all 174 nuclear
tests; the evaluated circulation example compiled and executed successfully.

## Reaction-rate-weighted fission observables

`endf::reaction_weighted` admits immutable nonnegative dimensionless incident-
energy curves and condenses them with y_g = integral(sigma*phi*y) /
integral(sigma*phi). The returned original sigma_g and absolute group flux
preserve the supplied pointwise production integral. Positive cubic basis
integration follows the union of histogram/linear knots and discontinuities;
scaled positive accumulation reuses the existing spectrum owner. A batch shares
its flux and reaction denominator passes across all yield curves.

Multiplicities above one remain unnormalized. Zero reaction groups return None;
a defined zero yield remains Some(0). There is no representative-energy lookup,
extrapolation or zero-flux fallback. Source IDs, processed cross-section metadata
and reaction identity remain explicit. Original work/range/cancellation bounds
apply, with separate output-curve/value limits. This is condensation of supplied
laws under one declared spectrum and processed state, not spectrum evolution,
resonance reconstruction, self-shielding or a data-accuracy certificate.

## Transported radionuclides on the original thermal clock

`thermal::BurnupTransportHeating` binds every compartment of the original
transport state to a distinct original thermal cell and material. Volumes must
agree with the thermal geometry. The shared source transaction evolves the
entire reaction/transport cohort once for each accepted native interval, then
credits each compartment's integrated locally deposited energy to its cell.
Radiation and escaped channels remain separate. Past emitted energy stays in
its original location; subsequent daughter events release energy wherever
their nuclei then reside.

The adapter reuses per-compartment `NeutronFluxSchedule` and
`BurnupProductionSchedule`. `NuclideTransportSchedule` supplies a bounded,
right-continuous history of the same directed transfer/outflow paths with
changed first-order rates. The union of independently authored flux, feed and
movement transitions clips the original thermal scheduler. Restart binds
current inventories at the actual absolute time and does not replay earlier
sources or movement. The thermal material continues to own its mass and
caloric law; moving tracked nuclei adds no implicit carrier enthalpy.

A small default `NuclearHeatSource::validate_sources` hook admits complete
shared collections at construction and idle/positive advances. Independent
sources retain their existing behavior without another traversal or allocation.
Transport and spatial diffusion require exactly one handle per represented
cell, one accepted cohort and one shared set of histories. Partial, mixed or
duplicate cohorts refuse before any candidate is published.

The existing adaptive controller compares individual local populations,
physical reactions, released energy and external supply. It also compares every
directed transfer and exterior outflow count in original nuclide order, once
per cohort. Thus equal total material movement cannot hide isotope or path
redistribution. Native numerical accounts, source discontinuities, thermal
stability limits, work ceilings and cancellation keep their original owners.
Local refinement does not establish a certified global error bound or solve
carrier flow.

Six coupled regressions cover analytic migration and heat location, continuous
supply ages, independent feed/flux/flow shutdowns and absolute-time restart,
complete-cohort/material/volume binding including idle calls, atomic budget/
cancellation/domain refusal, and adaptive temperature-dependent capture
against an independent nonlinear analytic solution. Existing spatial-source
tests now check earlier constructor-time refusal of incomplete or mismatched
cohorts. The transport constructor additionally checks the flattened shape and
volume length before allocating its retained volume array.

The combined native integration gate passed all 242 nuclear and thermal tests,
including the six new coupled transport regressions and the constructor guard.

## Evaluated fission-neutron multiplicities

`endf::neutron_yield::FissionNeutronYield` imports the distinct MF1/MT452 total,
MT456 prompt and MT455 delayed neutron yields from supplied ENDF-6 text. Exact
MAT, nuclide/isomer and metadata temperature select one section. The original
bounded reader admits histogram/linear TAB1 laws, preserves multiplicities
above one, converts incident energies to joules and rejects duplicate matching
sections, truncation, unsupported interpolation and excess work or storage.
MT455 with LDG=0 retains the ordered positive precursor decay constants in
inverse seconds. Family abundances and outgoing neutron spectra are separate
data; the importer does not infer either from these constants. Polynomial
LNU=1 and energy-dependent delayed constants LDG=1 explicitly refuse.

MF1 yield data can come from an unreconstructed zero-temperature evaluation.
The existing MF3 cross-section importer still requires reconstructed processed
data at its declared positive temperature. `collapse` admits their common
isotope, target mass and library/revision identity, retaining both metadata
records. It reuses the original reaction-weighted integration to preserve
`integral(sigma_f * phi * nu)`, with no representative-energy sampling or yield
normalization. A zero-reaction group has an undefined yield (`None`), not an
invented neutron multiplicity. The output retains prompt/delayed/total meaning.

Synthetic regressions exercise raw evaluation import, isomers, native spatial
prompt growth and delayed release against analytical solutions, weighted
production, incompatible data, malformed/unsupported records and late
cancellation. They establish format and equation behavior, not validation of
an evaluated isotope. The format follows ENDF-102 (2023), sections 1.2–1.4,
<https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf>.

## Evaluated independent fission products in burnup

`endf::fission_yields` imports identified ENDF-6 neutron-induced independent
fission yields from MF8/MT454, with MF1/MT451 establishing the parent isomer and
NSUB=11 source family. Product multiplicities, uncertainties and incident
energies are retained without normalization; cumulative MT459 yields cannot
enter the birth-to-decay bridge. Histogram/linear energy laws require an
unchanged explicit product set. Single-row LE=0 sources are energy independent;
multirow laws do not extrapolate.

`FissionYieldEvaluation::collapse` reuses the existing reaction-rate-weighted
owner to preserve product production integrals. Its `FissionYieldCollapse`
appends disjoint neutron-group fission channels to original `BurnupData`, with
explicit per-group prompt energy and represented or explicitly exported
positive products. It refuses duplicate fission sources and binds the fixed
processed cross sections to their exact temperature. Subsequent daughter
evolution, heat release and thermal coupling retain their original owners.
Reported uncertainties remain source metadata; covariance propagation,
temperature-dependent yield synthesis and evaluated-data accuracy are not
claimed. Format and analytical fission/decay/shutdown tests cover these paths.
The format follows ENDF-102 (2023), section 8.3.
