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
validation. Spectrum collapsing, self-shielding and temperature interpolation
remain with the supplying data/transport owner.

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
