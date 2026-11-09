# fs-nuclear contract

## Tritium production with an outgoing neutron

Two explicit processed reaction identities enter the original finite-target
and neutron-energy owners without a separate propagator:

| Identity | ENDF selection | Emitted neutrons per event |
| --- | --- | --- |
| `NeutronTritonProduction` | MT33, exactly LR0 | 1 |
| `ContinuumNeutronTritonBreakup` | MT91, exactly LR33 | 1 |

These generic identities describe neutron/triton emission, with the final
residual nuclide supplied separately. They are distinct from absorbing (n,t)
MT105. For Li-7, explicit products H-3 and He-4 give
`Li7 + n -> H3 + He4 + n`: each event removes one finite Li-7 target and
preserves neutron number while reducing neutron kinetic energy and consuming
the authored negative nuclear Q. The original diffusion binding requires
outgoing group yields summing to one, and the original activation owner retains
products, signed Q, residual deposition and later separately authored decays.
No tritium yield, outgoing spectrum or energy-routing fraction is inferred.

The MT/LR pair is checked before a curve is admitted. Generic MT91/LR0,
MF6-defined LR1 mixtures and discrete/pseudo-level breakup are not admitted as
these tritium channels. `PointwiseCrossSection::breakup_flag()` exposes the
exact validated flag retained by its immutable reaction identity. Source QM
and QI remain separate observable metadata, converted to joules; neither is
automatically used as prompt Q, a laboratory threshold or deposited heat.
This matters for breakup, where intermediate excitation and final mass
difference can have different Q values. Callers still provide the complete
representative-energy and threshold declaration.

The reaction definitions and Li-7 MT91/LR33 example follow
[ENDF-102 (2023), §§0.4.3.3–0.4.3.4, 3.2–3.3 and Appendix B](https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf).
ENDF names the largest-ZA final nucleus as the residual, making simple Li-7
neutron/triton emission MT33 rather than MT22. Evaluations can instead use
discrete levels or different breakup representations. This is an explicit
supported subset; alternate representations are never automatically added or
relabelled, and checking physical overlap before separately binding them remains
the caller's responsibility. Temperature-table assembly rejects mixed identities.

Manufactured processed records exercise both admitted representations separately
through group collapse, finite Li-7 depletion and H-3/He-4 production. An
independent nonlinear solution checks timestep convergence; neutron count,
baryon/charge and signed-energy accounts close independently. Refusal cases
cover wrong LR, silent absorption and incorrect outgoing multiplicity.
These fixtures establish the functioning pathway, not evaluated isotope data,
complete lithium reaction coverage or a total blanket breeding ratio.

## Nonfission secondary neutrons, signed prompt Q and finite targets

`DiffusionReactionKind::NeutronProduction` consumes one incident neutron and
emits an explicit expected multiplicity in each outgoing group. The original
positive reaction propagator advances those secondary particles and their
same-interval reactions; no fission identity, extra branch probability or
normalized spectrum is substituted. `neutron_producing_reactions` counts the
events and `secondary_neutrons` counts all emitted neutrons, including the
outgoing neutron in a one-neutron inelastic channel. Fissions, captures and
prompt fission neutrons retain their separate meanings.

`BurnupRate::NeutronProduction` and `TabulatedNeutronProduction` remove a
finite parent on the same original depletion operator. Supplied residual and
charged-product yields remain explicit. The processed ENDF importer retains
distinct MT16, MT17 and MT37 identities for (n,2n), (n,3n) and (n,4n); these
tables cannot enter the absorption/irradiation constructors. Spatial binding
requires a production reaction on both sides and checks that the authored
outgoing neutron yields sum to the table's declared multiplicity (one for the
explicit neutron/triton channels above, or 2/3/4 for the multiplication channels). It
never fills a spectrum, infers products or turns those events into captures.
Fixed background targets remain possible through unbound diffusion reactions.

`NeutronKineticModel::new_with_reaction_energies` admits an explicit signed
prompt Q and laboratory threshold per reaction. The event energy available
for local heat and outgoing radiation is
`E_incident + Q - sum(neutron_yield[g] * E_outgoing[g])`. The separately
supplied local heat fraction partitions that available energy. Negative Q
therefore consumes incident kinetic energy into nuclear rest mass; it is not
negative material heat or an omitted energy loss. Signed Q is retained in
per-cell, interval and lifetime accounts and enters the independent energy
residual. The original constructor retains its zero-Q capture/downscatter
admission and does not silently enable neutron production.

An active reaction must satisfy its representative-energy threshold and have
nonnegative available energy. Q and threshold are separate physical inputs:
`-Q` alone does not supply recoil-dependent laboratory threshold kinematics.
Zero-cross-section placeholders permit below-threshold groups, but active
temperature/composition rebinding rechecks the original threshold and outgoing
energy budget before publication. These checks require consistently collapsed
source data; they do not resolve an intra-group threshold or emission spectrum.

Finite activation uses the same accepted neutron/burnup history. A bound
reaction with Q in the neutron energy model, and every neutron-producing
reaction, must have zero prompt `BurnupReaction::energy_per_event`; competing
prompt-energy ownership refuses. Subsequent product decays retain their own
actual event energies. `nuclear_generated_energy_j` is the signed net of prompt
neutron Q and inventory releases. `nuclear_released` remains the original
inventory export, while neutron-associated exports already include prompt Q;
adding Q again to material heat would double count it. The conservative
neutron/material coupling-defect weight includes absolute signed Q.

Accepted source histories, shutdown/restart and existing material-heat adapters
retain these accounts. Adaptive source comparison includes each cell's signed
prompt Q separately from combined heat, so compensating changes cannot hide
temporal error. Cancellation, work/domain refusal and failed count/energy audits
publish no partial physical state. The numerical kernels, dense node/work
limits and original mean-flux/mean-target coupling approximation are unchanged.

The event distinction follows [NNDC ENDF reaction definitions](https://www.nndc.bnl.gov/endf/help.html);
the energy balance follows [OpenMC's energy-deposition methods](https://docs.openmc.org/en/stable/methods/energy_deposition.html).
Focused manufactured tests cover triangular secondary-neutron kinetics with
positive/negative Q, accepted fusion history and source-off release, finite
Be-9 to two He-4 products against an independent nonlinear solution, material
heating, wrong identities/multiplicities/thresholds and late cancellation.
Cross sections, effective speeds, residual yields, Q, thresholds and deposition
fractions remain authored. Angular transport, charged-particle stopping, fission
energy transport, evaluated blanket validation and a total tritium breeding
ratio are not supplied by this extension.

## Time-resolved fusion forcing across material cells

`fusion_neutrons::history::FusionNeutronHistory` retains the original L1
accepted intervals and their target cells. Every interval uses the existing
number/energy-conserving source projection; no new neutron spectrum is inferred.
Intervals in one cell cannot overlap, while different cells may emit together.
A bounded construction sweep prepares right-continuous original-order rates
without subtracting old intensities. Binary-search lookup borrows those rates;
the declared interval and rate-storage limits bound preparation and retained data.

`NeutronKineticState::trial_with_history` clips evolution at those source
transitions and an explicit maximum transport step, reuses the original
continuous-source neutron propagator and audits the whole requested interval.
Fast neutron residence/interaction times can be resolved independently of the
plasma recording intervals without changing speeds or material coefficients.
It resumes from the supplied checkpoint clock,
does not replay earlier births, and retains source-free slowing/capture after
the last emission. Per-cell heat, radiation, leakage, source moments and global
storage are accumulated independently. A late cancellation or exhausted work
or account budget publishes no intermediate physical state.

The material adapter's `bind_with_history` uses the same forcing through the
existing `NuclearHeating` cohort transaction. Pulse and accepted-history forcing
remain explicit variants; `pulse()` and `history()` return optional references
to the retained original data. Original thermal material and geometry admission,
source transitions and shared-clock checks remain mandatory.

Accepted intervals resolve the evolving burn only to the supplied temporal
resolution. Uniform emission within each retained interval and projection of
each channel to its mean energy remain declared approximations. Refining burn
steps is required when intra-step source variation affects neutron residence
times or material heating. G1/G3 tests compare piecewise capture with its
analytic solution, distinguish equal-yield sources with different timing,
and check concurrent cells, restart, source-off evolution and atomic refusal.
`fusion_neutron_history` runs two actual D-T histories through spatial neutron
transport and a source-off tail using explicitly illustrative coefficients.
`fusion_neutron_water` now drives the original IF97 receiver from accepted burn
intervals rather than a whole-burn constant source.

## Fusion-neutron transport and material heat

`fusion_neutrons` binds the original spatial multigroup `DiffusionModel` to
explicit representative kinetic energies and per-reaction deposition fractions.
`NeutronKineticState` reuses its positive continuous-source propagator, reaction
histories and reciprocal face transfers. Capture removes incident kinetic
energy; downscatter removes the incident/outgoing representative-energy
difference. Each removal is split explicitly between local heat and outgoing
radiation, while boundary leakage carries the remaining neutron energy.
Every cell and the full interval audit neutron storage, incoming source energy,
heat, radiation and leakage independently. Lifetime exports use the existing
compensated accumulators. Source/model/domain/work/account failure publishes
no candidate state.

`FusionNeutronSource::project_interval_mean` retains the original immutable
D-T and D-D-neutron emission. Each channel's count and energy are conservatively
projected onto adjacent authored group energies, with unsupported coverage
refused. The named approximation collapses a channel to its mean energy and
holds its rate constant over the original interval. Native continuous injection
preserves residence-time effects for that declared source; it does not resolve
the original subinterval burn history or a thermal/Doppler spectrum. Users
must split/refine emission intervals and the energy grid when those effects
matter. No channel branching fraction or automatic source normalization is
introduced.

With `thermal`, `NeutronKineticHeating` uses the existing `NuclearHeating`
transaction and original material caloric laws. One complete shared spatial
cohort advances once per accepted thermal step. `FusionNeutronPulse` retains
the original source and absolute start/end, clips native thermal steps at
source transitions, and refuses unrepresentable pulse duration. Shutoff does
not erase stored neutrons: their later scattering, absorption and leakage
continue. Complete distinct cell binding, shared accepted clock, original
reference volumes, cancellation and atomic material/source publication remain
mandatory. Only the local deposited channel becomes material heat.

The original zero-Q `NeutronKineticModel::new` is a bounded nonmultiplying
kinetic-energy interpretation with fixed authored group data. Explicit
prompt fission and other production, signed Q and thresholds use the constructor
described above; none is inferred. `NeutronKineticState` alone supplies no finite
activation/depletion, changes of material mass/composition, temperature feedback
to cross sections, abstract delayed-precursor energy or photon transport. Effective
speeds/cross sections and representative energies must have compatible physical
averaging. The original diffusion approximation does not resolve angular or
near-void streaming; see the primary
[Moltres theory](https://arfc.github.io/moltres/getting_started/theory.html)
and [OpenMC neutron-physics](https://docs.openmc.org/en/stable/methods/neutron_physics.html)
descriptions of diffusion/scattering/capture assumptions.

G1/G3 checks use independently solved constant-source absorption/leakage and
two-group slowing/capture equations, reciprocal spatial transfer, restarts,
projection count/energy, zero/source-off evolution and atomic refusal. Thermal
checks compare original material heating with direct neutron evolution and
analytic capture, including pulse boundaries and late deposited heat.
`cargo run -p fs-nuclear --features thermal --example fusion_neutron_water`
connects an actual D-T burn to the existing neutron owner and finite IF97
water receiver, checking the combined plasma/neutron/material energy account.
Its neutron coefficients are explicitly illustrative, not evaluated water or
blanket data; IF97 supplies the receiver's original caloric response.

## Finite fusion-neutron activation and delayed energy

`fusion_neutrons::activation::FusionNeutronActivationState` combines finite
target depletion and radioactive products with incident neutron kinetic
accounting. The original `DiffusionBurnupState` advances the neutron field and
inventories once; kinetic accounting consumes that same accepted neutron report.
Reaction coefficients bound to the original inventories may evolve with target
populations and authored temperature-dependent microscopic tables. Original
group identity, effective speeds, geometry, reaction order, fixed backgrounds,
face transfer and leakage remain unchanged.

Every authored capture-event energy must exclude incident neutron kinetic
energy and later product decay energy. Kinetic losses, actual nuclear-event
release, local heat, outward radiation and escaping energy retain disjoint
accounts. Products release decay energy only when they actually transform;
source shutoff does not erase stored neutrons or activated inventories. No
implicit target replenishment, isotope yields, capture Q, breeding reactions,
photon transport, bulk mass update or evaluated-data inference is supplied.

One combined allowance charges kinetic and native nuclear defects, channel
export arithmetic, the global energy residual and energy-weighted capture-count
inconsistencies without cancellation. Source, topology, temperature, work,
nonconvergence and cancellation refusal publish no candidate. G1/G3 regressions
compare the same native history, finite-depletion refinement against an analytic
solution, source-free radioactive cooling, evolving capture data and spatial
energy accounts. Interval-mean flux and target populations remain the original
owner's approximations; convergence of reaction counts is not a temporal-error
bound, and physical timestep refinement remains necessary.

## Fusion activation on the shared material clock

With `thermal`, `FusionNeutronActivationHeating` uses the existing
`NuclearHeating` transaction and original material calorics. A complete cohort
binds one accepted activation checkpoint and one immutable `FusionNeutronForcing`
to distinct original thermal cells. Pulse and accepted-history transitions clip
the native interval; restarting does not replay earlier emissions. Cell
temperatures are reordered into the original spatial order before the single
coupled trial, so authored capture tables see the correct material environment.

Only the combined deposited channel becomes heat. The accepted cell report
retains separate kinetic and nuclear reports; the original burnup counters
report finite captures, decays and supplied nuclei. Fixed-background captures
remain in the neutron report. Shared numerical energy and neutron-source
defects are charged once per cohort, independent of thermal cell ordering.
Radiation remains separate unless the existing explicit radiation-response
owner is attached. Source-off product decay continues on the same clock.

Complete coverage, original volume/material/mass binding, temperature domain,
accepted clock and atomic cancellation rules remain those of the existing
owners. The native thermal step limit also bounds frozen-temperature and
mean-flux intervals; source transitions alone do not establish temporal accuracy
or guarantee a representable stiff propagator. G1/G3/G4 thermal regressions
cover material heating against direct native activation, delayed heat with no
remaining neutrons, restart and whole-cohort refusal/cancellation.

## Adaptive fusion and material evolution

`NeutronKineticCellSource` and `FusionNeutronActivationCellSource` implement
the existing `AdaptiveNuclearHeatSource` interface. `advance_adaptive_to`
therefore refines their original `NuclearHeating` transactions against one
coarse and two fine physical intervals, publishing only the fine trajectory.
The shared driver retains source transitions, material events, native CFL cuts,
bounded work and atomic cancellation. No second adaptive controller or physical
source solver is introduced.

Comparisons use each original neutron group, finite nuclide population,
stored neutron energy, reaction count and current-interval energy channel.
Kinetic and nuclear releases are also compared separately, so opposing changes
cannot conceal a temporal defect in their combined heat. Full spatial transfer,
leakage and injection accounts are compared once per cohort, in original order;
small late releases are never reconstructed from lifetime-account subtraction.

These are local time-error estimates, separate from numerical energy/count
allowances. They do not certify a global error, increase the underlying order,
refine the authored source spectrum, or waive a native physical-domain refusal.
The maximum step must still admit the original stiff neutron propagator.
G1/G3/G4 checks cover analytic post-pulse capture with material cooling, finite
activation and temperature-feedback refinement, and unchanged accepted states
after exhausted work or late cancellation.

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
integrated-count maps. Direct branches into stable daughters use that owner's
closed-form maps, including uniform-production age factors, with no series or
matrix squaring. Genuine radioactive daughter chains keep the general positive
propagator. Both paths retain the same population and energy accounts.
Reusing a `PreparedDecay` avoids repeated matrix
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
with their existing owners. Same-node exponential survivors retain representable scaled tails through the
numerical owner's scalar attenuation kernel. Large inventories that could amplify
unrepresented off-diagonal or integrated-count coefficients retain its explicit
refusal; the physical wrapper never bypasses those range checks.

Ten focused physical regressions cover analytic repeated-rate ingrowth and
integrated emissions, distinct energy channels for duplicate daughters, activity
versus interval energy, preparation reuse and model binding, partition/restart,
tiny heat with charged unresolved depletion, stable and zero-duration identities,
input admission, early/mid/final cancellation, and energy/range refusals. All ten
passed, alongside the 19 fission-source tests, in a native Rust source build
using the unchanged production decay modules from `fs-la` and `fs-thermochem`
and the real `fs-math` and `fs-qty` crates. This isolated the applicable owners
from unavailable constellation dependencies; it is not a workspace Cargo gate.

## Finite-population stochastic radioactive decay

`decay::stochastic::StochasticDecayInventory` provides an integer realization of
the original immutable `DecayModel`. This complements the continuous expectation
in `DecayInventory`: a single nucleus survives or transforms, and daughters can
transform only after their sampled birth times. The supplied total population
and per-advance event limit are bounded by 2^53, retaining exact conversions to
the existing floating-point observation accounts.

Each occupied unstable nuclide retains an aggregate exponential clock with
hazard `lambda * N`. A transformation selects an original ordered branch,
decrements its parent, increments its daughter, and resamples only the affected
clocks. Unaffected clocks remain scheduled. Branch closure must round to exactly
one under compensated summation; probabilities are not normalized. Random draws
use the original `fs_rand::Stream` and its 53-bit sampling resolution. There are
no fixed time bins or per-event allocations. Event work scans the original
bounded node list and one parent's branches.

`new_at` accepts an explicit absolute time and caller-owned `StreamKey`. Cloning
the entire inventory preserves populations, pending compensated event times,
random-stream state, and per-event lifetime energy. Accepted observation
partitions therefore preserve the sampled trajectory and accumulated-energy
bits on the same ISA/toolchain. A random-stream checkpoint alone is not a
complete inventory checkpoint. Stable/empty inventories draw nothing, and
zero-duration observations are identities. Events at an accepted endpoint are
processed exactly once. `next_event_time_s` returns the first representable
endpoint containing the next event; positive sub-ULP waits retain their low
part instead of becoming zero-time loops. A wait beyond every finite endpoint
remains an unstable future event rather than changing the nuclide's law.

`trial_to` returns an independent candidate and report; `advance_to` publishes
only after every event, energy account, and final checkpoint succeeds.
Underflow/nonprogress, counter or RNG exhaustion, work or energy limits, and
callback refusal leave the original inventory and all future clocks intact.
Per-branch energy retains the original model's authored mean-energy
interpretation and its separate deposited, transportable, and escaped accounts.

This lane models demographic fluctuations under constant rates in the original
closed, acyclic, at-most-64-node model. It does not sample uncertain half-lives
or emission spectra, add production/transport, or certify distributional
forward error. Large populations should use the existing deterministic
propagator when individual transformations are not required. The runnable
`cargo run -p fs-nuclear --example stochastic_decay` uses explicitly
manufactured coefficients and demonstrates a whole-inventory restart.

## Scheduled Poisson production of discrete radioactive inventories

`decay::stochastic::production::StochasticDecayProduction::poisson` binds an
existing integer inventory to the original model-specific
`DecayProductionSchedule`. This constructor explicitly selects independent
Poisson arrival intensities from the authored nuclei/s; deterministic dosing
and correlated reaction products are not inferred from a mean rate vector.
Binding starts at the supplied inventory's current time without replaying past
source history. The exact model Arc must match.

One compensated pending birth clock per nuclide races the original decay clocks
and authored source transitions. Newly supplied nuclei and their daughters can
decay during that same interval, with their actual sampled ages. Source changes
at an exact tie precede old-clock events; event times differing only in their
compensated low parts retain their physical ordering. Endpoint changes execute
once. A no-op source change consumes transition work but neither changes pending
clocks nor draws randomness. Shutting off supply leaves existing decay active.

The original RNG, branch selection and per-event energy accounts remain shared.
Births add no inferred creation energy, bulk carrier mass or incoming enthalpy.
Integer per-nuclide birth counters and the original decay report remain separate.
Population, birth, decay and schedule-work ceilings are explicit. The complete
candidate, including source clocks and RNG, rolls back on any refusal. Whole
and partitioned observations and full-owner clone/restart preserve the sampled
history on the original ISA/toolchain; an all-zero schedule preserves the closed
inventory's trajectory exactly. There are no per-event allocations.

`thermal::StochasticDecayProductionHeating` binds these histories to the existing
material transaction. Births, source changes and transformations clip its native
clock; births alone deposit no heat. Only original branch emissions are credited,
using the original absolute thermal endpoint and separate energy destinations.
The existing optional radiation response can distribute the outward channel.
The example `cargo run -p fs-nuclear --example stochastic_production` demonstrates
explicit manufactured source pulses, shutdown decay and complete restart.

## Integer capture activation and radioactive decay

`burnup::stochastic::StochasticActivationState` evolves finite targets, capture
products and their subsequent radioactive daughters on one integer history.
`StochasticActivationSchedule` prepares the original `BurnupModel` reaction
rates for each explicitly authored flux and nuclear-temperature frame. Fixed
and tabulated multigroup cross sections, original nuclide identities, half-lives,
temperature domains and immediate event energy keep their existing owners.
New products can capture again or decay within the same observation interval;
their decay clocks start at their actual creation time. Capture deposits prompt
energy, while later decay events own their separate delayed energy. Shutdown
removes capture hazards and preserves the radioactive tail.

Each channel must remove one represented parent and retain exactly one daughter
with unit yield. Fission, missing products and fractional/multiple expected
yields refuse: mean yields alone do not define a joint discrete outcome. In
particular, a burnup model retaining multiple charged decay products uses the
existing expected-population owner, not this one-product stochastic lane.
Spontaneous branches close to exactly one without normalization. Total integer
population is conserved and bounded by 2^53; lifetime reaction counters are
checked against integer overflow.

The original compensated radioactive clock sampler and counter-based random
stream are reused. An event resamples only reactions whose parent population
changed; self-product channels resample only themselves. Rate changes redraw
only changed hazards; unchanged rates and observation boundaries preserve their
pending clocks. Authored forcing wins exact ties, while genuinely earlier
sub-ULP reactions retain their order. Cross sections are evaluated once per
prepared frame, and events allocate nothing. Reaction scans, frame work,
resampling, events and noncancelling energy roundoff have explicit trial bounds.
Failures and cancellation publish no partial populations, energy, clocks or RNG.
Cloned checkpoints preserve trajectories on the same ISA/toolchain profile.

Six focused activation regressions check independent competing-hazard and
equal-rate birth-age formulas, tabulated multigroup shutdown, partitions and
forcing replacement, exact/sub-ULP ordering, and atomic budget/cancellation
refusals. All 68 focused activation and existing burnup/decay/stochastic tests
passed the original-source native gate. The `stochastic_activation` example
exercises two capture generations, two half-lives, a changing incident spectrum,
shutdown and restart with exact 64-nucleus conservation. Its coefficients are
explicitly manufactured; they are not evaluated isotope data. This scope does
not add neutron feedback, automatic thermal-temperature feedback, or a joint
fission-product distribution.

## Integer activation heat and evaluated radiation

With `thermal`, `StochasticActivationCellSource` binds the accepted integer
capture/decay history to its original material cell. `StochasticActivationHeating`
uses the existing `NuclearHeating` transaction, clipping native thermal steps at
sampled events and authored forcing changes. It passes the original absolute
endpoint to the source, credits prompt capture and subsequent decay heat once,
and keeps outward radiation and nonrecoverable energy separate. Independent
cells require distinct stream keys; clock mismatch, material-domain refusal,
work budgets and cancellation retain the original material and nuclear states.

`BurnupDecayRadiation::from_stochastic_activation_step` projects the original
integer decay counts through the existing evaluated spectral partition. It
checks original model identity, interval and reaction totals; it neither
resamples photons nor re-evolves populations. The same source supports
`with_grouped_decay_radiation`, retaining group-dependent receiver deposition.
Prompt capture radiation has no invented daughter-decay spectrum: it remains
explicit ungrouped escape in this grouped-decay configuration. The original
mean-energy routing and spectral reconciliation still own the energy account.

The prepared activation schedule explicitly owns nuclear temperature and flux.
Material temperature is checked against the admitted model domain, but does not
silently redraw clocks or replace the authored irradiation history. Temperature
feedback requires an appropriately authored/refined history; this adapter does
not claim an implicit neutron/thermal solve.

Four focused coupled regressions pass: prompt-versus-daughter event calorimetry,
shutdown/restart with thermal subdivision, material/RNG rollback, and grouped
radiation with separate capture energy and receiver heat. The full thermal-enabled
nuclear production source compiled through the focused native dependency seam.
The six activation and four coupled regressions also passed after including the
concurrently published cross-section and radiation-transport source additions.
This is not a full workspace Cargo result.

## Stochastic transport and event-local radioactive heat

`decay::stochastic::transport::StochasticDecayTransportState` advances a complete
cohort of integer radioactive inventories using the existing `NuclideTransport`
description. Each directed path's coefficient is an authored per-nucleus
hazard in 1/s; the donor's current integer population determines the actual event
intensity. All cells share the exact original model Arc and absolute time, with
distinct logical decay streams and a separate explicit movement stream. This
owner models finite-inventory, well-mixed movement, without inferring geometry,
bulk carrier flow, sensible/latent enthalpy or external supply.

Original decay clocks and per-path, per-nuclide movement clocks compete at their
compensated event times. A nucleus can transfer, decay into a daughter, transfer
again or leave through an outflow within one observation interval. A change in
donor population resamples only the affected hazards at that actual event time.
Observation boundaries consume no randomness. Exact ties retain original
compartment/node order before authored path order. Whole/split observations and
cloned restarts retain complete population, clock, RNG and lifetime-account bits.

Transfers preserve the nuclide and integer cohort total. Outflows remove that
nuclide and explicitly count its departure; subsequent exterior decay belongs
outside the retained model. Movement releases no nuclear energy. Each original
branch's deposited, radiation and escaped channels remain at the compartment
where the decay actually occurred, including emissions before a later transfer.
The initial cohort is bounded by 2^53 nuclei; constructor bounds also cover
compartments, rate coefficients and clock storage. Separate trial budgets bound
whole-cohort decays, movements and clock scanning/rescheduling. Original per-cell
emission defects are rounded outward into one cohort energy allowance. Rejected
work, arithmetic, accounting or cancellation leaves the complete cohort intact.
Trial storage is allocated before event processing; events allocate nothing.

`thermal::StochasticDecayTransportCellSource::bind` creates one handle for every
original thermal cell. `StochasticDecayTransportHeating` requires that complete
cohort and advances it exactly once per grouped material trial, independent of
cell ordering. Handles share the accepted state and report; a single designated
handle advertises the cohort's next event without redundant all-clock scans.
Event-local emissions enter the original atomic material/radiation transaction,
so transfer alone produces no temperature change and rejected material updates
publish neither movement nor decay history. Scheduled Poisson supply composes
with this same transport history through the owner described below.

## Scheduled production and transport on one radioactive timeline

`decay::stochastic::transport::production::StochasticDecayProductionTransportState`
binds scheduled external supply to an existing integer transport cohort through
the explicit `poisson` constructor. Its schedule vector has one optional original
`DecayProductionSchedule` per compartment. Every supplied schedule shares the
cohort's exact original model Arc. Binding uses rates active at the cohort's
current time and preserves its existing decay/movement clocks; earlier source
history is not replayed. Poisson source clocks share their compartment's original
random stream, while movement retains its separate original stream.

The original transport event kernel races authored source changes, sampled
births, movement and radioactive transformations at their full compensated
times. This is one combined trajectory, without production/transport splitting
at rounded endpoints. A supplied parent can move, produce daughters, lose those
daughters through outflow or release their decay heat in another compartment
within the same observation interval. Source shutdown stops births while
retained nuclei continue moving and decaying. At exact ties, movement-rate
changes precede source changes and the original decay/movement ordering,
followed by births, with fixed
compartment/nuclide tie order. Unchanged source rates retain pending clocks.

The standalone production owner and the combined owner share source scheduling
and birth mutation; the original transport owner and combined owner share the
transport trial kernel. The original source-free transport path retains its
state and report, and a single supplied compartment without movement retains
the standalone production history. Observations consume no random draws, and a
complete combined-owner clone retains pending source/decay/movement events,
schedule positions, all streams and cumulative integer/energy accounts.

Every accepted interval checks the exact whole-cohort balance
`retained_end = retained_start + births - outflows`. Construction and trial
budgets include source clocks and schedule work alongside transport storage and
clock work. Birth, schedule, decay and movement bounds apply across the complete
cohort, and the retained population ceiling is enforced after each actual birth.
All per-compartment source counters, original directed movement accounts and
original local emission channels remain distinct. Supply and movement infer no
creation heat or carrier enthalpy. Event processing allocates nothing, and any
refusal retains the entire original history.

`thermal::StochasticDecayProductionTransportHeating` binds that same combined
history to every original material cell. The existing transport adapter now
uses one sealed `StochasticDecayCohort` interface for supplied and unsupplied
cohorts, sharing complete-cohort admission, the single event-boundary scan,
original endpoint handling, local heat mapping and atomic material publication.
Each cell's interval report exposes the shared complete production/transport
report through `report()` and its original local emission through `decay()`.
The manufactured `stochastic_production_transport` example demonstrates a source
pulse feeding a downstream compartment with daughter decay, outflow, integer
balance and exact cloned restart.

## Scheduled circulation, valves and extraction of discrete radionuclides

`decay::stochastic::transport::schedule::StochasticDecayTransportSchedule`
binds a bounded, right-continuous movement history to the original model,
compartments and ordered transfer/outflow paths. Each change replaces the
per-nuclide hazards at an absolute time. Zero rates close a path; reopening or
changing selective extraction needs no new inventory or scheduler. Opposite
directed paths describe circulation or reversal while retaining separate
cumulative movement accounts.

`burnup::transport::schedule::NuclideTransportTimeline` owns the common frame
representation and validation. The deterministic thermal transport schedule
reuses it and retains its existing public API and model/volume admission.
Histories reject nonfinite or negative coefficients, reordered paths, changed
coefficient shapes, nonpositive or unordered change times and excess frames.
Admission polls cancellation within bounded coefficient chunks; retained
frames are borrowed without copying rate vectors during event processing.

Both the standalone and Poisson-supplied integer transport owners accept
`with_transport_schedule`. Binding at a later time installs the currently
active rates without replaying earlier changes. Only movement clocks whose
rates actually change are resampled on the original movement stream at the
true change time. Nuclear clocks, source clocks, lifetime emissions and
unchanged path clocks are preserved. Closing an exactly tied movement clock
therefore prevents that old-rate movement; a movement strictly before the
boundary still occurs even when their outward floating-point endpoints match.

The complete event race includes movement changes, source changes, decay,
transfer/outflow and births. `transport_changes` records installed frames;
the existing clock-work budget charges the retained schedule slot, inspected
coefficients and resampled clocks. No-op frames consume bounded schedule work
without consuming randomness. Admission, arithmetic, work or cancellation
failure publishes none of the candidate clocks, counters or populations.
Whole/split observation and cloned restart retain the same accepted history.
The existing thermal cohort adapter clips at these same rate boundaries and
keeps all nuclear heat at the compartment of the actual decay event.

`NuclideTransfer::from_volumetric_flow` and
`NuclideOutflow::from_volumetric_flow` admit typed `VolumetricFlowRate` and
original `Volume` values into that same transport description. For a
nonselective, well-mixed tracer the per-nucleus hazard is `k = Q / V_donor`;
the recipient volume does not replace the donor volume. Opposite paths through
unequal compartments therefore retain distinct residence times even at equal
volumetric flow. These constructors are shared by deterministic burnup and
discrete transport, with no separate flow model or duplicate rate owner.

Referenced volumes must be positive and finite, paths must have valid distinct
endpoints, and the original nuclide count must fit the existing nuclear bound.
Zero flow closes a retained path exactly. A positive flow whose rate overflows
or rounds to zero refuses instead of becoming infinite flow or a silent closure;
representable positive subnormal hazards remain positive. These are fixed
control-volume tracer coefficients, without implicit carrier continuity, bulk
mass evolution or transported sensible/latent heat. Changed geometry requires
an explicitly rebuilt rate frame.

The manufactured `stochastic_scheduled_transport` example combines a finite
Poisson feed pulse, unequal compartment volumes, valve opening/closure/reopening,
downstream extraction and two radioactive half-lives. It checks integer
`supplied = retained + outflow` at each observation, no births after source
shutdown, every movement change once, and identical whole/partitioned/restarted
history. Thermal integration regressions separately exercise pure valve and
movement boundaries without heating, recipient-local decay and radiation,
source shutdown, arbitrary observation partitions, reversed cell mappings,
and atomic work/cancellation/material refusal.

## Stochastic decay events on the shared thermal timeline

With the existing `thermal` feature, `thermal::StochasticDecayHeating` binds
integer inventories to original material cells through `NuclearHeating<S>`.
`StochasticDecayCellSource` retains each cell/material identity and the complete
stochastic inventory. Its absolute clock must match the thermal clock at
binding; independent cells must have distinct logical random-stream keys.

The next sampled transformation is a thermal step boundary. The pre-event
material state supplies heat fluxes over that interval, and only the original
branch's deposited joules are credited at its containing endpoint. Radiation
and escaped energy remain in their separate nuclear accounts. Daughters retain
their sampled birth times; several sub-ULP events may share one representable
thermal endpoint, subject to the source's event budget. The original explicit
thermal timestep and global substep limits continue to apply.

The source receives the thermal owner's original absolute endpoint, rather than
reconstructing it as `start + duration`, which need not reproduce the same
binary64 time. Existing continuous sources retain their original integration
duration. All source, material, energy, and cancellation checks remain inside
the existing joint transaction. Failed heating publishes neither the material
candidate nor changes to inventory, pending clocks, or random-stream state.

This coupling uses the existing uniform-cell thermal model and the authored
constant half-life and mean-energy data. It does not infer environmental
half-life changes, sample radiation spectra, change bulk caloric mass, or add
particle transport. Stochastic sources do not implement the deterministic
adaptive-error estimator. Five focused regressions cover event-time heat,
energy destinations, restart/partition identity, exact endpoint clocks, stream
binding, and complete refusal rollback.

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
spectrum assumptions. The explicit `endf::doppler` operation below can supply
free-gas heated pointwise curves before collapse. Linear temperature
interpolation itself does not perform resonance reconstruction, Doppler
broadening, self-shielding or transport.
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

## Processed ENDF/PENDF fission, absorption and secondary-neutron data

`endf::PointwiseCrossSection::parse` reads caller-supplied text and selects an
explicit MAT, isotope/isomer, total-fission (MT18), radiative-capture (MT102),
charged-particle absorption (MT103–107), or nonfission neutron production
(MT16, MT17, MT37), neutron/triton (MT33/LR0), or continuum triton-breakup
(MT91/LR33) reaction and exact positive processed
temperature. MF1/MT451 supplies ENDF-6
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

The common `NeutronReaction` identity distinguishes (n,p), (n,d), (n,t),
(n,He-3) and (n,alpha). The same original tables enter tabulated irradiation and
`BurnupRate::TabulatedCapture`; `captures` counts nonfission parent-removal
events. Products, yields and immediate energy remain explicitly authored.
Original finite diffusion/burnup and fusion kinetic owners advance those
inventories together; no alternative propagator or isotope dataset is introduced.
The channel assignments follow the
[NNDC ENDF reaction definitions](https://www.nndc.bnl.gov/endf/help.html)
and ENDF-102 §0.4.3.2.

All imported MF3 channels require LR=0 except the explicit MT91/LR33 identity
described above; additional unrepresented breakup is refused. Nonfission
neutron-production tables use their distinct rate and
binding pathway described above. MF3 supplies no secondary spectrum, product
yield, threshold closure or deposited-energy model. MT105 identifies one
triton-producing reaction; it does
not imply a total breeding ratio. G1/G3 tests use manufactured processed tables
for all five channels, independent analytic depletion, and Li-6(n,t)He-4 with
explicit tritium/helium products, neutron/charge/baryon accounting, source-off
capture, and separate incident-kinetic and nuclear-event energy accounts.

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
of silently turning a physical coefficient into zero. Import and collapse do
not infer energy/temperature extrapolation, resonance reconstruction, Doppler
broadening, self-shielding, a spectrum or evaluated-data certification. Explicit
free-gas heating is a separate operation described next.

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

## Free-gas Doppler heating before native reaction collapse

`endf::doppler::PreparedDopplerBroadening` heats a complete processed
`PointwiseCrossSection` from its actual source temperature T0 to T >= T0. The
caller explicitly admits the free-gas model for this input; an arbitrary
processed, shielded or bound-atom curve is not assumed to satisfy it. The
operator uses the original target/neutron atomic-weight ratio and the increment
`theta = k_B (T-T0)/AWR`, in joules. Cooling is refused. Unchanged temperature
uses the original bitwise evaluations, domain and discontinuity conventions;
identity tabulation retains the original grid and interpolation laws.

The physical kernel is the Maxwellian average of **relative speed times cross
section**. It acts on `g(E)=sqrt(E)*sigma(E)`, with dimensionless relative speed
`x=sqrt(E'/theta)` and incident speed `y=sqrt(E/theta)`. Its normalized positive
weight is `x/(sqrt(pi)*y) * [exp(-(x-y)^2)-exp(-(x+y)^2)]`. Both image terms are
retained, with `expm1` and a regular small-y expression avoiding cancellation.
This preserves a 1/v law and produces the physical low-energy 1/v tail of a
constant cross section. It follows the primary
[NJOY BROADR free-gas formulation](https://github.com/njoy/NJOY2016-manual/blob/master/broadx.tex)
and its incremental-temperature use; it is not Gaussian smoothing in energy.

Both exterior source laws are mandatory: zero cross section or an inverse-speed
continuation matched to the corresponding endpoint. They are declared physical
assumptions, not inferred evaluated data. Point queries may extend outside the
finite source interval only under those explicit laws. Heated output energies
must be positive because the cross section generally diverges at zero speed.
The prepared owner caches a global bound on g. Truncating `|x-y|` at L omits no
more than `g_max/sqrt(E) * [erfc(L)+2L exp(-L^2)/sqrt(pi)]` in cross section:
the excluded relative speeds require a target thermal speed above L. This
separate analytic kernel-tail bound does not bound errors in the exterior laws.

Integration partitions at every original histogram/linear segment and jump,
then into panels at most one thermal width long. Each visited linear interval
prepares its endpoint roots, thermal coordinates and interpolation law once;
Gauss nodes borrow that data without repeating endpoint roots or source-array
lookups. Deterministic Gauss 8/16 comparison controls adaptive refinement; its
reported error is a numerical estimate, not an interval certificate. The estimate plus omitted-kernel bound
must fit the caller's absolute/relative point tolerance. Original zero segments
need no quadrature. Source scanning, output allocation, segment visits, kernel
evaluations and refinement depth are bounded; output knots share one work
allowance. Callbacks poll cancellation throughout preparation and integration,
and failures publish no partial curve or mutation of the original source.

`tabulate` produces the original `PointwiseCrossSection` type with the actual
heated temperature, unchanged evaluation/nuclide/reaction/AWR/Q identities and
an explicit new data ID. Its result separately retains the source ID,
temperature, exterior laws and maximum knot error estimates. The **caller owns
the output energy grid**: linear interpolation of those knots adds a distinct
approximation that is not bounded by the quadrature estimates. Broadening moves
resonance wings away from the input knots, so reusing a coarse source grid is
not automatically adequate. Grid refinement remains necessary before claiming
accuracy of collapsed reaction rates.

The generated curve feeds both original flat-in-energy and explicit-spectrum
collapse APIs. Their temperature tables reach the existing irradiation,
capture, depletion and nuclear/thermal owners without a second reaction or heat
solver. This supplies no resonance reconstruction, self-shielding, transport,
bound-atom scattering law, probability table, secondary-particle kinematics or
evaluated-data accuracy certification.

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

Discrete and continuum radiation records retain their declared FD/FC
normalizations, relative intensities, energy uncertainties, transition types,
conversion coefficients and interpolation laws. Continuum densities convert
from probability/eV to probability/J together with their energy coordinates.
Beta endpoint records remain distinct from monoenergetic particle lines;
relative intensities are not normalized to one and already include branching.
Unknown spectrum mean energies remain unknown. Covariance payloads are
structurally traversed through SEND and their presence is retained, without
uncertainty propagation. No radiation transport or deposition is inferred.

Input bytes, all scanned records, retained modes, spectra, discrete records,
continuum points and interpolation regions have explicit ceilings checked before
data-dependent allocation. The latter three storage bounds are cumulative across
radiation types. Aggregate-energy storage remains bounded by seventeen pairs;
checkpoints preserve the caller's original error.

The reader regressions cover the pinned published ENDF/B-VIII.0 tritium,
helium-3 and cobalt-60 files, unknown energy components, spectral SI integrals,
photon multiplicity and conversion coefficients, malformed data and
identity/count/padding/SEND, exact resource and cancellation limits, and
unrepresentable input. The default-feature actual-root native nuclear gate
passed all 183 tests after spectrum retention. The fixture source and exact original Git blobs are
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
origin codes and spontaneous fission (which has no single residual) refuse in
this assembly. Multi-product spontaneous fission uses the sibling native burnup
assembly described below. Branch closure, topology and model work limits are
admitted by the original `DecayModel`, without branch normalization.

`with_burnup_product_retention(DecayProductRetention { alpha_particles, protons }, ...)`
optionally retains emitted He4/H1 in that same burnup inventory. Both switches
default to false. Enabled species require explicitly supplied stable ground-state
evaluations; absent products or nonstable laws refuse. The original RTYP pass
caches charged-particle multiplicities once. Lowering merges coincident identities,
so Be8 alpha decay yields two He4 nuclei, and repeated emission codes retain each
charged nucleus. The residual-lineage `DecayModel` and parent half-lives stay
unchanged. Generated helium/hydrogen can use the existing selective compartment
transport and extraction; their birth times and subsequent residence are resolved
within the same positive reaction solve. Spectral binding compares the exact same
retained products and multiplicities. Retaining more nuclei never multiplies a
parent event's heat or photon spectrum.

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
solver tolerances or treated as validation evidence. Charged-product retention
does not infer kinetic stopping, electron capture, material composition or mass,
neutron inventory, or a complete nuclear mass/charge ledger. Energy routing
remains an independently supplied assumption.
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

Four further regressions cover retained alpha inventories after one half-life,
coincident and repeated He4/H1 yields, selective helium movement/extraction during
its production interval, and exact spectral binding with one energy credit per
decay. The focused original-source native gate passed all 24 assembly/emission
tests. This gate compiled the full default-feature nuclear production source and
used the original selected numerical dependencies; it is not a workspace Cargo
or constellation integration result.

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

## Energy-resolved emission from integrated decay histories

`endf::decay::model::emissions::PreparedDecayEmission` prepares one selected
radiation type on an explicit energy grid. It retains per-parent particle and
first-energy-moment coefficients, the original evaluated assembly and spectral
coverage. The same map applies across interval durations and histories without
reintegrating the spectral records. Application scales by original integrated
parent transformations, including daughters produced and decayed within the
interval; endpoint activity is not used as an interval rectangle.

Absolute intensities use FD times RI for discrete records and FC times RP(E)
for continua. Parent-wide branching is already included and is not applied
again. Histogram, linear and logarithmic continua (ENDF INT=1 through 5) use
positive analytical count and energy integrals over their intersections with
the group grid. Segment logarithms are reused across group intersections;
small-contrast limits and separately scaled energy factors avoid cancellation
and unnecessary intermediate overflow. Logarithmic domains require positive
energy and/or density as specified by the supplied law. INT=6 remains outside
this spectral-density interpretation. No midpoint energy or
silent renormalization replaces the supplied distribution. Internal grid edges
belong to the upper group; the last upper edge is inclusive. Positive spectral
support outside the grid refuses instead of silently dropping radiation.

Admitted monoenergetic lines include gamma/X-ray photons, alpha particles,
neutrons, protons and discrete electrons. Beta endpoints do not become emitted
line energies. Supplied beta-minus continua can be grouped, while unresolved
electron-capture/positron intensity, fragment/unknown radiation semantics and
unsupported positive interpolation laws refuse. Covariance propagation, beta
shape reconstruction, angular transport, material attenuation and dose remain
separate work. A missing spectral record is explicit `Unreported` coverage;
the result lists any actual transformations with missing radiation data.

`bind_burnup` checks nuclide identities/laws and the original evaluated decay
branches, source labels, products and energy routing. Additional capture/fission
channels and reordered reactions are allowed. The bound response consumes only
physical decay reaction counts, including each compartment's original local
report from `BurnupTransportStep`. Movement does not emit particles or relocate
past emissions. Both entry points require the producing owner's original model
Arc; source reports themselves are caller-supplied values, not authenticated
receipts. These source observables never add heat to the original average-energy
calorimetry or change material mass and caloric laws.

Preparation and application retain bounded storage/work and cooperative
cancellation. Scaled positive accumulation avoids intermediate product overflow.
Unrepresentable positive per-decay coefficients refuse; amplification of a
subnormal prepared coefficient by more than one decay also refuses because its
lost preparation precision cannot be recovered by later multiplication. These
are numerical admission limits, not nuclear-data uncertainty bounds.

G1/G3 regressions cover analytical logarithmic moments, the inverse-energy and
inverse-square limits, neighboring coordinates, flat spectra, mixed regions,
discontinuities, group refinement and energy-unit rescaling by 1e-180 and 1e180.

The `evaluated_decay_emissions` example imports the original cobalt-60/nickel-60
records, evolves two physical half-lives and emits a two-group gamma history.
It independently checks the quarter-parent endpoint and the six-line evaluated
multiplicity of 1.99848902 photons per decay. `fs_couple::radioactive` re-exports
the original evaluated-data and emission owners under `nuclear-thermal`.

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

## Moving fuel and products with spatial neutron feedback

`DiffusionBurnupState::new_with_transport` binds the original
`BurnupTransportState` to exactly matching spatial volumes and compartment
order, retaining complete neutron/reaction bindings. One private inventory
owner holds the populations and accepted history. Every Picard iterate
restarts the original neutron and material histories, then invokes the existing
joint reaction/transfer/outflow/production operator once. Fuel or products
arriving in a compartment can react there during the same interval.

The coupling obtains mean parent populations from actual integrated residence,
including continuously populated zero-removal receivers. Movement counts,
nuclear reaction counts and external supply remain separate. Event heat stays
in its event compartment; previously emitted energy does not move with later
product transfer. Optional movement, residence and cohort-balance accounts are
returned by `DiffusionBurnupStep`, and the accepted transport state retains
its original preparation and lifetime movement accounts.

The independent-cell API and existing coupling budget remain available. Bound
transport path/coefficient ceilings stay fixed; the current burnup propagation,
work and accounting allowances apply to the entire joint inventory graph.
Native thermal binding and idle admission check that graph. Adaptive comparisons
include individual material-transfer and outflow counts once per shared cohort.
Analytical same-interval fission, return-product heat, zero-flow equivalence,
continuous supply/outflow, nonlinear ODE refinement and transactional refusal
regressions accompany native thermal and adaptive integration checks.

The joint inventory remains limited to 64 nuclide nodes across compartments,
with a separate existing neutron/precursor limit. Transport rates and geometry
are explicitly fixed. Mean-flux/mean-inventory time correlation still requires
timestep refinement. This addition does not solve carrier-fluid dynamics or
bulk material/enthalpy advection, and does not establish a certified global
error bound or evaluated-reactor qualification.

## Nonlocal radiation energy responses

`radiation::RadiationResponse` applies caller-supplied absorbed-energy fractions
to integrated outward-radiation joules. Each source row names explicit receivers
and an explicit escape fraction. Source and receiver IDs are opaque physical
identities, independent of the optional thermal crate. Immutable rows are shared
through `Arc`; canonical sparse ordering provides deterministic application
without a dense all-cell response matrix.

Fractions are finite and nonnegative. Each row must close to one within a fixed
16-machine-epsilon admission allowance; coefficients are never renormalized.
Any admitted row-closure defect is retained and charged in joules at the actual
emission magnitude. Duplicate source rows or receiver paths, missing/extra
emission coverage, work-limit violations and unrepresentable transfers refuse.
Exact zero emission and complete escape remain valid physical limits.

The result separates per-receiver deposition, total emission, deposition,
escape, signed residual and noncancelling numerical defects. Product rounding,
row closure, receiver reductions and scalar exports cannot cancel one another
to satisfy the caller's energy allowance. Positive products below representable
range refuse rather than disappearing. Checkpoints bracket work, each sparse
row/path and final publication. The response has no mutable physical history.

These coefficients describe absorbed ENERGY, not photon extinction or collision
probabilities. The distinction between production, local heating and separately
transported photon deposition follows the public
[OpenMC energy-deposition documentation](https://docs.openmc.org/en/stable/methods/energy_deposition.html).
In particular, attenuation is not automatically absorption; the underlying
quantities are distinguished in [NIST's photon energy-absorption description](https://physics.nist.gov/PhysRefData/XrayMassCoef/chap3.html).
No evaluated response data, radiation spectrum, transport solution, flight-time
delay, dose or data-uncertainty certificate is inferred. The response must be
valid for the complete emitted energy mix, material state and geometry to which
the caller applies it. Nuclear emission histories remain with their source owners.

Eleven focused G0/G3/G4/G5 regressions exercise unequal partitions, shared
receivers, deterministic ordering, complete escape/deposition, coefficient and
coverage admission, bounded work, small coefficient defects amplified by large
emission, scalar-export loss, subnormal arithmetic and cancellation.
All eleven passed in a standalone build of the actual response module using
the pinned `nightly-2026-07-06` compiler. This validates the response equations
and numerical boundaries, not external transport data or a workspace Cargo gate.

## Radiation deposition on the shared material thermal timeline

`NuclearHeating::with_radiation_response` binds the response to the complete
original source collection and explicitly named `RadiationReceiver` materials.
Response IDs use the original `CellId.0`. Receivers need no nuclear inventory;
multiple emitters can heat one receiver, and a source can reabsorb its own
outward radiation. Exact source/destination coverage is required. Receiver mass,
geometry, body/segment identity, material, available caloric identity and storage
basis retain the same correspondence checks as the original source bindings.
Changed correspondence refuses, including on an idle call.

Each final native interval first trials the original nuclear owners. The
response then partitions their integrated emitted joules. Local deposition and
all radiation contributions are compensated together per physical receiver,
then converted once to mean watts over that same interval. Authored boundary
heat remains additive. Only accepted thermal steps publish nuclear histories;
response, material-domain, energy-budget or checkpoint refusal discards the
entire candidate. The caller's exact thermal checkpoint error is preserved.
Source shutdown retains radiation released later by existing decay reservoirs.

The source `EnergyChannels` and lifetime emission accounts remain unchanged.
`FissionHeatingAudit` additionally reports `radiation_deposited_j`,
`radiation_escaped_j` and `radiation_roundoff_j`. Thermal credit now includes
local deposition plus the received radiation. Nonrecoverable source export is
never deposited. Without a response, all emitted radiation remains exported and
the original direct heat-credit path is retained without extra response maps.
Response, aggregation and heat-credit defects share the original whole-advance
numerical energy allowance; accepted substeps cannot renew that allowance.

All existing `NuclearHeatSource` families use this same response seam, including
prescribed fission, irradiation, burnup, capture activation, radioactive decay,
spatial diffusion and transported inventories. Adaptive execution compares all
original thermal cells, including radiation-only receivers. It reuses the same
response and accounting transaction on each trial. Immutable response/binding
data are shared across clones. The `fs-couple/nuclear-thermal` fission facade
re-exports `radiation`; the radioactive facade additionally exposes the shared
`NuclearHeating<DecayCellSource>` route and its receiver/audit/error types.

The fixed response assumes its full emitted-energy mix and material-domain
applicability remain valid. No changing spectrum, scattering solution,
flight-time storage, dose or automatic response refresh is inferred. A new
externally supplied response may be attached for the same still-valid cell
bindings at an accepted clock boundary without replaying earlier emission.
Spatial thermal transport remains the
original explicit uniform-cell approximation, with its original refinement and
material-domain requirements.

Seven focused integration regressions cover analytic prompt/delayed calorimetry,
complete escape and self-return, multiple scheduled emitters and restart,
receiver correspondence, rollback/cancellation/accounting refusal, original
radioactive decay and adaptive receiver cooling against an independent ODE.
All seven passed against the complete production `fs-nuclear` library compiled
with `thermal` enabled, using the original applicable numerical and thermal
production modules. The adaptive balance separately retains gross nuclear
input, receiver storage and convective export. The broader test-library build
was terminated by the environment (exit 137), so this is a focused public-API
result, not a completed full-suite or workspace Cargo gate.

The original IF97 water example now demonstrates explicit synthetic self-return,
nonlocal absorption and radiation escape during irradiation and shutdown. Its
CSV separates emitted, received, escaped and nonrecoverable energy; its closure
checks include received radiation once. The full water example and `fs-couple`
Cargo checks remain unexecuted because the constellation dependencies are absent.

## Adaptive energy grids for Doppler-heated reactions

`PreparedDopplerBroadening::tabulate_adaptive` constructs a native
`PointwiseCrossSection` over an explicit positive energy interval. Its
`endf::doppler::adaptive::DopplerGridAccuracy` controls interpolation separately
from `DopplerAccuracy`, which still controls each original kernel integration.
`AdaptiveDopplerTabulation` retains the existing tabulation's source identity,
temperature, exterior laws and numerical-error reports, plus a sampled
interpolation-error estimate and the count of actual point evaluations.

Three interior probes compare the heated law with each candidate interval's
linear-in-energy interpolant. The comparison includes both the probe's
quadrature/tail uncertainty and the interpolated endpoint uncertainty. Cached
endpoints and interior values are reused when intervals split. Wide energy
ratios split geometrically to resolve inverse-speed tails without repeatedly
bisecting a large linear energy range.

Sampling alone could miss an off-center narrow resonance. During the existing
bounded source scan, preparation caches potential feature energies. It omits
only strictly interior knots whose two positive-width adjacent intervals have
exactly the same constant cross section. Both admitted histogram and linear
laws are then unchanged across that knot. Endpoints, duplicate-energy jumps,
and every nonzero ordinate change, including a one-ULP change, remain.
The original source data and every point-quadrature interval remain intact.

A bounded binary search checks this shared immutable feature cache against the
cell's retained thermal support. Where a feature is present, the accepted cell
must span at most the configured amount of `sqrt(E/theta)` (default 0.5,
maximum 1). Rationalized square-root differences retain narrow thermal windows
at high incident energy. Cached features are not copied into output grids or
rescanned at each probe; cache allocation is bounded by `max_source_points`.
Redundant plateau knots no longer force thermally fine output or exhaust an
otherwise adequate output budget. Identically zero curves retain two endpoints.

The interpolation estimate is **not an interval certificate**, a bound at every
unprobed energy, or a certification of nuclear data or arbitrary-spectrum
reaction rates. Integration error and grid error remain distinct. Both point
and grid absolute tolerances should admit physically negligible or vanishing
wings; a relative tolerance alone cannot resolve an exactly zero ordinate with
a nonzero omitted-kernel bound. Tightened-grid convergence remains available
through the original flat and spectrum-weighted collapse APIs.

Source-knot searches, cell refinement and every point evaluation share one
`DopplerLimits` work allowance. Output allocation is bounded by
`max_output_points`, refinement depth is at most 64, and unrepresentable
subdivision is refused. Cancellation and exhaustion return no partial curve
and leave the source unchanged. At unchanged temperature, the requested domain
must match the complete original domain, whose knots, interpolation laws and
jumps are retained exactly without kernel work.

Public-path regressions parse manufactured ENDF sources and exercise
off-center resonance discovery, inverse-speed tails, refinement of native
flat/spectrum-weighted reaction rates, exact identity behavior, shared work,
bounded output, cancellation and deterministic retry. Plateau-mesh regressions
admit the same flat law with 3, 17 and 257 input knots within a two-output-knot
budget for both histogram and linear interpolation. A separate G3 control
retains the thermal guard for one-ULP features and duplicate-energy jumps.
These introduce no evaluated-isotope accuracy, resonance reconstruction or
shielding claim.

## Energy-group radiation deposition

`radiation::grouped::GroupedRadiationResponse` binds a radiation kind and an
explicit ascending SI energy grid to independently closed sparse response rows
for each physical source and group. Changing isotope mixtures can therefore
change receiver heating even when the total emitted energy is unchanged.
Coefficients are supplied absorbed-energy fractions for the declared geometry,
materials and within-group spectrum; attenuation, thermal emissivity and
particle counts do not supply them automatically.

Canonical physical-source/group pairs reuse `RadiationResponse` for all
receiver reduction, coefficient closure, positive-product admission, numerical
energy accounting and cancellation. Physical identifiers remain independent of
internal row indices. Source coverage and group lengths are exact, including
zero emitters. Preparation checks explicit source, group, source/group, receiver
and transfer limits before the corresponding expansion. Application neither
regroups nor rescales supplied integrated joules.

G1/G3/G4 regressions cover changing spectral mixtures, split histories, sparse
multiple sources, self-return, receiver overlap, deterministic ordering, exact
coverage, per-group closure, allocation bounds, numerical refusal and early/late
cancellation. Spectral completeness and reconciliation with a nuclear source's
original outward-energy account remain obligations of the coupling owner.

## Evaluated fission-energy components

`endf::fission_energy::FissionEnergyRelease` reads identified ENDF-6 MF1/MT458
records through the original bounded reader. All nine components, native
polynomial coefficient/uncertainty pairs, reference energies and component
TAB1 overrides retain their source identities. Histogram/linear overrides
preserve LDRV and IFC; duplicate components, malformed shapes, unsupported
laws, inconsistent identities, excess work and cancellation refuse.

Incident-energy queries return J/fission. Signed polynomials use the original
electronvolt coefficient convention; no historical unit-repair heuristic or
source-total normalization is applied. NPLY=0 without a component table is a
thermal/threshold reference, and cannot silently become a constant law at all
incident energies. Reference access remains explicit. Outgoing fragment,
neutron and photon energies, delayed components and reported ER/ET totals stay
distinct; importing these records does not add heat or invent decay lifetimes.

Nine focused regressions passed in the pinned native build of the actual
`fs-nuclear` crate root using the original applicable numerical owners. They
cover format, component meaning, SI conversion, interpolation, reference-only
refusal, numerical range, resource limits and cancellation. This is not a full
workspace Cargo gate or independent validation of an evaluated isotope.
Format and physical distinctions follow [ENDF-102 (2023), section 1.5](https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf).

## Reaction-weighted fission energy and original burnup assembly

`endf::fission_energy::model::PreparedFissionEnergy` prepares all nine original
components over an explicit incident-energy domain. Reference-only components
require an explicit `ReferenceEnergyTreatment`: refusal, a caller-declared
constant approximation, or the ENDF-102 (2023) Sher--Beck prescription. The last
requires a same-evaluation prompt/total neutron-yield law and an independently
supplied nu(0); no missing zero-energy value is extrapolated. Component tables
take precedence. The 2023 neutrino correction is -0.100 E. Source ER and ET
remain independently observed totals, without normalization to component sums.

Histogram/linear tables, linear polynomials and piecewise Sher--Beck laws retain
their original interpolation and jumps. Higher polynomials use bounded adaptive
linearization with the h^2 max|q''| / 8 interpolation bound; signed coefficients,
nonnegative physical values, work limits and subdivision limits are enforced.
The reported per-component bound concerns the interpolation in exact arithmetic,
not evaluated-data uncertainty or a certified floating-point error bound.
Preparation and later condensation are separate, allowing immutable curves to
be reused. All energy ordinates passed to the original dimensionless integrator
are explicitly energy / 1 eV; the public physical results restore joules.

`PreparedFissionEnergy::collapse` reuses `collapse_reaction_weighted` to retain
integral(sigma_f phi E_component) / integral(sigma_f phi) for each component and
a separate incident-energy mean. No representative group energy substitutes for
that integral. Zero-fission groups have undefined means (`None`). The uniform
interpolation bound also bounds its contribution to a positive weighted mean.
Group flux amplitudes may subsequently vary; within-group spectral shape and
the processed cross-section state must remain applicable.

`CollapsedFissionEnergy::prompt_products` routes fragment, prompt-neutron and
prompt-photon energy to explicit deposited, radiation or escaped destinations.
It separately reports outgoing product energy, incident energy, signed prompt Q,
interpolation error and an outward-rounded numerical routing account. Outgoing
kinetic energy includes incident-neutron kinetic input. A coupled neutron-energy
owner must retire that input and avoid crediting transported outgoing neutrons
twice; a prescribed population-only flux has an external kinetic-energy supply.

`append_prompt_products_to_burnup` assembles these coefficients through the
original `FissionYieldCollapse::append_to_burnup`. Delayed neutron, photon, beta
and neutrino components are not released at the fission event. Product decay,
shutdown heat, inventory histories and thermal coupling retain their original
owners; MT458 cannot supply missing decay lifetimes. Both condensations must
match full original pointwise cross-section/spectrum data, interpolation laws,
processed metadata, grids and numerical group results. Matching labels and
group averages alone do not establish matching spectra. Shared original source
arrays use a pointer fast path; independent arrays are checked in cancellable
batches. A caller-specified routing allowance is in J/fission, separate from the
original runtime accounting budgets and the interpolation approximation.

The native pinned-rustc build of the actual complete nuclear core test root
compiled 281 owning tests without diagnostics. The focused run passed all 23
new fission-energy tests (nine importer, eight preparation/condensation, six
burnup integration), plus six affected original fission-yield tests. Coverage
includes analytic weighted integrals, discontinuities, nonlinear refinement,
source-shape mismatch despite equal group averages, irradiation/shutdown heat
counted once, bounded work and cancellation. Original applicable numerical
implementations were used; unavailable constellation dependencies prevent a full
workspace Cargo/DSR result. This increment does not establish thermal-feature,
facade/example compilation or empirical evaluated-isotope accuracy.

## Evaluated decay radiation energy partitions

`radiation::decay::PreparedDecayRadiation` prepares composition-dependent group
energies from an original `PreparedDecayEmission` and its original decay-energy
routing. Every parent with outward energy must have the selected spectrum and
no positive outward energy in another decay-heat category. The selected
spectral energy must agree with the original outward mean energy within an
explicit `DecayRadiationClosure` data-discrepancy allowance. Accepted group
weights are reconciled to that original energy; the absolute adjustment is
reported separately and never creates heat or enlarges a numerical tolerance.

Preparation retains immutable per-parent coefficients and model identity.
Application uses integrated original decay counts, including newborn daughters,
or matched original burnup decay-channel counts. Capture and fission channels
retain their independent energy accounts. Near-one branch sums are handled
explicitly when converting branch counts to parent events. Scaled positive
products, storage/work limits and cancellation reuse the original numerical
owners; missing spectra, incompatible models and unrepresentable amplification
are refused.

Stochastic application uses the actual integer parent and branch counts,
checking exact count consistency and the original 2^53 event ceiling before
binary64 conversion. It preserves the trajectory and consumes no randomness.
The returned spectrum is expected emission conditional on those decays; photon
sampling, angular correlations and branch-resolved cascades are not supplied.
Ordinary caller-provided step reports are not authenticated receipts.

Seven regressions cover daughter ingrowth and subdivision, different per-parent
routing, explicit spectral reconciliation, missing-energy refusal, reordered
burnup with competing capture, bounded work and cancellation, and stochastic
trajectory/replay/count consistency. Spectral reconciliation records data
agreement; it does not propagate evaluated covariance or establish data accuracy.

## Composition-dependent radioactive heating

`NuclearHeating::with_grouped_decay_radiation` attaches original evaluated
partitions and a `GroupedRadiationResponse` to the shared thermal transaction.
Source identities, original model identity, radiation kind and group edges must
agree exactly. Receiver binding retains original material, caloric law, mass and
geometry. Direct deposition and received radiation are combined once per cell
before conversion to thermal power. Receiver-only cells and self-return are
supported.

The source owner evolves each trial once. The projection reads its integrated
event report, and the existing thermal owner commits inventories, heat, clocks
and audits together. Failure or cancellation publishes none of those trial
histories. Immutable bindings are shared across trial clones. An aggregate or
grouped response can replace the previous response only at an accepted boundary.

Adapters cover deterministic decay and continuous production, capture
activation, burnup and production, migrating burnup compartments, diffusion
burnup, fusion-neutron activation, compatible kinetic source wrappers, and
stochastic decay, production, compartment transport and supplied transport. They reuse the respective original report and state
owners rather than replaying inventory evolution. Prompt capture/fission or
other outward energy without the selected radioactive response remains escape
and is identified by `radiation_ungrouped_j`. This is a subset of escaped energy.
`radiation_reconciled_spectral_j` records the evaluated mean/spectrum adjustment
separately from heat and arithmetic error. Actual group/account closure defects
remain subject to the existing advance-wide numerical energy budget.

Eight integration regressions cover independent parent/daughter calorimetry,
whole versus split histories, continuously born daughters, exact rollback and
cancellation, binding refusal, stochastic replay and thermal credit, competing
capture with separate prompt escape, and moving-parent/daughter residence with
nonlocal heating. Together with the pure partition and grouped-response tests,
all 438 tests passed in a native build of the actual nuclear library with
`thermal` enabled and the original applicable dependency implementations.
The pinned constellation siblings are absent, so this is not a completed
workspace Cargo, distributed-runner or full dependency-suite gate.

Energy-group response fractions remain independently supplied absorbed-energy
responses for their declared geometry, material and within-group spectrum.
Composition can change the mixture across those groups; changing material or
geometry still requires an applicable response. This does not infer attenuation,
solve radiation transport, compute dose or automatically update response data.

## Evaluated fission products across material temperatures

`CrossSectionTable::temperature_component` exposes one knot/group contribution
of an existing processed table. It retains that coefficient at its own knot and
zero at all other knots/groups, with the original linear-Kelvin interpolation
and temperature domain. All components share the original temperature and
coefficient arrays: a K-by-G table retains K*G numerical coefficients even when
all K*G components are used. Each component's flux hazard visits only its one
incident group. Nested decomposition, invalid indices, extrapolation and
unrepresentable positive coefficients refuse through the original owner.

`endf::fission_yields::temperature::append_temperature_series_to_burnup` builds
these partial `TabulatedFission` channels from identified processed-temperature
condensations. Each channel retains its knot's independent product yields and
prompt energy. Consequently the summed rates interpolate sigma, sigma*Y and
sigma*Q together; conditional yields and energy-per-event means are not
independently interpolated. Subsequent daughter decay retains its original
inventory, transport and thermal ownership.

The original MF8 evaluation, processed isotope/library revision, named group
grid and complete incident spectrum must agree across strictly increasing
knots. Full source content is retained and checked, with shared-array fast
paths. The explicit burnup temperature domain must lie inside the supplied
range; a single knot remains valid only at that exact temperature. Zero-reaction
knots remain interpolation anchors. Group flux amplitudes may change later,
but within-group spectral feedback and temperature-grid accuracy are not
inferred. All original reactions are preserved. The existing per-knot product,
export and duplicate-fission admission is reused, with global reaction/product
limits over the complete series and one shared cross-section payload budget.

Thirteen cross-section tests and thirteen evaluated-product tests passed in the
pinned native build of the actual nuclear core. The new cases cover partition
closure, independent production/energy moments, shared storage, one-group hazard
work, positive numerical range, zero anchors, incompatible source content,
aggregate limits and cancellation. This is numerical/assembly evidence for
manufactured sources, not empirical isotope validation or a workspace Cargo gate.

## Prepared planar photon transport

`radiation::transport::prepare_slab_response` calculates a native
`GroupedRadiationResponse` from ordered homogeneous layers and independently
supplied gamma/X-ray collision kernels. Analytic characteristic integration
retains uniform-source self-absorption, interlayer attenuation and two-sided
vacuum escape. A positive angular quadrature explicitly represents isotropic
emission. Same-group isotropic scattering is resolved by the original
`fs_la::prepared_lu::LuWorkspace`, factored once per interacting energy group
and reused across all physical sources. Groups without same-group scattering
bypass factorization. Downscatter is solved in descending energy order.

The transported coefficients describe ENERGY: local heat plus outgoing photon
energy must close the total collision coefficient. Photon-number cross
sections, total attenuation or a mass energy-absorption coefficient alone do not
supply this kernel. `PhotonMassMaterial` retains explicit material identity,
radiation kind, group grid and per-mass coefficients; conversion uses a supplied
original reference density and refuses lost positive products. No material
mixture, secondary-electron stopping law or evaluated data is invented.

Preparation checks resource bounds before expansion, retains the original
linear equations for absolute residual checks, and refuses nonclosing or
negative solutions without clipping or normalizing them. Admitted floating-point
overshoot of a generated unit fraction is bounded to one with its complete
correction charged, together with observed preparation closure/solve defects,
through the original response's numerical energy audit on every application. The sparse response
retains every optical layer, including zero-heating receivers. Arithmetic
closure, solve residual and charged work are reported separately. Subsequent
decay intervals reuse the immutable response and the existing energy audit.

This is steady infinite-planar transport with flat isotropic scattering
sources inside each layer. Ballistic characteristic integration is analytic
for each supplied ordinate; spatial scattering, angular and energy-group
discretizations still require refinement. There is no flight-time storage,
anisotropic phase function, upscatter, electron transport, dose or automatic
material-property refresh. Unrepresentable positive transmission tails refuse.
The formulation follows the collision-probability equations in SCALE/CENTRM
section 7.4.3.3; its use of photon energy kernels is explicit here.

G1/G3/G4 regressions cover analytic uniform-source absorption, multilayer
shielding, elastic-scattering zero heat, lower-group energy redistribution,
ballistic subdivision, convergence toward an independently derived two-stream
slab solution, density/depth scaling, numerical/resource refusal and
cancellation. Manufactured kernels establish these equations, not accuracy for
an evaluated shielding material.

### Composition-dependent photon opacity

`PhotonMassMaterial::from_mass_fractions` prepares a homogeneous mixture from
original constituent tables and explicit mass fractions. Total extinction,
local heat and every outgoing photon-energy coefficient are combined by
compensated mass-weighted summation. The result retains the target material ID,
each original source table, authored fractions and an explicit
`kernel_applicability` statement. Fractions must sum to one within the original
16-epsilon closure allowance; no fraction or energy coefficient is normalized.
All sources must share the photon kind and exact energy grid. Storage/work
admission precedes expansion, positive product underflow/overflow refuses, and
bounded checkpoints retain cancellation. Repeated transport uses the prepared
table without constituent loops. `macroscopic_cell` applies the actual mixture
density, so existing shielding and nuclear heating adapters consume the result.

NIST's [mass attenuation additivity](https://physics.nist.gov/PhysRefData/XrayMassCoef/chap2.html)
does not establish additivity of
[mass energy-absorption coefficients](https://physics.nist.gov/PhysRefData/XrayMassCoef/chap3.html).
The caller must justify independent interactions, common within-group spectral
weighting, and each supplied local-deposition/redistribution kernel IN the
mixture. This constructor does not derive those kernels from attenuation,
model secondary-electron transport, infer compound absorption-edge changes, or
qualify measured data. Resolved stratified materials remain distinct slab cells.

Focused G1/G3/G4 tests cover analytic mass-weighted multigroup coefficients,
thin/thick oblique Beer–Lambert calorimetry, an independent infinite
scattering/reabsorption series, partial-density and regrouping equivalence,
retained source identity, numerical range, incompatible inputs, budgets and
cancellation. Their sources are manufactured laws, not empirical validation.

## External collimated and diffuse photon absorption

`radiation::transport::prepare_slab_response_with_sources` admits explicit
uniform layer emission and collimated or diffuse gamma/X-ray irradiation through either
exterior face. Sources have independent identities; every physical layer remains
a heat receiver. A beam's supplied `mu` is its positive inward direction cosine,
so its path through a layer is `thickness / mu`. The grouped response accepts
incident joules actually crossing that face, with no inferred area or extra
cosine weighting.

`PhotonSourceGeometry::Diffuse` supplies constant incident radiance over the
inward hemisphere. Each ordinate carries the fraction
`mu * weight / sum(mu * weight)` of the supplied face-crossing joules. Normalizing
this angular moment preserves source energy even for a mirrored full-range
quadrature whose positive half does not integrate `mu` exactly. It does not
renormalize material coefficients or transport responses. The resulting angular
approximation must be refined separately from layer subdivision. G1/G3 tests
check convergence to the analytic absorption-only transmittance `2 E_3(tau)`,
vacuum and thin absorption, orientation reversal, multigroup scattering energy
closure, and diffuse heating on the original native thermal cells.

Analytic Beer–Lambert first collisions partition the beam into direct local
heat, outgoing photons and uncollided transmission. Outgoing energy, including
same-group scattering and downscatter, enters the original flat-source
multigroup operator; its factorization is shared with all internal sources.
Escape includes both transmitted and scattered energy leaving either face.
Total attenuation is never treated as absorbed heat. Beam closure and kernel
roundoff are retained in the original per-joule deposition audit; work,
storage, cancellation and positive-tail refusal remain enforced.

Isotropic transfer operators are prepared only for groups with internal
emission or an incoming scattering channel. Purely absorbing boundary beams
therefore use linear layer work and do not fail on unilluminated grazing
directions. Actual diffuse, internal and scattered paths retain the original
angular law and positive-tail checks; no physical scattering is omitted.
Focused regressions compare a thick collimated path with Beer--Lambert, retain
the refusal for actual grazing illumination, and check the deterministic work
account under layer doubling. The existing conservative storage admission is
unchanged; the work regression is not a wall-clock performance claim.

G1/G3/G4 tests cover normal/oblique thin/thick absorption, ordered multilayer
calorimetry and reversal, an independent infinite collision-series solution,
pure scattering with zero heat, lower-group absorption, simultaneous internal
and external sources, and incidence/identity/resource/cancellation refusal.
Scattered emission retains the isotropic flat-source approximation and requires
spatial/angular convergence. Automatic spatial refinement presently accepts
uniform internal sources through its original entry point. The explicit-source
`refine_slab_response_with_sources` also supports boundary irradiation and mixed
source families. Internal sources divide among sublayers; boundary sources retain
their full incident joules, physical side and angular distribution on every grid. Receiver
energies reduce to original physical cells. The same work/storage allowances,
consecutive-grid convergence rule and charged arithmetic audit remain active.
G1/G3/G4 controls add beam convergence against an independent two-stream
solution. Surface reflection, flight-time storage and dose remain unmodeled.

## Evaluated prompt fission energy on the changing material timeline

`endf::fission_energy::temperature::append_prompt_temperature_series_to_burnup`
accepts matched product/MT458 condensations at each processed temperature and
uses the same temperature-component assembly. Explicit fragment, prompt-neutron
and prompt-photon routing is retained at every knot. The original fixed-state
bridge and the series share source matching, routing and the J/fission routing
allowance; neither path adds delayed MT458 components to product-decay heat.

`CollapsedFissionEnergy::source_preparation` retains the actual prepared law.
A series requires identical original MT458 coefficients, uncertainties and
tables, the same continuation policy and Sher--Beck neutron-yield source where
applicable, and the same prepared curves, domain and approximation bounds.
Changing a reference policy or numerical preparation at a temperature knot
cannot masquerade as temperature dependence. Different group means resulting
from the processed fission cross section remain valid.

The result is ordinary `BurnupData` with original `TabulatedFission` channels.
Existing finite inventories, prescribed or kinetic fluxes, transport and
material heating consume it through their established owners. There is no new
inventory, thermal integrator or delayed-energy reservoir. Only actual native
thermal substeps publish accepted nuclear histories; out-of-domain temperatures
and cancellation discard the candidate. The original frozen-temperature step
remains a first-order partitioned approximation, with native/adaptive timestep
refinement available. Processed temperature knots must be supplied explicitly;
Doppler processing, self-shielding and within-group spectral feedback are not
automatically inferred by this assembly.

The actual nuclear root with `thermal` enabled compiled 478 owning tests without
compiler diagnostics. All 28 evaluated-energy tests, three new evaluated-fission
thermal integration tests and three existing thermal-temperature tests passed.
Together with the 26 cross-section/product tests for the temperature-component
assembly, the focused scope contains 60 distinct passing tests, including 19
new regressions. The new thermal cases exercise self-heating against an
independent continuous ODE with first-order timestep refinement, shutdown and
restart with daughter heat counted once, and complete domain/cancellation
rollback. Native original numerical and thermal implementations were used;
missing constellation dependencies still prevent full workspace Cargo/DSR
verification. No empirical isotope or release-performance claim is introduced.

## Spatial convergence of photon deposition responses

`radiation::transport::refinement::refine_slab_response` doubles homogeneous
sublayers while retaining the original energy groups, angular rule and material
kernels. An original uniform source is divided by the same exact power-of-two
weights. Deposits are reduced back to original physical receiver IDs before the
response is exposed; original thermal cells and nuclear inventories are never
subdivided or copied into new physical owners.

Every source/group/receiver coefficient and the escape fraction must agree
between consecutive grids within explicit absolute-plus-relative tolerances.
Redistributing kernels require two consecutive agreeing comparisons; a failed
comparison resets the count. This avoids accepting the known symmetry plateau
between one homogeneous scattering layer and two equal half-layers. Ballistic
kernels need one comparison. The last-grid change and accepted subdivision
count remain visible. This is an observed spatial convergence test, not a
certified error bound or a test of angular/energy resolution. Failure to
converge, cancellation or exhausted cell/storage/work limits returns no coarse
fallback.

A single work allowance covers every grid solve, expansion, reduction and
comparison. The core's storage admission is reused before sublayer allocation.
Fine-grid linear solves can apply at most two residual corrections through
their existing LU factors. Compensated residuals retain signed FMA product
tails; the original absolute unit-energy residual threshold remains mandatory.
Original coefficient-closure and preparation defects, plus aggregation
arithmetic, survive the reduction and are charged by the existing response
audit in emitted joules. Three regressions verify analytic ballistic invariance,
improved agreement with an independent scattering ODE at tighter tolerance,
original-ID mapping, zero-heat receivers, retained numerical charges and refusal.

## Material-bound photon transport on the decay heat clock

`NuclearHeating::with_photon_slab` and
`NuclearHeating::with_refined_photon_slab` derive a planar photon problem from
the original material cells. The shared input builder checks original material
identities, radiation kind, energy grid, source coverage and ordered collinear,
contiguous geometry. Reference density is stored mass divided by original
volume; layer depth is volume divided by an explicitly declared common area.
The geometry tolerance must resolve each layer and cannot admit a reversed
physical segment. The core's dimension/storage admission precedes material
coefficient expansion.

The prepared response enters `with_grouped_decay_radiation`, retaining the
same nuclear reports, original caloric laws, receiver binding, accepted clocks,
heat credits and atomic rollback. All optical layers remain bound, including a
transparent or scattering-only layer receiving zero heat. Refinement changes
only the temporary photon calculation. Its final response addresses the
original material cells, so accepted decay heat is still credited once.

Both routes apply to the existing deterministic and stochastic decay,
production, activation, burnup, migration, diffusion and compatible kinetic
source adapters. The photon kernels are fixed material data with explicit
applicability; there is no inferred thermal emissivity, temperature-dependent
opacity, expansion or replacement mass/material law. The infinite-planar,
isotropic-scattering, steady-transport and local charged-particle deposition
assumptions remain explicit.

Six thermal regressions cover independent parent/daughter shielding and
receiver calorimetry, original-mass effects, zero-heat optical-layer binding,
geometry/material and exact callback refusal, refined original-ID heating,
nuclear-history preservation, whole/split histories and early resource refusal.
The radioactive
`fs-couple` facade exports the transport/material/refinement types. The runnable
`evaluated_decay_emissions` example now combines the original published Co-60
decay records with explicitly manufactured photon kernels and a spatially
refined response; its CSV separates direct decay heat, shield heat and photon
escape, and checks the original radiation account at every interval. It does
not claim that its illustrative shielding coefficients are evaluated data.

## Prescribed photon illumination on the native heat timeline

`thermal::ExternalPhotonHeating` drives the existing enthalpy timeline from
explicit gamma/X-ray group-power histories without requiring nuclear inventories
or source identities that name thermal cells. `PhotonPowerSchedule` retains its
authored identity, exact radiation kind and energy grid, and continuous
piecewise-affine absolute-time knots. Powers are already total energy per second
entering the declared source geometry; no extra area or cosine factor is applied.
The final knot is held, so a finite pulse ends with zero power.

`from_slab` and `from_refined_slab` reuse the original material/geometry admission
and photon transport owners. Original mass/volume supplies density, volume/area
supplies layer depth, and every response destination remains bound, including
transparent or scattering-only cells. Refinement changes the temporary photon
grid while retaining original source energies and thermal receiver identities.

The observer clips native intervals at source knots, integrates their powers
analytically, and applies the existing grouped response. Absorbed joules become
native mean heat rates once; incident, absorbed, escaped and credited energy
remain separate. Immutable schedules and responses are shared across trial
clones; accepted clock and accounts publish only after the original thermal step
and binding checks succeed. Other sources compose through `CoupledObservers`.

The default constructors retain fixed material coefficients; the temperature
feedback constructor below supplies explicit changing opacity. Transport remains
steady. Photon flight-time storage, thermal expansion and charged-particle
transport are not inferred. Exact source integration does not eliminate native thermal
timestep or photon spatial/angular/energy discretization errors. Focused kernel
and owning integration regressions cover affine energy, spectral pulse
calorimetry, independent identities, refinement, restart, zero-heat scattering,
and atomic refusal; their numerical checks do not qualify measured source data.

## Temperature-dependent photon absorption and shielding

`TemperaturePhotonMassMaterial` retains complete original photon mass kernels
at ordered temperature knots. Total collision, local heat and outgoing photon
energy coefficients interpolate together on one original material identity,
radiation kind and exact energy grid. Exact knots retain their source objects;
interior samples pass the original closure/downscatter admission without
renormalization. The closed kelvin domain never extrapolates.

`TemperaturePhotonMassMaterial::from_mass_fractions` prepares homogeneous
mixtures at constant authored mass fractions. It retains the original
temperature laws and their complete-kernel applicability in the target mixture.
The union of constituent knots on the intersection of their closed domains
preserves intermediate opacity peaks and changes of slope; mixing endpoints
alone is insufficient. Every constituent, including zero-weight sources, must
share the photon kind and exact energy grid. Prepared samples preserve local
heat and outgoing energy separately. Attenuation additivity does not establish
the required mixture-dependent deposition kernels. Retained source/sample
values and preparation work are bounded; subsequent sampling needs no new
constituent evaluations. Analytic transmission and native piecewise-opacity
heating regressions exercise both interpolation and the existing thermal path.

`ExternalPhotonHeating::from_temperature_slab` evaluates these kernels from
the original native caloric temperatures. Stored mass and geometry still
supply reference density and layer depth. Changed illuminated start states
rebuild the original multigroup slab response, optionally through its existing
spatial refinement. Unchanged temperatures and dark intervals reuse the
prepared response. `prepared_temperatures_k` identifies the response's actual
start-state temperatures.

The existing source schedules, absorbed/escaped/credited accounts, thermal
storage and atomic observer transaction remain authoritative. Pending response,
cache and work updates publish only with an accepted thermal endpoint.
Endpoints and events validate every table's supported temperature, including
final endpoints with no subsequent source query. Aggregate table storage,
preparation count and cumulative transport work have explicit bounds.

Coefficients are frozen over each native explicit interval. Thermal timestep
refinement remains necessary; accounting tolerance and spatial response
agreement do not certify temporal error. This extension does not infer
temperature-dependent density, phase mixtures, expansion, photon flight time,
electron transport or measured material qualification.

## Temperature feedback for internally emitted nuclear photons

`NuclearHeating::with_temperature_photon_slab` applies the same prepared
temperature-dependent material kernels and feedback cache to evaluated internal
photon emission. Original decay/activation/burnup histories supply integrated
group energy; local nuclear deposition and ungrouped outward energy keep their
existing owners. Every emitter remains uniform and isotropic in its original
slab layer. Changed illuminated start temperatures refresh the transport
response; unchanged temperatures and dark intervals reuse it.

The response, preparation account and temperature cache remain pending until
the original nuclear/thermal transaction accepts its endpoint. A domain,
thermal, transport or cancellation failure preserves accepted inventories,
response and cache together. All endpoints retain the original temperature
support and material/mass/geometry bindings. Fixed-response rebinding disables
temperature feedback. Native timestep refinement remains necessary for this
first-order coupling; the original spatial/angular transport limitations and
the exclusion of inferred density/composition changes remain unchanged.

`NuclearHeating::advance_adaptive_to` includes transported radiation in its
existing componentwise time-error comparison. Each receiver's absorbed joules
and the emitted/absorbed/escaped interval totals use the existing
`NuclearTimeTolerances::nuclear_energy_j` allowance. Temperature, stored energy,
nuclear inventories, reaction counts and original released-energy comparisons
remain active. Loose thermal scales therefore cannot conceal an inaccurate
radiation heat partition between receivers.

`last_radiation_deposition` retains the original response result from the last
accepted native interval, published with its nuclear/thermal transaction after
crediting radiation once. Direct nuclear heat and ungrouped outward energy stay
with their original accounts. The adaptive single-interval owner discovers
CFL/event cuts before source evolution; its coarse result is compared with the
sum of two actual half-interval results, never continuing cumulative histories.
These remain local temporal error estimates, separate from conservation
roundoff and spatial/angular transport accuracy.


## Temporal convergence of external photon absorption

`ExternalPhotonHeating::advance_refined` compares complete native trajectories
from one accepted input using at least three physical time partitions. Each next
maximum step halves the prior largest actual native interval, so existing CFL
or event clipping cannot leave the refinement axis inactive. Native conduction,
calorics, latent storage, source knots, photon transport, domains and accounting
gates remain the original owners; the returned result is the finest original
trajectory, never an extrapolated state.

Componentwise agreement includes specific stored energy and caloric temperature
of every native cell, each original receiver's accepted photon thermal credit,
and incident/absorbed/escaped/credited photon accounts. Relative scales use
changes from the input and this call's energy, not the arbitrary caloric datum
or old cumulative totals. The finest difference must meet tolerance; successive
comparisons must show decreasing resolved discrepancy, or both be
roundoff-limited. Observed order is reported only for resolved differences with
actual dyadic interval maxima. These are measured endpoint/history-energy
comparisons, not rigorous global bounds, pathwise extrema, spatial/angular/
spectral certification or measured-law qualification.

Substeps, response rebuilds and original transport work are bounded across all
candidate trajectories, including discarded work and nested spatial refinement.
Temporary remaining-work caps preserve the original per-advance limits and are
cleared on the returned source. Cancellation, physical-domain and numerical-
account failures refuse the entire call immediately. No failed trial publishes
thermal state, source accounts or cached response. Zero-duration calls preserve
native event semantics without inventing a refinement ladder.


## Radioactive neutron sources on the original spatial history

`diffusion::decay::DecayNeutronState` binds original radioactive inventories and
explicit parent-wide neutron group yields to existing spatial neutron cells.
One prepared reaction operator advances the complete block
`[A_neutron, Y diag(lambda); 0, A_decay]`, including daughter ingrowth, emission,
transport, capture, multiplication and leakage during the same interval.
Integrated decay events and neutron residence come from that same propagation;
neither endpoint activity nor uniform injection of an interval's total emission
replaces the actual birth history. Stable-nuclide yields must be zero, and
branching already included in a parent-wide yield is not applied again.

The existing decay owner retains source populations, branch events and raw
energy exports. The original diffusion owner retains reaction, interface,
fluence, precursor and leakage accounts. Their extensive defects are summed
without cancellation, including aggregate source allowances. Prepared maps are
shared only across identical intervals and propagation limits. Refusal or
cancellation publishes neither component's proposed state or cumulative history.

`diffusion::decay::energy::DecayNeutronEnergyState` transfers emitted neutron
kinetic energy out of the parent's original outward radiation allocation. It
attaches the existing neutron kinetic owner to the accepted spatial report;
there is no second neutron propagation. Stored neutron energy, deposited heat,
remaining outward radiation, escape and signed prompt reaction Q stay disjoint.
The raw source release is diagnostic and must not be added again to the combined
cell exports. Binding begins new combined export counters while preserving
supplied source histories and resident neutrons.

G0/G1/G3/G4 tests cover analytic parent and daughter emission with simultaneous
capture/leakage, true birth-time versus uniform-source behavior, split intervals
and restart, stable/zero-duration/source-free limits, disjoint energy balances,
invalid model/energy bindings, work limits and late cancellation. These prove
the admitted equations and accounting, not evaluated-isotope calibration.

The joint system retains the existing 64-node bound. Sources are stationary
closed radioactive chains; incident neutrons do not deplete these source
inventories, and target coefficients remain the prescribed diffusion background.
The original diffusion approximation and frozen-coefficient assumptions apply.
The kinetic energy layer admits explicitly funded prompt fission as described
below; abstract delayed-precursor energy models remain excluded. Parent-wide expected spectra do not infer
sampled branch/energy correlations, angular transport or nuclear-data uncertainty.

## Evaluated radioactive neutron emission

`PreparedDecayEmission::bind_neutron_group_means` binds original evaluated
neutron spectra to an explicit `NeutronKineticModel`, retaining both original
Arcs, projection identity, original-order yields and signed per-parent
count/energy residuals. Each source-bin mean projects onto neighboring authored
representative energies; parent-wide branching is not applied again. Both
moments must close within the fixed floating-point allowance, and emitted
kinetic energy must fit the original heavy-particle outward-radiation allocation.
Gamma radiation, deposited heat and escaped energy cannot fund neutron creation.

Missing unstable spectra require explicit nonemitting indices and cannot
contradict supplied positive neutron-emitting modes. Source storage/work limits
and cooperative cancellation apply; extrapolated means, incompatible species
and unrepresentable coefficients refuse. Bin-mean condensation does not resolve
within-bin reaction thresholds, angular transport or evaluated-data uncertainty.

Five manufactured G0/G1/G4 regressions cover parsed parent-wide delayed-neutron
emission with analytic birth/absorption history, continuum count/energy
condensation, explicit missing-spectrum declarations, incompatible species and
energy funding, bounded storage and late cancellation. These establish format
and coupling behavior, not evaluated-isotope accuracy.

## Shared material heating from radioactive neutron emission

`thermal::{DecayNeutronCellSource, DecayNeutronHeating}` binds a complete shared
radioactive-neutron checkpoint to distinct original thermal cells, requiring
matching reference volumes, material identities and accepted clocks. One joint
trial advances each final native thermal interval; every cell receives only
its own disjoint deposited energy. The accepted physical state and interval
report are shared across all cell views, including when thermal IDs order
cells differently from their spatial indices.

Radioactive transformations, neutron captures, source continuation, resident
neutron storage and numerical defects retain their original owners. Shared
energy/count defects are charged once. The original thermal transaction admits
source and material endpoints together; cancellation, numerical refusal or a
material-domain failure publishes neither nuclear nor thermal history.
Material calorics and physical mass remain the original thermal owners, and no
temperature-dependent half-life, cross section or deposition fraction is
inferred from a temperature change.

Five G1/G3/G4 integration tests cover analytic decay/capture calorimetry without
double credit, shared spatial cell mapping, restart, clock/volume/cohort
refusals and transaction rollback. A separate two-source/two-cell analytic
regression checks both directed neutron transfers, each cell's capture and
fluence, original inventories and cell-resolved heat/radiation/escape.
These checks retain the source, transport and spectral limits stated above;
they do not establish temperature-dependent nuclear feedback or a full
constellation/thermal-extension gate.

## Finite radioactive emitters under neutron irradiation

`diffusion::burnup::decay::BurnupDecayNeutronEmission` binds explicit parent-wide
neutron yields to the exact original `BurnupModel`. Enabling these emissions on
`DiffusionBurnupState` makes radioactive decay, capture, activation products and
neutron emission consume the same finite stock. Each existing Picard iteration
freezes the admitted microscopic-rate contraction and macroscopic neutron
coefficients, then advances one positive neutron/nuclide block. Continuous
neutron, precursor and nuclide feed, daughter birth and residence all retain
their ages within that block. Neutron-source edges use only spontaneous-channel
fractions; capture cannot create a fictitious radioactive neutron.

Acceptance checks both original neutron-versus-nuclear event agreement and
actual mean target population times actual integrated fluence. The latter
prevents a newly born target and neutron field from falsely converging at zero.
This closes the frozen-mean equations, not their unresolved nonlinear time
correlation; physical interval refinement remains necessary. Existing explicit
work, node, edge, iteration and cancellation limits apply to the entire joint
system. Independent inventories and the original finite-compartment transport owner
are admitted. The combined neutron/precursor/nuclide node ceiling remains 64.

The existing `FusionNeutronActivationState` accounts that same accepted neutron
trajectory. Emitted kinetic energy is an internal transfer funded exclusively
by original, rate-weighted decay radiation, never capture energy. Raw nuclear
release remains diagnostic; combined heat, remaining radiation, escape and
neutron storage are disjoint. External injection is reconstructed from authored
external rates, with per-cell/group reconciliation rather than subtraction of
large totals. Parent-count defects are converted by actual neutron yield before
being charged to a neutron allowance.

`FusionNeutronActivationCellSource::bind_without_external_source` uses the
existing shared thermal cohort and an explicit empty forcing history. Resident
neutrons and radioactive stock continue heating original materials after
external forcing ends. Material temperature feeds existing microscopic data;
original volumes, calorics, physical mass, clock and atomic publication remain
owned by the thermal transaction.

Thirteen new manufactured regressions exercise analytic birth ages, continuous
parent/daughter feed, capture of the emitting stock against an independent
nonlinear ODE with physical refinement, zero-target startup, rare emission,
energy/source separation, native material heating, restart and refusal. The
affected native test selection passed 114 tests, including the separately
bound evaluated-spectrum extension; whole-workspace and full thermal-extension
gates remain unclaimed.

## Evaluated spectra on finite emitting inventories

`EvaluatedDecayNeutrons::bind_burnup` retains the original evaluated emission,
kinetic projection and strict `BurnupDecayEmission` mapping while supplying
the finite stock's `BurnupDecayNeutronEmission`. Original nuclide order,
identities, half-lives, products, branch fractions, energy routes and neutron
grid must match. Reordered decay channels and appended neutron-induced
reactions are admitted; neither new isotope data nor a missing-spectrum zero
is inferred. Projected parent-wide yields retain their original count and
energy moments without a second branch factor.

Two additional native regressions pass for evaluated emission with competing
finite capture, reordered decay channels, changed identity/energy/grid refusal
and late cancellation. They test the admitted mapping and manufactured
equations rather than evaluated-isotope calibration.

## Transported finite radioactive emitters

The finite-emitter binding also admits the original `BurnupTransportState`.
Its physical channel compiler and accepted-evolution ledger are shared with
the standalone transport path. One neutron/nuclide operator includes decay,
irradiation, donor-proportional transfer, exterior outflow, daughter ingrowth
and continuous supply. Every nuclear and neutron-source edge is scaled to the
complete removal hazard; departure is never counted as a radioactive event.
Physical reaction and movement counts use the same joint residence integrals.

A daughter may be born in one cell, move and emit in another during the same
interval. Original per-cell reaction/heat histories stay at their event
locations; moving a remaining nucleus does not move previously released heat.
Exterior outflow retains its original nucleus ledger and excludes subsequent
emission from the represented domain. It creates neither neutron leakage nor
heat. External feed remains distinct from internal transfer and products.

Transported compartments retain the exact shared nuclear model and original
nuclide order, matched neutron-cell volumes and unchanged neutron yields for
each transported nuclide. An external joint solve clears incompatible
nuclear-only prepared maps; native standalone advances restore their own
compatible map and reuse it under the original rules. Transfer-path and
coefficient limits remain active alongside the complete joint node/edge bounds.

The existing activation energy and material transaction consume the resulting
local events directly. Seven additional native regressions cover analytic
moving-parent and daughter emission, outflow/restart, continuous-feed reduction,
nonlinear moving-emitter capture against an independent ODE with physical
refinement, joint-budget/cancellation rollback and local material calorimetry.
The affected suite passed 129 tests, including the original transport cache,
circulation, continuous-feed and lifetime-account regressions.

Compartment rates remain prescribed. Carrier fluid dynamics, transported bulk
mass and sensible/latent enthalpy are separate owners; this coupling does not
infer them. Frozen-mean reaction closure, diffusion-model limits and the
existing kinetic-energy model exclusions still apply.

## Evaluated spontaneous-fission yield input

`endf::fission_yields::spontaneous::SpontaneousFissionYieldEvaluation` admits
identified ENDF-6 NSUB=5, MF8/MT454 independent products through the existing
bounded record reader. The source must declare one energy-independent LE=0
set. Product identities, isomers, multiplicities, uncertainties and source
energy-field metadata remain unchanged; yields are never normalized to one or
two. Cumulative MT459 yields cannot substitute for immediate product births.

The opaque type cannot enter incident-neutron collapse. The original induced
yield API still requires NSUB=11. Fragment yields establish neither half-life,
branching, event energy nor neutron multiplicity/spectrum; these remain
separate explicit inputs. Existing allocation/record limits and cancellation
apply, including a final checkpoint before the evaluation is returned.

Three manufactured regressions cover data retention, cross-sublibrary and
cumulative/ambiguous/energy-dependent source refusal, explicit bounds and late
cancellation. All three passed in the recovered original-source native suite,
with the original induced-yield boundary retained.

## Prompt fission, finite products and decay energy on one trajectory

`NeutronKineticModel::new_with_reaction_energies` admits the original diffusion
owner's prompt fission channels. Supplied prompt yields retain their original
multiplicity and group distribution. The available immediate non-neutron energy
is E_in + Q - sum(yield * E_out); Q is prompt nuclear conversion before that
subtraction, not total delayed product energy. Active reactions must satisfy
the original threshold and nonnegative available-energy checks, including after
finite-composition or temperature rebinding. The convenience kinetic-only
constructor still requires explicit energy inputs for every neutron producer.

`FusionNeutronActivationState` uses the same accepted reaction trajectory and
original finite fuel/product inventories. Bound fission channels must carry
zero prompt event energy in burnup, even when explicitly declared Q is zero,
so prompt release cannot be credited twice. Actual radioactive products retain
their own half-lives, branches and later event energies. Declared product decay
neutrons can re-enter that same finite neutron/inventory block and consume only
the original radioactive radiation allowance. No abstract precursor energy,
missing neutron spectrum, fission Q or evaluated isotope property is inferred.

Original material-cell heating consumes disjoint prompt and decay deposition
through the existing shared thermal transaction. Nuclear identities, material
mass and calorics retain their owners; source-free continuation, restart and
late refusal preserve the established atomic state/history semantics. The
original diffusion approximation, representative group energies, fixed supplied
prompt Q and frozen-mean nonlinear closure still require physical refinement.

Four new kinetic/activation regressions cover analytic multiplication and
continuous feed, signed conversion with distinct outgoing groups, finite fuel
and product aging against an independent nonlinear ODE with step refinement,
and energy/precursor/cancellation refusal. One material regression follows an
initially empty neutron field through radioactive emission, prompt fission,
beta and delayed-neutron product decay, local heat, radiation and leakage.
It checks original material temperature/mass, restart and late rollback.
All five passed in the recovered original-source native suite with `thermal`
enabled. The material fixture also received an independent particle-balance
and caloric-domain check.

## Evaluated multi-product spontaneous decay in native burnup

`endf::decay::model::burnup::EvaluatedBurnupDecayAssembly` lowers original
MF8/MT457 modes directly into an original `BurnupModel`. Pure spontaneous
fission (RTYP=6) requires an identified `SpontaneousFissionProductBinding` with
independent NSUB=5 products for the exact parent/isomer. Each positive fragment
yield is retained in an explicitly supplied nuclide evaluation or named in the
explicit exported-product boundary. No missing evaluation becomes a stable
sink, no fragments are pruned or renormalized, and per-fission yields receive
no second branching factor. Fully exported fragment sets are supported.
Compound fission modes and nonzero SF RFS without a supported single-residual
interpretation refuse. The residual-lineage `DecayModel` is unchanged.

The assembly retains original evaluations, energy routing, immutable fragment
yield sources (including uncertainties), nuclide order and exact mode/channel
mappings. Parent-average emitted energy uses the shared evaluated routing;
branch Q is not additional heat. Ordinary residual resolution and optional
retained He4/H1 coproducts use the same helpers as the original decay assembly.
Native branch closure, allocation/topology limits and cancellation remain
active. No new evolution solver, inventory or thermal energy ledger is added.

`burnup_data` copies native declarations for explicitly appended irradiation
channels. `bind_burnup` admits reordered decay channels/products and added
neutron reactions while requiring unchanged identities, laws, groups, domain,
source labels, branches, multiplicities and energy. Product comparison is
bounded O(P log P), with cancellation during inspection. Original inventories,
finite irradiation, compartment movement and material adapters consume the
same native model and accepted trajectory.

`BoundEvaluatedBurnupDecay` projects accepted `BurnupStep` and
`BurnupTransportStep` reaction counts into source-ordered spontaneous-fission
events and exported fragment births. Reports retain each original event cell;
movement is never counted as fission. `SpontaneousFissionReportBudget` bounds
output values, work, per-cell container count and additional noncancelling
count-arithmetic defects. Native trajectory defects and source uncertainties
keep their original owners. Exported products have no represented future decay
or heat until explicitly modeled in their destination.

Fragment yields do not determine neutron multiplicities or spectra. Those
remain separately supplied inputs; ENDF spontaneous-fission neutron data in
MF1/MF5 are not synthesized from an MF8 average or a missing spectral record.
Seven manufactured regressions cover same-interval fragment decay and exact
parent-average energy, partial/complete exports, reversed source/channel/product
order with competing capture, ordinary charged-product parity, moving parents
with local exports, unsupported source declarations, bounds and late refusal.
All seven passed in the recovered original-source native suite alongside the
existing finite-inventory, movement and material regressions. Existing 64-node
native propagation limits remain unchanged.

## Evaluated spontaneous-fission neutron input

`endf::spontaneous_neutrons` imports identified NSUB=4 neutron multiplicities
and outgoing distributions using the original bounded ENDF record reader.
MF1/MT452,456,455 retain distinct total, prompt and delayed meanings; only
constant LNU=1, NC=1 multiplicities are admitted. Delayed LDG=0 decay constants
retain their original family order. An absent section is not a zero yield,
and total minus delayed is not silently substituted for an explicit prompt
source. The opaque types cannot enter the induced NSUB=10 importer.

MF5/MT18 prompt and MF5/MT455 delayed spectra admit energy-independent LF=1
tabulations with the ENDF-102 §5.4.2 INT=1,2,3 laws and zero-density endpoints.
Original component probabilities, source order, incident grids, interpolation
and outgoing points remain available after SI conversion. Every supplied
LF1 incident row must describe the same distribution. The LF5/LF7/LF9 extensions
below retain original parameters, an explicitly supplied reference and a finite
outgoing cutoff. Other LF laws and actual incident dependence refuse explicitly.
Neither probabilities nor spectra are
renormalized; their physical normalization belongs to a separately bounded
consumer. These quantities are per fission, so a mixed-mode parent's SF branch
must be applied exactly once when constructing its transformation source.

The source containers alone do not evolve delayed families, assign energy
from an MF8 average, or add aggregate delayed neutrons to explicit fragment
decays. Family lifetimes are not family abundances. Record, point, component,
family and cumulative interpolation limits bound allocation and source work;
cancellation remains active through the final publication checkpoint.
Format: [ENDF-102 (2023)](https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf),
§§0.3.2.4, 1.2–1.4, 5.1.1 and 5.4.1–5.4.2.

Five new input regressions, six existing induced-neutron input regressions and
the preceding 159 selected native regressions passed: 170 tests, zero failures,
with `thermal` enabled and no compiler diagnostics. The source selection uses
the original fs-la decay/prepared-LU and thermal timeline/advection owners.
It is not a full Cargo/constellation gate; the seven pinned sibling checkouts
and DSR/RCH executors remain unavailable in this workspace.

## Evaluated SF neutrons, radioactive fragments and material heat

`endf::decay::model::burnup::neutrons::EvaluatedBurnupNeutrons` prepares
original MF1/MF5 prompt sources and ordinary MF8 neutron emissions against one
`EvaluatedBurnupDecayAssembly` and an explicit `NeutronKineticModel`. Each
positive SF parent requires a matched prompt multiplicity and spectrum;
total or delayed multiplicities cannot substitute. Source Arcs and exact
MF1/MF5 metadata remain attached, and the original MF8 parent isotope/isomer,
MAT and atomic weight must agree. Compatibility checks do not certify the
physical consistency of separately supplied evaluations.

Every original component integral, component-probability sum and weighted
distribution integral must lie within the caller's explicit normalization
allowance. Admitted discrepancies are retained: diagnostics expose original
nu, integrated particles per fission, their signed difference, component
integrals and both original branch sums. Native decay-event counts include all
admitted parent branches, so the SF coefficient uses B_SF/B_all to produce nu
times the original native SF event count. This conversion neither changes
source probabilities nor applies another factor to independent product yields.
Ordinary MF8 spectra remain per parent transformation without a second branch
factor. Missing ordinary-neutron spectra and contradictory zero declarations
refuse; empty SF average-only STYP=5 markers do not establish spectral coverage.
Positive MF8 spectra on a mixed SF parent must identify disjoint supplied
ordinary decay modes. Aggregate delayed SF sources are never added over
explicit radioactive-fragment neutron emission.

MF5 and MF8 share the original analytical continuum moments and conservative
representative-energy projection. Preparation is bounded and cancellable;
mixed-mode attribution uses a prepared RTYP set rather than rescanning every
decay mode for each spectral record. The native projection retains particle
count and first energy moment, reports signed rounding residuals, and refuses
unsupported energy coverage. Only the original routed heavy-particle outward
energy can fund neutron kinetic energy. The same accepted reaction trajectory
transfers it once from radioactive radiation into transported neutrons.
Local heat, electromagnetic radiation and escape are not additional funding.

`bind_burnup` keeps the exact original nuclide order, laws, decay channels,
products and energy routes while admitting reordered channels and appended
neutron reactions. It returns the original `BurnupDecayNeutronEmission` for
the same finite inventory; no extra inventory, time solver, precursor history
or thermal ledger is introduced. Existing spatial transport, activation and
material heating consume that emission through their shared accepted step.
Fragment births can decay and emit neutrons during their birth interval;
capture acts on the evolving neutron population, with its deposited energy
reaching the original material calorics. Existing diffusion, group-condensation,
finite-inventory closure and native propagation limits still apply.

Seven new binding regressions and one material regression passed with the
preceding 170 tests and 15 original emission/interpolation regressions: 193
native tests, zero failures, with `thermal` enabled and no compiler diagnostics.
The manufactured mixed-mode chain checks same-interval prompt SF and daughter
neutron births against independent closed-form inventory and capture histories,
the actual native SF counts, disjoint energy transfer, complete material heat,
temperature and mass, split continuation and late cancellation/retry. Refusal
cases cover missing/overlapping sources, total-versus-prompt identity, original
heavy-energy limits, normalization discrepancies, nearly closed branching,
changed decay data and preparation bounds. This is original-source numerical
and integration validation, not empirical isotope calibration or a full
Cargo/constellation gate.

## General-evaporation spectra on the original decay trajectory

MF5 LF5 general evaporation retains the original dimensionless shape g(x),
positive temperature theta, restriction energy U, reference energy E and
interpolation metadata. Exactly one supplied temperature row establishes the
finite outgoing bound R=E-U; a constant temperature at different incident
energies does not establish a constant truncated spectrum. The parser admits
INT=1 through 5 with their actual logarithmic domain requirements. Singleton
metadata has no interpolation interval and may retain a zero reference energy.

The consumer streams E_out=theta*x and density g(x)/theta through the original
MF8/LF1 continuum moment integrator, without allocating a transformed spectrum.
Positive LF5 endpoints are valid. Positive support beyond R or the supplied
energy grid refuses, and source normalization is checked without clipping or
repair. The retained first energy moment enters the same conservative neutron
projection, finite radioactive inventory, capture and material heat transaction.

Seven additional regressions cover source retention and refusal, exact LF1/LF5
finite-history parity, probability conservation under shape/temperature
rescaling, original material heating and cancellation, and the six delayed LF5
families from the primary ENDF/B-VIII.0 Cf-252 evaluation. The unchanged MF1/MF5
fixture preserves its small normalization discrepancies; independently computed
50-digit decimal integrals check each original count and energy moment. This
checks evaluated spectral input, not aggregate delayed-family evolution or
empirical calibration of the manufactured inventory/material examples.

The focused original-source suite passed 200 native tests, zero failures, with
`thermal` enabled and no nuclear compiler diagnostics. All nine changed source
and fixture files matched the compiled snapshot. This extends the preceding
193-test seam; it is not a full Cargo/constellation or DSR gate.

## Finite Maxwellian and evaporation spectra in native decay coupling

MF5 LF7 Maxwellian and LF9 evaporation retain the original positive temperature,
single supplied reference energy and restriction energy. Their normalized
densities are respectively sqrt(E_out)*exp(-E_out/theta) and
E_out*exp(-E_out/theta), restricted to the declared 0<=E_out<=E-U. This finite
normalization defines each analytic law; component probabilities and evaluated
multiplicities retain the original explicit normalization allowance.

A shared bounded moment kernel integrates counts and first energies directly.
Small cutoff/temperature ratios use scaled coordinates, narrow groups use a
local series with an explicit complete-tail bound, and broad tail groups use
upper incomplete-gamma values. Original fs-math deterministic functions,
double-double exposure and scaled positive products preserve narrow positive
tails without subtracting nearly equal cumulative probabilities. Unsupported
energy coverage or lost positive prepared moments refuse explicitly. No outgoing
tabulation, quadrature mesh, new population solver or thermal ledger is created.

Prompt analytic sources pass through the same conservative group projection,
original SF event counts, finite parent/fragment inventories and material
heating transaction as LF1/LF5. Different neutron mean energies change capture
deposition and radiation while retaining the original radioactive energy total.
Delayed source parsing still does not establish aggregate-family evolution or
permission to duplicate an explicitly represented fragment's later emission.

Ten additional regressions were written: source metadata/refusal, independent
high-precision gamma references, tiny cutoff ratios, narrow/adjacent-float tail
bins, refinement and unit rescaling, bounds/cancellation, finite SF histories,
and full material heat/continuation/rollback. The recovered focused original-source
suite subsequently passed 30 native tests, zero failures: 10 source, 5 analytic
moment, 11 binding, 1 evaluated Cf-252 and 3 material regressions. All nine
feature-file hashes matched the committed source. The compile used sequential
code generation and disabled only unrelated test modules; its unused test-helper
warnings do not establish a full Cargo/constellation gate. No empirical or
performance claim is made.

## Delayed radioactive families, stored energy and material heat

The original joint neutron/precursor/nuclide propagator admits radioactive
parent-to-family births through BurnupDecayNeutronEmission::new_with_precursors.
Immediate and delayed yields retain exact nuclide and family order. Only native
spontaneous reaction fractions create these births; capture, material transfer
and exterior removal do not. Accepted family births enter native precursor
injection accounts, while native decay, transfer and outflow retain their
within-interval ages. Transported stocks require unchanged source identities,
spectra and yields in each compartment.

NeutronKineticModel::new_with_precursor_storage explicitly assigns each family
the future-neutron reserve sum(spectrum[g]*E[g]). Initial families and external
family feeds carry that reserve. Family decay converts it internally into the
original neutron spectrum, with no additional Q, heat or radiation. Transfer
moves the reserve; outflow exports it separately from neutron leakage. Existing
model constructors still reject delayed families, and positive induced-fission
delayed yields still require a separate physical funding declaration.

Finite activation reserves both immediate neutron energy and future family
energy from the original rate-weighted radioactive radiation allowance. It
reports internal births/reservations separately from external neutron/family
feeds. Neutron and precursor storage, internal conversion, original event
release and exported reserve close in the same accepted step. Cumulative
precursor outflow survives continuation and restart. The history entrypoint and
material transaction retain these accounts; only accepted deposited energy
heats the material, including capture after delayed release.

This raw native API requires caller-authored compatible energy and cohort
ownership. It does not infer aggregate ENDF-family funding or permission to
also retain the same radioactive fragments; the evaluated source adapter owns
those additional admissions.

Thirteen focused native regressions passed after runtime recovery, covering
analytical family birth/release, continuous feeds, capture competition and
physical refinement, precursor and parent transport/outflow, source and energy
refusal, histories, cancellation/retry, and original material temperature and
energy. The multi-interval capture regression composes the actual accepted
numerical defects for its complete-history balance; per-step and independent
reference-refinement checks remain unchanged. Production source matched the
committed precursor foundation. This selected original-source validation does
not establish a full Cargo/constellation or DSR gate.

## Evaluated delayed SF families on the finite material trajectory

`EvaluatedSpontaneousFissionDelayedNeutrons` prepares ordered MF1/MT455
lifetimes and multiplicity with the original MF5/MT455 component laws. Each
component creates nu*p*A precursors per fission, where A is its unmodified
integrated density. Dividing the projected component by A defines the
conditional native emission spectrum; it does not repair the original
multiplicity, component probabilities, or normalization discrepancies.
Original records, family order, integrals, projected moments and source energy
remain available. A tighter attachment allowance rechecks all original
normalizations. The existing analytical LF1/LF5/LF7/LF9 integrators and
conservative moment projection own preparation.

`with_delayed_spontaneous_fission` attaches those prepared families to unique
exact native family indices. Original nuclide/MAT/mass, group identity,
representative energies and all family metadata must match. The previously
admitted prompt source supplies the same B_SF/B_all coefficient, applied once
to convert per-fission births into native parent-decay births. Preparation
allocations and work are bounded; cancellation remains active through final
publication.

`DelayedSpontaneousFissionFunding::ParentHeavyRadiation` is an explicit caller
declaration of the maximum future-neutron energy per SF event. Both the
original spectral moment and the actual native family reserve must fit it.
Immediate emission plus delayed storage must also fit the original parent's
routed heavy-particle radiation. No local heat, electromagnetic radiation,
future fragment decay or additional Q supplies missing energy.

A positive aggregate delayed source is admitted only when every positive
independent SF product lies outside the retained nuclide cohort. The immutable
assembly already requires those omitted products to be explicitly exported.
Even a retained stable fragment refuses this aggregate boundary: ownership is
declared for the whole fission-product cohort, not inferred by selectively
disabling its neutron channels. Explicit fragment decay remains the alternative
for a represented cohort; the two histories cannot be superposed.

Binding returns the original native radioactive emission owner with both
immediate and family birth coefficients. The same accepted joint propagator
depletes the parent, ages each family, emits neutrons, captures them and reserves
or releases energy. Material heating receives only deposited energy at the
actual accepted capture time. Restart retains the original inventory,
precursor reserve and prior released energy without replaying its heat.

The focused native consumer suite passed 20 tests, zero failures, with
`thermal` enabled: five new adapter regressions, the new full material history,
11 existing prompt/LF5/LF7/LF9 binding regressions and three existing material
regressions. Exact changed-source hashes matched the compiled snapshot. The
new cases check original component discrepancies, family ordering and mapping,
exclusive product ownership, explicit energy funding, preparation bounds,
late cancellation, finite source residence times and split continuation. The
material case checks independent closed-form populations, both neutron groups,
reserved energy, deposited heat, radiation, temperature, mass and restart.

The actual ENDF/B-VIII.0 Cf-252 six-family LF5 fixture also passed preparation
through the new adapter. Independent high-precision decimal moments check its
original yields and energy; its small published normalization discrepancies
remain visible and a tighter allowance refuses them. That preparation-only
case does not infer an MF8 energy budget or a fission-product cohort.

The corrected 13-test native foundation and this 20-test consumer run use
selected original-source modules and the existing native numerical owners.
Unused test-helper/import warnings arise from the narrowed selection. These
runs are not full Cargo/constellation, DSR, empirical calibration, or performance
benchmarks.

## Explicit future-neutron conversion in induced fission

`NeutronKineticModel::new_with_fission_precursor_funding` admits original
positive delayed-fission yields under reaction-ordered declarations of
`InducedFissionPrecursorFunding::AdditionalDelayedConversion`. The declaration
bounds the actual sum of each family yield times its reserved neutron energy.
An unused part of the cap is never generated or deposited. Positive products
that underflow, insufficient/nonfinite caps and funding on nonfission reactions
refuse, including channels whose current cross section is zero.

The existing `NeutronReactionEnergy::nuclear_q_j` remains prompt Q. Prompt
threshold checks and E_in + Q_prompt - E_prompt_out deposition are unchanged.
Additional future-neutron conversion enters precursor storage at the accepted
induced-fission event, independently of prompt Q and external source injection.
Family decay subsequently transfers that reserve to the original neutron
spectrum; transport and outflow retain their existing separate energy accounts.
This partition follows the separate prompt and delayed neutron components in
[ENDF MF1/MT458](https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf) and the
[OpenMC energy-deposition description](https://docs.openmc.org/en/stable/methods/energy_deposition.html).
No delayed beta/gamma energy, evaluated-data match or calibration is inferred.

Native family births are reconstructed from the accepted original fission
reaction counts and yields, including their actual rounded products. Each
family balance is audited independently. Whole-step energy includes the
additional conversion once, and state/history retain separate cumulative and
per-cell accounts through continuation and restart. Prepared positive-source
reaction indices avoid rescanning unrelated capture/scatter reactions for this
account. Evolution remains with the original native diffusion propagator.

Legacy kinetic constructors retain their original refusals. Finite and
radioactive-inventory consumers must additionally establish exclusive ownership
of the fission-product cohort; a raw energy declaration alone does not do so.

The focused native core suite passed 15 tests, zero failures: five new induced
family tests, four existing prompt-fission tests and six existing precursor
tests. Independent analytical chains cover prompt versus delayed neutron
residence, capture and heat, unused funding, original thresholds, outflow,
history/restart and late cancellation. All six changed Rust files matched the
compiled snapshot; narrow formatting passed. This original-source selection
is not a full Cargo/constellation or DSR gate, an empirical isotope calibration,
or a performance benchmark.
