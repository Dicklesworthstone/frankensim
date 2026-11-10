# fs-thermal contract

## The enthalpy timeline as a composition-kernel Component

`component::ThermalTimeline` puts the explicit `DepositionSimulation` behind
`fs_ports::component::Component` (contract v0), so the kernel in `fs-system`
can deliver heat to its cells from any other component and route what leaves
it. Nothing is integrated in the wrapper: a step is the owner's own
`advance_observed`, with the delivered heat given to the owner as its
`step_sources`.

Ports, in this order: `cell<id>` per declared cell (a heat port that receives
the heat rate into that cell and returns the mean of its temperatures at the
two ends of the step), and `environment` if declared (returns the mean rate at
which heat entered through the timeline's own authored boundaries: ambient
films, reservoirs, authored sources; the temperature it receives is not used
by the timeline).

Boundary heat is on that wire or it does not exist. Without the port, a
timeline in which any active cell has a film coefficient, a reservoir or an
authored source is refused at `ThermalTimeline::start` and at every step
(`TimelineRefusal::Boundary`), and so is a step in which an event gives a cell
such a boundary. The test is on the declared boundary, not on the heat it
exchanges, so exchanges that cancel or are still zero are refused too. While
no boundary exchanges, the environment port returns exactly zero.

The state is `TimelineState`, the owner's checkpoint, made only by `start`
(which applies the events scheduled at the current time: the births a run
starts with) and by a step.

Interval and energy: the owner's clock is one `f64`. A step of `dt` from time
`t` ends at the representable `t + dt`, and the owner integrates over its own
substep durations, whose sum is not exactly `dt` once `t` is large beside
`dt`. The kernel's energy is kept, not its rate: each port's `rate * dt` is
handed to the owner substep by substep, in proportion to the clock interval
each substep covers, and what is left goes into the substep that ends the
step. A step that does not move the clock is refused
(`TimelineRefusal::Interval`). A step of exactly the distance `next_event_s`
reported ends on that event.

`limit_dt_s` returns the smaller of the proposal and the owner's explicit
bound for a substep from that state (which includes its `max_dt_s`), read from
the owner without taking a step. A kernel that honours it normally takes one
thermal substep per step, so the rest of the system advances on the thermal
clock, as the hand-written couplings of this timeline do; a longer step is the
owner's several substeps. `next_event_s` is the time to the next scheduled
event.

The receipt's stored change is the owner's `EnergyAudit`; numerical loss is
what entered through the ports less that. The energy tolerance is
`(64 + 2 cells) eps` times the sum of the sizes of the port energies, of the
terms `|m e|` of the stored energy at both ends (they cancel when laws have
datums of opposite sign) and of the owner's totals.

G0, in `crates/fs-system/tests/thermal_timeline_slice.rs` (run with
`cargo test -p fs-system --features slices`); the list is in
`crates/fs-system/CONTRACT.md`.

No claim: the entropy generated inside the timeline is reported as zero. That
is a lower bound for conduction between cells; for the boundaries it is one
only if the environment wire is at the authored ambient temperature, which is
not checked. A birth after the start is refused, because its mass and enthalpy
would enter with no port to carry them. The class is `Smooth`; a law with a
latent interval has kinks at the interval's ends, where a kernel that
differences the port map gets a one-sided slope. No analytic tangent. The
implicit and staged variants of the owner are not wrapped. Beside port heat,
the environment heat carries that heat's rounding (about 1e-16 of it), since
it is the owner's external heat less the port heat. Against the owner driven
directly with held rates, a step agrees to rounding, not bit for bit: the
wrapper hands the owner energy. Not pinned by a test: the sign of the
numerical loss (it is at rounding level); every term of the tolerance but the
stored-energy terms, and its per-cell factor (the plain stored-energy sum
needs it from about 1e6 cells); equality's use of the budget, and its
blindness to the schedule past the next event and to contact coverage; that
the step limit reads the bound without a trial step (a trial step would
return the same bound). fs-thermal depends on fs-ports and fs-exec for this
module only.

## Shared immutable inputs across thermal candidates

Every clone of the original `DepositionSimulation` shares its immutable
resolved event schedule, borrowed caloric-law table and admitted global
temperature-slope bounds. Copying these authored inputs takes constant work
independent of future event count or registered material count. Implicit
residual evaluations, line-search trials, adaptive stages and phase-transfer
candidates all retain that sharing through the existing timeline owner.

Each candidate still owns its active cells, enthalpies, geometry, links,
surface coverage, physical time and event cursor independently. Accepting one
branch or refusing a cancelled trial cannot mutate another branch or consume
its future births. The original borrowed law identities and lifetimes are
unchanged; laws remain subject to the existing immutability contract. This
single-threaded storage introduces no new thread-safety guarantee.

The focused regression follows explicit split/replay trajectories through a
birth, contact and boundary change, cancels after candidate work, retries from
the original state and accepts an implicit enthalpy endpoint. It checks input
allocation identity together with independent clocks, events and heat storage.

## Stationary operating points of original conduction networks

`ImplicitConductionChain::solve_operating_point` and the graph counterpart
solve each original cell's net heat rate in watts, using the same analytic
energy tangent and direct tridiagonal or sparse Krylov owner as native implicit
evolution. Original calorics, nonlinear conductivity, contact resistance,
radiative/composed boundaries, density geometry and local/nonlocal source
feedback remain bound. Physical clocks and scheduled inputs are frozen at the
queried instant. Bounded Newton/backtracking checks original material domains
and both cellwise and total heat residuals; no fictitious timestep, clipped
temperature, apparent latent capacity or regularized null space supplies a root.

Every currently conducting component must reach an active temperature-dependent
bath. Isolated inventories require conserved-energy constraints and are refused
by this lane. A stationary latent state retains its supplied enthalpy; a singular
off-equilibrium plateau is not assigned a phase fraction. Driven steady states
can have nonzero heat throughput and entropy production. Residual acceptance
does not imply uniqueness, dynamic stability or continuum accuracy.

The immutable `ThermalOperatingPoint` exposes ordered original enthalpies,
temperatures, actual cell residuals and work counts. It cannot be consumed as an
evolved timeline or heat credit. Its `prepare_thermal_impedance` method feeds
the exact retained material/source state into the existing harmonic owner,
whose resource, interior-domain and residual gates still apply. The input
timeline and future event cursor remain unchanged on success and refusal.
Advection is outside this stationary lane. G1/G3/G4 tests cover analytical
heated chains, nonlinear conductivity, sparse nonlocal feedback, impedance
handoff, latent/null-space behavior, source-domain limits, cancellation and
work ceilings.

## Closed-component thermal equilibrium with conserved original storage

`ImplicitConductionChain::solve_isolated_equilibrium` and the graph counterpart
solve common temperatures with a separate conserved inventory for each original
connected component. They preserve every mass, material binding, geometry,
physical clock and future event. The conserved row is
`sum m_i*(e_i-e_i,initial)=0` in the original stored coordinate: enthalpy,
explicit fixed-volume internal energy, or declared sensible storage. No
pressure-work conversion or assertion that every coordinate is internal energy
is introduced. Original-coordinate differences and their mass products use the
existing `fs-math::Dd` owner, with independent physical residual acceptance.

The reusable `equilibrium::isolated::caloric` core queries each original
continuous nondecreasing `T(e)` on its declared domain. Two bounded generalized
inverse searches retain both sides of an exact latent plateau, including when
a material selects a sensible one-sided derivative at a corner. A scalar
temperature solve uses original slopes and safeguarded bisection; latent heat
never becomes an apparent heat capacity. The source-evaluation ceiling is
shared across disconnected components. Original source and checkpoint refusals
propagate.

An already equilibrated component retains every original energy exactly.
Otherwise one latent parcel is determined by the remaining conserved energy.
Several independent latent allocations with a strictly interior inventory
refuse; a loose energy tolerance cannot select their phase distribution.
Conservation can uniquely determine an all-solid or all-liquid endpoint. An
unrepresentable original-coordinate balance refuses rather than being rounded
into a successful root.

The native owner reuses its original endpoint and heat-rate evaluators. Every
authored link must conduct initially and at the returned state, and exact
component membership is rechecked before each cell's original watt residual is
accepted. Disconnected components and isolated cells with no connecting link
are supported. Authored insulating links, changed conducting topology, external
baths/sources, advection and custom flux links remain outside this lane.
Temperature-dependent native Fourier/contact coefficients retain their
original material meaning.

The result contains a read-only `ThermalOperatingPoint` and per-component
conservation reports. It supplies neither a time-integrated heat/work credit
nor a phase-evolution or reachability claim. Scalar regressions cover unequal
masses and capacities, nonlinear calorics, large independent energy datums,
latent allocation and endpoints, domain/work/cancellation refusal and
unrepresentable exchange. Native integration regressions cover disconnected
inventories, retained event history, latent refusal, the independent watt gate
and closing material conductance.

The exact production caloric core compiled against the real `fs-math` and
passed all eleven focused native regressions. The five complete graph/chain
integration regressions are authored and syntax-checked; their native execution
remains unverified because the original runtime dependency build was killed.


### Closed-component continuous response at fixed energy

An accepted isolated equilibrium now exposes its own
`prepare_thermal_impedance` path and
`harmonic::isolated::IsolatedThermalImpedanceModel`. The ordinary preparation and
its unanchored-DC refusal are unchanged. Both paths retain the original material
caloric slopes, graph tangent, complex LU, and complete original cell-power
residual checks.

For peak phasors with exp(-i omega t), each closed conducting component obeys
-i*omega*sum(delta E)=sum(delta P). The isolated factor replaces one dependent
power row by the component's extensive-energy sum. At nonzero frequency its
right side is i*sum(delta P)/omega, retaining the physical storage pole under
net heating and a well-conditioned balanced low-frequency limit. At DC every
component's net heat must vanish within the explicit watt tolerance. The
default response fixes its energy increment to zero;
`response_with_component_energy` admits an explicit increment per component
only at DC. `heat_transfer_impedance` supplies a compatible unit transfer within
one component, without pinning any cell temperature.

Component energy is the original stored coordinate (including enthalpy for a
material that declares enthalpy), not a silently substituted internal energy or
fitted capacity. Original component membership is retained and the represented
tangent is audited for component isolation and conservative column sums.
Compensated original heat/energy sums retain small residuals beside large
inventory responses. Independent absolute joule and component-watt gates supplement every
unmodified native power row, including each replaced row. A separate absolute
cell-watt gate prevents a large storage response from hiding a finite power
defect behind a small relative backward residual. Their tolerances never
scale silently with the 1/omega pole; an unrepresentable requested response
refuses.

A single latent parcel's DC energy can be determined by its component inventory.
Several latent parcels in one component leave unspecified phase allocations and
refuse DC factorization. Their finite-frequency energy response remains valid
with zero temperature response on the admitted sharp plateaus. One-sided
caloric corners do not establish two-sided harmonic derivatives.

This is a local linear response about the exact accepted isolated operating
point. It adds no bath, leakage, artificial heat capacity, pseudo-time step,
finite-rate chemistry, trajectory, heat/work credit, phase-selection policy,
stability or reachability claim. The caller keeps amplitudes inside the admitted
differentiable material branch. The dense numerical payload allowance is
64*n*n+512*n bytes, with the same 1024-cell cap and native work/cancellation
limits.

Ten production-core analytic regressions and five native graph/material integration
regressions are included for conservation, DC heat transfer, finite-frequency
storage, unequal capacities, latent allocation and failure atomicity. These
new tests remain unexecuted after the execution workspace disconnected.

## Continuous prescribed bath-temperature histories

`BathTemperatureSchedule` and `with_bath_temperature_histories` on native
implicit chain/graph owners prescribe bounded continuous piecewise-linear
positive-Kelvin histories on existing reservoir and ambient channels. Original
calorics, geometric conductivity, material-density half-distance derivatives,
integrated/radiative boundary laws and native exposed-area Robin coefficients
remain active. Each BE endpoint and SDIRK stage samples its physical time;
fixed, observed and adaptive intervals clip at history knots. Scheduled
single-bath power replaces the original channel through the existing boundary
exchange owner, with one heat credit and entropy at the actual bath temperature.
Passive SDIRK entropy admission retains those actual stage exchanges. A
composed law's separately declared bath remains law-owned.

`GrayTwoBathBoundary::with_radiative_temperature_history` now binds a second,
independently prescribed radiation-temperature schedule through that original
exchange owner. The native chain/graph queries its physical endpoint or SDIRK
stage time, clips fixed/observed/adaptive intervals to its knots, and retains
both actual bath temperatures in the entropy quadrature. Convection and
radiation still resolve one shared exterior temperature and half-cell/coating
drop; conductivity, source-emissivity and density-driven geometry derivatives
remain in the same analytic enthalpy tangent. The borrowed schedule introduces
no new surface storage or duplicate heat credit. Timeless queries on a
scheduled law, exhausted time support and cancellation refuse, including zero
emissivity cases. Autonomous exchange laws retain their original arithmetic
through default time-aware trait queries. Focused G1/G3 regressions compare
both stage balances with independent scalar fourth-power solutions and check
actual-bath entropy, knot replay, adaptive refinement and rollback.

Accepted states
retain authored boundaries; conflicting native temperature events, outside-domain
queries, invalid samples, budget exhaustion and cancellation refuse
transactionally. Histories never imply extrapolation or additional heat sources.

Four exact-production native Rust checks passed for schedule interpolation,
support/admission, floating-point extremes and one-credit actual-bath entropy.
Four pipeline regression groups cover physical stage sampling and refinement,
nonlinear/radiative/composed tangents, latent-plateau Robin cooling, and rollback.
Those pipeline regressions were authored but not executed against the complete
Cargo graph in this environment.

## Sparse coupled heat-generation feedback in implicit conduction

`ImplicitCoupledHeatSourceLaw` and the chain/graph
`with_coupled_heat_source` attachment evaluate a declared set of original
cell calorics at one simultaneous enthalpy state and physical stage time.
The source returns signed cell powers plus complete derivatives in its
immutable sparse dependency order. A declared derivative budget bounds
additional sparse work; original caloric identities, cell coordinates and
reference masses remain bound. Non-neighbor feedback requires the existing
graph solver and is refused by the tridiagonal chain before a source query.

The original enthalpy residual adds each coupled power once alongside
independently owned authored/observed/local heat. Its Newton entry is
-dt*(dP_row/dh_column)/m_column because the native unknown is an energy
increment. Both graph forward and transpose actions consume the complete
entries. BE, both physical SDIRK stages and the adaptive time owners retain
their existing storage, quadrature and rollback owners; declared continuous
forcing transitions clip the same native time grid.

This enables original domain owners to supply coupled feedback such as
fixed-voltage Joule heating without freezing neighboring temperatures or
dropping cross derivatives. It does not introduce a thermal-side electrical
solver, infer reaction inventory, or credit upstream work twice. Calibrated
source validity and upstream energy ownership remain explicit. Source
entropy stays unclassified even when separately active contributions cancel
in their net cell power, including retained first-stage SDIRK heat.

Three exact production-leaf native Rust tests passed for ordered sparse
binding, cross-cell Joule residual derivatives with reversed coordinates and
unequal masses, and shape/domain/budget/cancellation admission. Full native
integration regressions cover independently eliminated circuit/conduction
endpoints, original graph forward/transpose actions, scheduled fixed/adaptive
BE/SDIRK and rollback, and source-cancellation entropy retention across both
stages. These integration definitions are unexecuted in this workspace's
incomplete native build.

## Crystallization in spatial heat exchange and finite molds

`enthalpy::crystallization::CrystallizationSimulation` can attach before the
first scheduled birth or adopt active all-amorphous cells from the original
`DepositionSimulation`. Bound materials receive a private
kinetic exposure; unbound molds/substrates keep their original caloric laws,
storage basis and thermal slope bounds. All cells exchange heat through the
existing link, boundary and CFL operator. Mass, geometry, cell/material IDs and
the original extensive thermal-storage account remain with that owner.
Attachment requires the original thermal law and kinetic parent branch to
carry the same complete caloric identity. Agreement at one starting state is
insufficient; unidentified or crossing-but-different caloric curves refuse.

Birth, contact, separation and boundary events use that original timeline's
admission and coverage owner, in authored order. Every thermal/kinetic interval
ends at the next event or target; target-time births activate without aging.
Each newly deposited bound cell defaults to explicitly zero exposure and its
unchanged incoming enthalpy. Existing cells retain their exposure, and ordinary
molds receive no crystal history. `new_with_incoming_states` instead accepts
original `CrystallizationState` histories keyed by active cell or unique pending
birth. The complete law identity and caloric reconstruction must validate and
stored enthalpy must match exactly; no age or fraction is inferred from
temperature. Missing keys, repeated births, ordinary-material targets, stale
enthalpy and another law refuse the complete attachment. The immutable feed map
is shared across cloned and implicit candidates, avoiding per-trial copies of
the incoming phase inventory.

Supplied feed retains its original exposure and crystal mass without aging
before birth, including a birth exactly at the advance target. Its validated
mixture temperature goes through the original birth admission owner, retaining
all geometry, material, boundary, coverage and work-budget rules. This permits
cold crystal mixtures whose enthalpy is outside the all-amorphous inverse chart;
the parent chart never reinterprets their phase energy. Imported crystal mass
is already part of incoming total mass, so `crystallized_mass_kg` counts only
subsequent signed conversion. No extra latent heat is applied at ingress.
Focused regressions cover active-state adoption, cold seeded feed and ordinary
molds, exact endpoint birth/restart, implicit continuation, source and enthalpy
mismatches, and rollback on late original admission failure or cancellation.
Empty intervals consume no thermal
steps; events do not renew the whole-advance substep budget. Incoming mass and
birth enthalpy enter their original accounts once, while heat residuals retain
cell-relative differences. A late birth-domain refusal or cancellation rolls
back the entire candidate, including its schedule cursor. Manual endpoint
boundary changes require due events to have been processed first.

Each candidate resolves crystallinity and temperature together through the
shared constitutive core. Crystal mass is part of the existing specimen mass,
and latent release is part of the existing enthalpy. Compensated audits count
both signed global residuals and cumulative absolute local caloric/roundoff
defects. Late material failures, cancellation and exhausted budgets return no
accepted history. Cloning checkpoints time, enthalpy and exposure; identical
time partitions replay bitwise. Endpoint boundary/conductance changes use the
original physical ownership and coverage checks.

`with_equilibrium_remelting` binds a prepared `EquivalentAgeRemelting` to
one original material. The full source-derived capability must match that
material's kinetic law; selection survives cloning and scheduled births.
Other materials retain their irreversible growth behavior. Existing heat
links and boundary channels can now remelt earlier crystals when later hot
material is deposited. The signed crystal-mass observation becomes negative
during net melting; it is neither an added mass source nor another heat debit.
Three timeline regressions cover latent heating followed by superheat and
recrystallization, hot incoming material remelting an older deposit through
finite contact, and original-history preservation on binding/cancellation
failure. These timeline regressions require the complete native dependency
graph to execute.

Focused regressions cover adiabatic and spatial calorimetry,
recalescence while cooling, different-Cp finite-mold exchange, timestep
refinement, source/identity/budget refusal, initial-chart roundoff, ordinary-cell
domain failures and atomic cancellation/restart, as well as independently aged
scheduled births, finite-mold contact/separation and continuation across events.
This fixed-reference-geometry,
uniform-cell reduction has first-order heat transport and explicit kinetic
timestep limits. Conductivity evolution requires a declared supplied law;
front morphology, moving geometry, mechanical unlocking and work feedback retain
their separate physical owners.

### Implicit kinetic heat transport on the existing network

`enthalpy::crystallization::implicit::CrystallizationBackwardEulerStep` extends
the original backward-Euler energy residual with complete original-state
calorics. Each joule-vector candidate starts from the same accepted cell
histories and advances the shared adaptive phase kernel once over that interval.
The existing thermal operator still owns links, covered boundaries, sources,
storage coordinates and residual assembly. Ordinary finite molds and substrates
retain their own enthalpy or explicitly fixed-volume internal-energy laws.
No fictitious crystallinity is assigned to those materials, and no parent-law
entropy integral is applied to a kinetic mixture.

`solve_implicit_step` reuses the existing fallible finite-difference and
Newton--Krylov owners for simultaneous network heat exchange. Endpoint source
callbacks see every actual proposed phase state. A bounded two-entry residual
cache avoids repeated phase reconstruction at identical numerical base points;
acceptance reconstructs the complete physical state and source payload afresh.
Both endpoint caloric errors, each original local heat residual and independently
reduced storage/transfer roundoff consume noncancelling energy budgets. Numerical
convergence alone cannot publish a physically inadmissible candidate.

`advance_implicit_to` splits at original birth/contact/separation/boundary events
and processes them through the existing timeline in authored order. New material
enters with its original mass/enthalpy and zero exposure. Whole-call nonlinear
work, physical interval count and accumulated energy defects remain bounded;
late source, solver or checkpoint failures return no partial timeline.

`advance_implicit_with_coefficients_to` uses that same event-aware march with
original state-dependent powers, link conductances and boundary-channel
conductances queried inside every nonlinear endpoint. Source callbacks must
keep their physical history immutable during speculative calls. The return
retains the payload used by the last accepted physical interval, or `None`
when no interval was solved, without a later source requery or an unbounded
per-step log. Exact endpoint events follow that interval, so its payload
describes the pre-event endpoint. A later source refusal rolls back the whole
request, including earlier accepted intervals and events.

The implicit path requires explicitly selected phase adaptivity for every
active kinetic material. It removes the explicit heat-transfer CFL restriction;
maximum timestep, phase tolerances and the joule-valued numerical-difference
scale still require independent resolution. It retains fixed reference geometry
and mass between births, with no moving-interface or pressure-work model inferred.
Four integration regressions specify a three-body kinetic/finite-bath analytic
balance, deposition/contact event replay, original refusal/cancellation, and
nonlinear contact/source continuation with original heat accounts and rollback.
The complete thermal library and this public integration-test target pass native
metadata compilation against the original dependency artifacts. That verifies
type integration; it does not execute the new network regressions.

### Instantaneous heat transfers on accepted phase history

`CrystallizationLaw::apply_specific_enthalpy` resolves a caloric transition at
exactly zero elapsed time. Irreversible exposure is frozen during heating and
cooling; explicitly selected equilibrium remelting can reduce surviving
crystals and their equivalent age. An unchanged enthalpy preserves the accepted
state bitwise after source, state and control validation. The original inverse
and melting projection own sensible and latent energy together.

`CrystallizationSimulation::apply_heat_impulses` applies one signed joule amount
per original active cell in deterministic physical-ID order. Bound cells use
that same phase closure; ordinary cells retain their original caloric law and
storage basis. The atomic candidate preserves time, event cursor, mass,
geometry, contacts and boundary sources. Due events must already be processed;
conduction and kinetic growth resume only with actual elapsed-time advancement.
No arbitrary impulse duration is introduced.

Cell-relative storage differences and compensated sums separately bound lost
heat increments and original caloric defects. Opposing local errors cannot
cancel admission failures. Caloric iteration and original cell budgets bound
work; a late source refusal, energy failure or caller checkpoint error returns
no partial batch. The audit's `external_heat_j` is relative to the thermal
subsystem: a coupled mechanical owner must cancel that internal transfer with
its own storage debit and account for any nonthermal export. Three timeline
regressions cover finite-mold calorimetry and later contact, latent remelting
and recooling without kinetic aging, and atomic refusal/cancellation including
unrepresentable heat increments. They require the full native thermal graph.

### Kinetic phase evolution against a constant pressure load

`crystallization::pressure::IsobaricCrystallization` binds the original phase
law to a declared constant absolute pressure. The optional original-density
capability on each caloric branch must be available and identified; no density
or expansion coefficient is inferred. The caller admits the caloric, density
and kinetic data at that same pressure and composition. This is not a
variable-pressure EOS or a rigid sealed-volume caloric law.

Signed heat changes the original enthalpy once. Instantaneous transitions reuse
the original frozen-age/remelting path; finite-time transitions reuse its
adaptive history integrator. Original phase densities and the evolving crystal
mass fraction determine volume. A symmetric density-difference calculation
retains sensible and phase-conversion contributions without subtracting rounded
mixture volumes. Internal-energy change is defined by delta_H-p*delta_V on
that same source convention; it is not independently measured EOS evidence.

The published-volume boundary work is reported separately. Its difference from
the source-volume increment, heat quantization and caloric defects consume the
energy budget without cancellation; a separate volume budget also applies at
zero pressure. In adiabatic growth, constant H can therefore coexist with
contraction work and an increase of U. No pressure work is subtracted from the
enthalpy target a second time.

`pressure::parcel::IsobaricCrystallizationParcel` owns the continuing accepted
mass, pressure, phase history, actual volume and elapsed time. Its immutable
mechanical reference volume remains distinct from the evolving free volume.
Clones retain retry state; a late source/work/accuracy refusal or cancellation
publishes no partial parcel. Instantaneous heat does not age it, while positive
elapsed intervals can crystallize even with zero heat. This uniform quasistatic
parcel supplies neither a moving mesh nor piston inertia, pressure evolution,
spatial continuity or fluid momentum.

### Implicit heat contact between finite kinetic parcels

`pressure::exchange::exchange_heat` solves one shared signed heat transfer
between two original isobaric parcels. Each nonlinear candidate advances both
accepted inputs with opposite heats through their original adaptive phase
owners. The backward-Euler contact equation uses both proposed endpoint
temperatures; recalescence can reverse the initial thermal ordering. A finite
donor cools as a receiver heats or melts. Zero conductance still advances both
adiabatic histories, without creating a heat source.

The caller admits a finite total conductance and a heat bracket that covers a
root within both original domains. A bounded safeguarded scalar solve does not
infer a conductivity mixture rule, contact area, moving geometry or uniqueness
of the nonlinear root. Source failures propagate, comparison states never
replace accepted inputs, and final cancellation publishes neither body. The
contact residual is distinct from the caloric and combined noncancelling energy
budgets. Opposite enthalpy transfers cancel while each original pressure load
retains its own volume work; latent heat and pressure work are not additional
enthalpy sources.

Six native leaf regressions cover stiff unlike-capacity equilibration,
recalescence-driven reversal, latent remelting and pressure work, insulation,
represented time, source/domain/work refusals and atomic cancellation. Two
actual-material integration cases use original caloric, density and entropy
owners for contact-driven melting and crystallization; they require the full
native material graph to execute. The single endpoint step has first-order
thermal time error, separate from its inner phase-history estimates.

`pressure::exchange::adaptive::advance_exchange_to` adds local thermal time
control around that same endpoint owner. One full contact step is compared
with two half steps, and only the conservative half steps are retained.
Temperature, specific enthalpy, crystal mass fraction and unsaturated exposure
all participate; a melting plateau cannot pass on temperature agreement alone.
Relative tolerances scale physical changes instead of absolute caloric datums.
Representable unequal half spans receive the first-order error correction.

Rejected comparisons spend bounded work but contribute no heat, phase history
or volume work. Double-double accumulation retains accepted physical accounts;
final endpoint storage/work and public-total rounding consume a separate
whole-call noncancelling energy budget. Exact target times, deterministic
replay, maximum comparison work and final cancellation are explicit. The
caller must admit the heat bracket over every candidate. Original source,
domain and nonlinear failures propagate instead of being disguised as time
refinement. These local discrepancies are not global error bounds, and the
maximum timestep must resolve the original source's features. Three native
regressions cover convergence to finite-body exponential cooling, latent
inventory-driven refinement and atomic work-budget/cancellation refusal.

`exchange_heat_with_conductance` and
`adaptive::advance_exchange_to_with_conductance` evaluate an admitted original
conductance source on BOTH proposed endpoint parcels. Temperature, kinetic
mass fraction, original material identity and actual volume remain available
to that source; crystal fraction is never substituted for equilibrium liquid
fraction. The source must be deterministic and physically immutable during
speculation. Its finite nonnegative result enters the same contact residual,
and its exact refusal propagates. The fixed-conductance and source-query APIs
share one scalar owner and one adaptive controller, with the existing fixed-G
insulation shortcut retained.

These APIs do not invent conductivity mixing or interface geometry. Source
continuity and the admitted heat bracket remain caller obligations; successful
endpoint samples cannot prove global source validity or root uniqueness. Two
native endpoint regressions exercise a nonlinear radiation-shaped conductance
against an independent cubic oracle and source/cancellation failures. An
adaptive regression compares nonlinear finite-body cooling with its analytic
ODE solution and verifies constant-source/fixed-G replay equivalence.

## Finite-rate crystallization and latent thermal feedback

`crystallization::CrystallizationLaw` retains an independent Nakamura exposure
`z = integral K(T) dt` and absolute crystal mass fraction
`X = Xmax * (1 - exp(-z^n))`. Its two original caloric branches share one datum:
`h = (1-X) h_parent(T) + X h_crystal(T)`. Distinct heat capacities and the
original temperature-dependent enthalpy gap enter this one account; latent heat
is never added again as an independent source. Relative conversion and maximum
absolute mass crystallinity are separate parameters.

Isothermal stepping returns the required thermostat heat. Nonisothermal
stepping integrates exposure with explicit midpoint while reconstructing both
midpoint and endpoint temperatures on the same enthalpy account. Work and
caloric tolerances are bounded, cancellation publishes no candidate, and
out-of-domain or unsupported remelting trajectories refuse. Immutable source owners may
supply positive whole-domain heat-capacity lower bounds to accelerate caloric
inversion; generic laws retain safeguarded bracket convergence. Sixteen focused
native core Rust tests pass, covering analytic JMAK, adiabatic recalescence with
unequal heat capacities, second-order time refinement, accelerated affine and
nonlinear inversion, datum shifts, endpoint behavior and cancellation.

The declared source domain may extend above melting for completely amorphous
feed. Both original caloric branches and the kinetic source must cover that
full domain; no chart is extrapolated. The supplied Nakamura rate must be
exactly zero at and above melting, and positive queried rates refuse instead
of being clamped. Once any nonzero crystalline mass fraction exists, direct
state queries and enthalpy inversion restrict temperature to the original
domain at or below melting. Cooling superheated amorphous feed enters the same
conservative crystallization account. The default growth path refuses
remelting; the separate capability below explicitly selects that model. The core regressions include superheated
cooling, source-domain/rate refusal and crystal-preserving melting-boundary
checks. A growing-timeline regression deposits 550 K feed against a finite mold
and follows sensible cooling into crystallization without a second latent heat
source; it requires the full native thermal dependency graph to execute.

The default is an irreversible isokinetic constitutive model with explicit
caller parameters and applicability. It supplies no material calibration, stochastic
nucleation, grain/dendrite morphology, segregation, remelting/reset or mechanical
locking. Time integration error remains distinct from caloric closure.

### Explicit equilibrium remelting and thermal cycling

`CrystallizationLaw::equivalent_age_remelting` admits and caches the original
melting-temperature calorics in an immutable `EquivalentAgeRemelting`.
Both branches and the original zero-at-melting kinetic source must cover that
temperature. `advance_specific_enthalpy_with_remelting` preserves the original
enthalpy at its midpoint and endpoint while latent heat is absorbed on that
plateau. Surviving crystal MASS fraction is obtained from the original forward
Avrami relation after a bounded equivalent-age reconstruction; complete melting
returns zero effective age and the original amorphous branch. Cooling resumes
the original growth law from that age.

This separately chosen history approximation supplies no retained-nucleus,
lamellar-population, dissolution-rate, molecular melt-memory or mechanical
unlocking model. Smooth unmelted paths retain the original midpoint method;
transition steps require refinement and carry no second-order event-accuracy
claim. A projected midpoint that would suppress continuing undercooled growth
refuses. Zero heat and zero growth rate preserve an accepted plateau state
bitwise, including an admitted caloric residual, without another inverse solve.
Five native core regressions cover partial/complete melt, superheat and renewed
growth, conserved energy, source/domain/history refusal, large datums and
saturated ages, false-arrest prevention, cancellation and stationary holds.

### Adaptive phase history on a prescribed enthalpy path

`advance_specific_enthalpy_adaptive` reuses the original irreversible or
explicitly selected remelting kernel. One full step and two half steps follow
the same linear enthalpy path; acceptance retains the two half steps without
extrapolating temperature, crystal mass fraction, exposure or energy. Absolute
temperature, fraction and exposure discrepancies each have an explicit local
tolerance. Exposure remains part of the comparison when a small crystalline
fraction would otherwise hide unresolved early history. The original final
enthalpy is retained exactly, and latent heat remains in that one account.

The controls bound trial count and timestep size; accepted and rejected trials
both spend that budget. Local caloric/enthalpy-stage overshoots can request
refinement, while source, state, temperature-domain and cancellation errors
retain their refusal. A zero-growth coarse/refined pair with active endpoint
kinetics requests refinement instead of accepting false arrest. Minimum-step,
represented-time or work exhaustion returns no partial accepted interval.

`CrystallizationSimulation::with_adaptive_phase_integration` selects those
controls per original material, including existing cells and scheduled births.
The existing heat-transfer owner still chooses the outer step and supplies its
linear enthalpy path. Each material-cell advance has its own bounded local
trial budget; the original cell and outer-step budgets bound their count.
The timeline reports accepted/rejected work and maximum accepted local
discrepancies separately from heat and caloric residuals. Clones preserve the
selection, and any late refusal preserves the original complete timeline.

Local discrepancy estimates do not bound global error, phase-event timing,
spatial error or the outer first-order heat-transfer error. The maximum local
timestep must resolve source-rate features: finite samples cannot detect every
arbitrarily narrow kinetic peak. Thermal/CFL and timestep refinement remain
independent requirements. Focused timeline cases cover analytic adiabatic
recalescence inside one outer step, latent reheating and a stationary plateau,
energy accounting, deterministic replay, budget exhaustion and cancellation;
these integration cases require the complete native thermal dependency graph.

`CrystallizationPopulation::advance_specific_enthalpy_adaptive` applies the same
bounded step-doubling controls to thermally equilibrated, distinct-age cohorts.
Every candidate uses the existing population midpoint/caloric owner and shares
its immutable cohort buffer. Acceptance checks temperature, the common exposure
increment, and every cohort's crystalline mass fraction as well as its aggregate;
a low-mass history cannot hide unresolved conversion in a mass-weighted mean.
The exact final enthalpy, original masses and separate ages survive. Active
endpoint kinetics prevent accepting false zero-growth onset. Failed refinements,
source failures and cancellation return no partial population. Stationary
zero-heat/zero-growth holds preserve accepted caloric residuals bitwise, avoiding
artificial undercooling and renewed growth at melting. This irreversible
population path adds no remelting, spatial-mixing or molecular-memory model;
its local discrepancies carry the same global/source-resolution limitations.
Focused regressions use an analytic cohort recalescence solution, a rare-cohort
error that an aggregate test misses, cooling onset, replay, and rollback.

Layer L3. Shared conservative thermal transport and solid/vapor exchange.
Runtime dependencies: std, fs-la, fs-math, fs-solver. Spatial discretizations and material
constitutive data remain owned by their existing crates.

## Drifting Maxwellian incident molecular moments

`sublimation::molecular::drift::incident_maxwellian` supplies incoming mass,
conditional vector velocity and crossing energy relative to the original
species enthalpy at a locally moving planar surface. Normal and tangential
drift use the same truncated three-dimensional Maxwellian. The caller supplies
partial density, original specific gas constant, gas temperature, relative
bulk velocity and an already-unit outward normal. Empty species yield exactly
zero mass and no conditional moments; unrepresentable nonzero outputs refuse.

Central moments use deterministic Gaussian functions and compensated arithmetic.
Away-flow tails use positive integral-ratio recurrences with checked terminal
bounds and at most 496 backward iterations. Scaled logarithmic amplitudes retain
representable low-probability incidence; dimensional full moments handle cold
inward beams once the omitted tail is below floating-point resolution. There
is no physical Mach cutoff. The signed mass change from zero normal drift is
evaluated independently, preserving small net deposition near saturation even
when the gross beams round equally. Energy includes tangential slip and the
flux-weighted normal second moment; adding bulk kinetic energy to a stationary
beam would not give these moments.

The kernel adds no state, thermochemical datum, film resistance, continuum
stress, velocity-dependent sticking or Knudsen-layer model. Original material
owners must combine these incident moments with emission/reflection and the
laboratory momentum/energy transformation. Focused tests specify independent
quadrature, opposite-half reconstruction, rotations, tiny drift, scaled tails,
cold beams, vacuum and representability refusals.

## Conducting ambient faces and separate physical baths

`ConductiveAmbient` and `ConductionObserver::with_ambient_faces` bind a selected
cell's complete original exposed ambient area to the same finite half-cell
material and declared surface-resistance law as conducting reservoirs:
`G(h) = A / (d/k(h) + R_surface)`. The original ambient and reservoir baths
remain separate and may change temperature through boundary events. Geometry,
reference mass (except explicitly bound planar phase layers), material identity
and coverage stay fixed; changing a bound
face refuses. Zero conductivity insulates exactly. Zero surface resistance
prescribes temperature at the physical face while retaining the material
half-cell resistance.

The explicit observer protocol carries total ambient conductance through the
original boundary heat and CFL owner. Authored htc is preserved, ambient heat
is not reclassified as an applied source, and the entropy owner retains the
exact applied coefficient and actual bath temperature. Observer composition
rejects duplicate channel ownership. Existing advection, flow, electrical,
rendering and solid-feedback wrappers forward this channel to their original
inner observer.

Implicit chains and graphs use the same boundary-aware residual and analytic
conductivity derivative for each bath, including nonzero dk/dh on a latent
plateau with zero dT/dh. Existing unselected channels retain their original
behavior. Source-only implicit observers explicitly refuse ambient overrides
rather than silently dropping physical heat. Source, topology, work-budget
and cancellation failures publish no changed thermal history.

## State-dependent contact in implicit thermal networks

`enthalpy::conduction::implicit::contact::ImplicitContactLaw` supplies nonnegative
area-specific contact resistance and its two specific-enthalpy derivatives on
the original ordered caloric objects. Chain and graph `with_contacts` attachment
requires unique active links and zero authored resistance on those links, so
source and fixed contact resistance cannot be counted twice. Unselected links
retain the original fixed-resistance path. At zero resistance the source may
select finite one-sided tangents; the native Jacobian uses the finite bulk-limited
`-G*(G/A)*dR''/dh_i` limit. Every source query still requires nonnegative resistance,
and isolation still checks the source domain.

The existing endpoint supplies `G=A/(d_a/k_a+R''+d_b/k_b)` to the same conservative
backward-Euler residual. Each Jacobian column includes both its bulk-k response
and `-G²/A * dR''/dh_i`, then divides by that cell's original mass for the native
energy unknown. The resistance-fraction form avoids explicitly squaring G.
Opposite heat rows retain opposite signs. No additional heat source, enthalpy
storage, clock or nonlinear solver is introduced. The existing nonsymmetric
graph solve, backtracking, entropy, energy gates and transactional refusal own
acceptance. Interface geometry, pressure and interfacial storage remain fixed.

## Existing thermal timeline

### Integrated conductive links in the implicit heat solve

Chain and graph `with_link_laws` attach original two-sided `ImplicitLinkLaw`
closures to selected faces. Their nonnegative secant supplies the existing
equal/opposite heat transfer; their complete dQ/dh replaces the bulk-series
tangent in the same energy-coordinate Newton matrix. Ordered caloric objects,
contact coverage and geometry stay explicit. Source-contact overlap and duplicate
ownership refuse. Unselected links retain their original arithmetic. This adds
no thermal store, heat source or time integrator; original entropy, energy,
nonlinear work limits and cancellation still govern publication.

The enthalpy, advection, contact, observation and event APIs retain their
documented control-volume semantics. A cell stores total energy as mass times
specific enthalpy. Every material chart in one balance uses a compatible energy
datum. The shared timeline is re-exported by fs-conduction; consumers must not
implement independent copies of its storage update.

`advance_observed_continuing_limited` continues the accepted observer from an
earlier `advance_observed_limited` interval within one enclosing physical
advance. It uses the same thermal kernel and all original observer hooks,
without calling `begin` again. Coupling owners can split heat and mass exchange
without resetting a source's per-window energy or physical history. The caller
must supply the remaining shared substep allowance and retain the original
deep-clone/rollback contract; continuation does not renew that work budget.

`DepositionSimulation::apply_heat_j` adds signed endpoint heat to those same
cells without advancing time or re-integrating sources and contacts. The
original material law resolves `h_new = h_old + Q/m`, including latent storage.
Same-time scheduled events must be processed first. The method checks actual
cell-relative stored energy against the supplied joules, with separate local
and aggregate absolute-error tolerances; opposite lost increments cannot hide
behind zero net heat. Compensated sums preserve small changes beside large
caloric reference energies. Domain, arithmetic and cancellation failures return
no candidate. The energy audit identifies the input as external to the thermal
subsystem; a coupling owner must account for the corresponding upstream debit.
No mass, geometry, boundary, event cursor or accepted time changes, and no
mechanical-to-heat conversion fraction is inferred. Five focused regressions
cover latent/signed inputs, unchanged clocks/sources, large reference datums,
local/aggregate loss gates and transactional refusal.

### Enthalpy transported by balanced mass flow

The existing upwind graph evaluates transport relative to one canonical
start-state enthalpy per connected flow component. Fixed mass makes the common
reference cancel analytically; its numerical value cannot create heating or
erase a small source by subtraction of gross formation-energy powers. Component
identities are prepared once, and the same relative fluxes feed thermal and
optional entropy balances. Latent enthalpy still moves at equal temperature.
The gross absolute transported-enthalpy diagnostic retains its original datum.

Signed boundary streams retain the low parts of their enthalpy differences
and mass-rate products through the existing compensated energy reduction.
Species-specific formation contributions can therefore cancel before the
remaining stream heat is combined with an independent wrapped source. The
same reduced powers feed explicit stepping and native implicit endpoint
residuals; the separate advection observation contains only boundary enthalpy,
while the thermal account credits all applied heat once.

The original explicit entropy observer reuses those boundary terms and
compensates its per-cell transported-power reconstruction. Stationary-species
mixing likewise retains those small terms before division by temperature, so
recovered stream heat remains classified as advection rather than another
source or a spurious entropy-balance defect.

Focused caloric regressions exercise opposite species formation offsets,
connection-ID permutations and fractional rates, followed by both explicit
heating and implicit endpoint decay. A coupled heater verifies that stream
cancellation neither erases its contribution nor duplicates its energy credit.
Entropy-enabled advances and endpoint queries preserve the corresponding
source classification and reference-independent entropy balance.

The existing cell mass-balance admission remains 64 machine epsilons of the
incoming/outgoing rate scale. Inside that allowance the equations enforce zero
continuity; the time-integrated absolute cell mass residual and its
`|h_reference * mass_residual|` energy equivalent are reported explicitly as
roundoff diagnostics. Neither is a physical heat source, mass loss, or second
energy debit. G1/G3 regressions cover a flowing calorimeter at large positive
and negative reference offsets, small heating in a uniform circulating loop,
disconnected latent-heat loops, canonical edge orientation, admitted mass
roundoff, and matching energy/entropy reference invariance. This does not add a
pressure, velocity, variable-mass or chemical-reaction closure.

`Advected::with_selective_outlets` explicitly binds a carried composition and
specific enthalpy to existing outlet IDs before stationary species admission.
The original boundary owns the mass rate. On acyclic graphs the same species
solver subtracts those component fluxes, checks availability without clamping,
and propagates the remaining mixed stream. Every selected outlet declares the
same caloric reference as the registered cell and inlet laws. Its carried phase
enthalpy is not inverted through the remaining gas law. The existing thermal
source, continuity, CFL and external-enthalpy accounts include that outlet once;
no separate latent-heat source is invented. Replacing boundary or selective
data invalidates the previous stationary species solution and outlet-law bindings.

`with_selective_outlet_laws` optionally binds every selected outlet to its
original registered `ThermalLaw`, without changing the existing
`SelectiveOutlet` descriptor. Bound laws must have the exact selected
composition/species order, a common enthalpy reference, enthalpy storage and
support for the actual selected h. Entropy additionally requires compatible
source `CaloricEntropyAnchor` values. Energy-only consumers can retain the
original explicit-h descriptors without adding law bindings.

The explicit and implicit entropy owners debit the selected stream's original
entropy at the boundary. Their existing signed control-volume account includes
`mdot*((h_cell-h_selected)/T_cell+s_selected-s_cell)` together with incoming
mixing. Formation-energy terms cancel in each cell's compensated watts before
division by temperature; the native power account remains independently checked.
The selected stream can have a different phase/temperature and lie outside the
remaining mixture's caloric domain. Its own source law supplies h/s, never an
inferred phase reference or a duplicated latent heat. An isothermal separator
can show negative total production: this exposes the missing work/reservoir
closure rather than turning a prescribed selection into a passive mixing law.

The constitutive caller owns selection rates, phase equilibrium and source
alignment. No condensation kinetics or new thermodynamic data are inferred by
this transport primitive. Negative remaining species flux, absent remaining
mixed flow, an unknown/duplicate boundary, a mismatched caloric reference,
selective recirculation and unbound selective entropy refuse. Focused G0/G1
regressions cover component conservation, downstream propagation, phase
enthalpy outside the gas inversion domain, the original thermal energy account,
signed isothermal separation, phase-reference offsets, preserved large
formation references, atomic source/cancellation refusal and these unsupported
conditions. This does not add transient species storage.

### Liquid enthalpy through a stationary solid skeleton

`Advected::with_liquid_carrier` retains the original total-cell `mass*h`
storage while a locally balanced mass graph carries only the liquid phase's
enthalpy. For an admitted pure-substance isothermal latent plateau, liquid
carries the liquid-side plateau enthalpy during coexistence and its original
specific enthalpy when fully liquid. The component reference is evaluated on
that same carried-enthalpy function. Circulation at the melting temperature
therefore preserves stationary solid fractions, while incoming superheat can
melt material through the existing energy account without a second latent
heat source or an independent phase-fraction integrator.

The descriptor requires the complete original phase chart and retained
caloric identity, affine liquid mass fraction on the plateau, equal constant
phase densities and fixed cell mass consistent with density times volume.
The public generic constructor requires the caller to provide the original
complete chart; selected temperature samples do not establish that binding.
Every flowing donor and receiver must contain liquid. Local mass
balance, conservative pair transfers, the existing combined explicit CFL and
implicit energy residual remain authoritative. The carried-enthalpy slope is
zero on the plateau and one in the liquid branch, so the original total-mass
outflow bound remains conservative.

Original `MassBoundary` inlets and outlets can supply balanced open liquid
streams. Every inlet uses that same material and enthalpy datum, lies within
the original source chart and is entirely liquid at or above its liquidus
enthalpy. Entry and discharge cells require liquid support and the same
fixed-density inventory. Outlets carry the current liquid enthalpy, including
in partially solid cells. The original external-enthalpy balance counts the
net stream heat once. Isothermal flushing preserves the stationary solid;
inlet superheat fills its latent-energy deficit before sensible heating.
Boundary and carrier attachment order does not change admission or arithmetic.

The native `MaterialFlowFeedback::with_liquid_carrier` adapter passes the
original pointer-bound phase curve's full knots and retains the descriptor
through explicit graph rebuilds and implicit Stokes endpoints. It requires
the existing physical creeping-flow heating owner. Native velocity represents
liquid volume flux per total face area; homogeneous retained inertia refuses
in both attachment orders. Open ordinary mobile-mushy flow uses the existing
opposed pressure-mouth work owner and explicit pure-liquid enthalpy on actual
inflow, including backflow. Its fixed thermodynamic reference pressure remains
separate from the hydrodynamic pressure gauge: adding hydrodynamic `p/rho` to
the inlet enthalpy would count the physical pressure supply twice. Exact-solid
pressure-mouth topology remains unsupported. Selective streams, stationary species
and nonempty bulk-enthalpy reconstruction still refuse until their own
phase-specific closures are supplied. Alloy melting ranges, unequal phase
densities and moving solid skeletons are outside this model.

The optional explicit-trajectory entropy audit now follows the same moving
liquid. Face and boundary entropy use the original law at
`h_liquid=max(h_bulk,h_liquidus)`, while stored entropy remains the actual
mixture's `mass*integral(dh_bulk/T)`. At pure-substance coexistence,
`s_liquid-s_bulk=(h_liquid-h_bulk)/T_m`; local balanced liquid mass flux
cancels that stationary-phase offset. Isothermal liquid circulation and
flushing therefore create no artificial mixing entropy or latent heat.
The existing component energy reference uses the same carried-enthalpy
function. Source heat/work, spatial upwind mixing and temporal entropy error
remain separate. Instantaneous implicit endpoint queries do not integrate
an entropy history.

Five entropy regressions cover circulation and flushing, independent
superheat-to-latent balances, open-stream entropy, freezing with separately
accounted cooling, large positive/negative caloric datum shifts, fully liquid
equivalence and cancellation.

Five thermal regressions cover isothermal circulation, large datum shifts,
implicit latent melting by superheat, fully liquid equivalence and admission
failures. Two native flow regressions exercise repeated physical wall-heating
epochs and the implicit Stokes energy equations. Five additional open-stream
thermal regressions cover phase-specific discharge, reversed flow, large
enthalpy datums, latent-to-sensible timestep refinement and refusal. Two
pressure-driven flow regressions retain analytic Darcy heat and melt-dependent
resistance across explicit rebuilds and implicit endpoints, including pressure
gauge shifts. These integration tests
require the full native thermal/flow dependency graph to execute.

## Sublimation and deposition

### Deposited coatings on curved substrates

`phase_geometry::RetreatGeometry::{spherical_coating,cylindrical_coating}`
extend the existing mass-owned front with a fixed inert substrate. The phase
mass and `FrontGeometry::volume_m3` describe only the coating; `extent_m` is
its thickness and the exchange area is the outer shell area. Cylinder end
caps remain inert. `occupied_volume_m3` additionally includes the substrate
and is the volume needed for dilute-body/gas-displacement admission.

The geometry follows the shell identities
`V = 4*pi*t*(r*r + r*rc + rc*rc)/3` and
`V = pi*L*t*(2*rc+t)`, where `r=rc+t`. Inversion and recession use factored
forms so a thin deposit is not lost by subtracting nearly equal radii or
total volumes. Each query uses constant storage and bounded scalar work;
there is no measured hardware-speed claim.

The original thermostatted front integrates coating thickness at
`dt_thickness/dt = -j/rho`, refreshes area from the same conserved mass, and
reports extinction exactly at the substrate. After extinction, phase area
is zero and unused time is returned to the caller. This does not remove the
substrate or assign it a replacement boundary law. Existing finite-energy
moving-chamber consumers can use the same geometry without a second front
integrator.

G0/G1/G3 tests in `phase_geometry` cover analytical shell mass/area/tangents,
thicknesses below one ULP of the substrate radius, deposition heat and mass
closure, subdivision invariance and coating-only extinction. These are
uniform constant-density, shape-preserving coatings with an existing seed;
they do not infer nucleation, porosity, substrate thermal storage, curvature
corrections to equilibrium pressure, or gas-volume work.

`sublimation::SublimationLaw` supplies a bounded immutable saturation-pressure,
latent-enthalpy, density and molecular-mass chart. The material crate supplies
reference ice/CO2 correlations and an explicitly parameterized
Clausius–Clapeyron approximation. A numerical correlation is not a sourced
specimen certificate, and supplied density/accommodation/film coefficients
retain their caller-declared applicability.

Positive mass flux means solid becomes vapor; negative flux deposits vapor.
The input pressure is the volatile partial pressure. Outgoing molecules use
surface temperature and incoming molecules use gas temperature. One common
accommodation coefficient preserves detailed balance at equal temperature.
The optional dilute gas film is combined with interface kinetics by series
resistance. Zero accommodation or film coefficient seals the surface.

`surface_flux` returns signed mass and latent-energy fluxes, their analytic
surface-temperature derivatives, and outward-normal interface velocity.
`isothermal_mass_transfer` integrates an ideal, thermostated finite vapor
reservoir exponentially at fixed surface temperature. It uses the shared
deterministic expm1 implementation to preserve short timesteps and does not
overshoot its frozen-temperature equilibrium. Its caller supplies the thermal
work needed to hold the surface at that temperature.

`advance_surface` couples that transfer to the caller's original T(h) chart
with a bounded implicit endpoint solve. The endpoint vapor carries h_end + L.
The energy equation is:

    m_initial * (h_end - h_initial) + transferred_mass * L(T_end) = 0

This is algebraically the balance of remaining solid energy and transported
vapor enthalpy. The latent account is part of vapor energy, not another
independent energy debit. All mass limits are applied before deriving heat,
vapor inventory and recession. A shifted common caloric datum leaves the
physical temperature and mass solution unchanged.

The nonlinear solver requires a continuous nondecreasing caloric chart and a
bracketable endpoint in the declared constitutive domain. It refuses
unsupported temperatures and failed convergence rather than clamping them.
The default budget permits 96 evaluations; callers may select 1 through 4096.
Each evaluation polls the supplied checkpoint. Functions return candidate
values and mutate no caller state on failure.

`sublimation::advance_surface_with_heat` solves the same original chart and
mass-transfer law with explicit signed supplied energy:
`m_initial*(h_end-h_initial) + transferred_mass*L(T_end) = supplied_heat_j`.
Heat can sustain sublimation and cooling can sustain deposition without a
sensible-only intermediate state outside the admitted solid domain. The
receipt includes supplied heat, full vapor enthalpy and their common residual;
the caller debits that heat once. Relative solve tolerance retains the physical
supplied-energy scale when heat nearly cancels the latent transfer.
Zero-duration impulses change only caloric storage, and the original
`advance_surface` remains the zero-heat entry point. Coating/chamber callers
retain their existing zero supplied-heat convention.

Five focused regressions cover independently calculated isothermal
Hertz-Knudsen exchange under heating/cooling, shifted caloric references,
sealed/zero-area/zero-duration limits, heat-driven flux reversal with changing
latent heat, finite-reservoir closure, and invalid/cancelled/exhausted solves.
These are numerical model checks, not new material calibration or a spatial
moving-interface claim.

### Deposition on an initially bare thermal substrate

`sublimation::coating::advance_supported_surface` separately conserves the
substrate and deposited-solid inventories while sharing one interface
temperature. The deposit may start at exactly zero mass when the caller
explicitly admits heterogeneous deposition. Only the deposit can sublime;
exhaustion leaves the substrate intact and permits later regrowth. The
substrate's original T(h) chart and deposit's original h(T) chart remain
authoritative, with independent enthalpy reference levels. The bounded solve
closes `M_s*(h_s1-h_s0) + m_d0*(h_d1-h_d0) + dm*L(T1) = 0`, where positive
dm enters the vapor. The vapor carries the deposit's `h_d1+L`; substrate,
deposit and vapor energy changes close one common account.

`advance_supported_surface_with_heat` includes an explicit signed heat input
in that same substrate/deposit/latent endpoint. The original two caloric charts
retain independent enthalpy references; the substrate is not consumed and the
deposit remains thermally active. Heating can sublime an existing coating,
while cooling can sustain deposition on an exactly bare substrate without a
sensible-only temperature excursion. The phase receipt records supplied heat
once; its full vapor enthalpy already includes latent energy.

`CoatingSimulation::sublimate_cell_with_heat` applies the accepted endpoint to
one actual substrate, coating and independent vapor bath. It checks represented
inventory changes, includes their energy roundoff in the receipt, validates the
original substrate law and commits after the final checkpoint. Exact coating
exhaustion leaves the substrate intact and supports later seed-free regrowth.
The direct operation retains the clock, authored sources and reference geometry;
shared chambers require their joint phase owner. Seven focused regressions
cover coupled heating/cooling, independent datums, finite vapor, extinction and
regrowth, invalid input, bounded work and first/final cancellation.

This is a thin, isothermal deposit approximation with frozen exchange area
and the existing prescribed-temperature reservoir semantics. It supplies no
nucleation barrier, porosity, resolved coating temperature gradient or
contact-resistance model. The coating timeline below applies this transaction
with the original heat-flow owner. Other spatial consumers must explicitly own
their geometry and contacts; whole-cell sublimation remains a separate model.
Six focused tests cover seed-free growth, finite-vapor inventory and energy,
extinction/regrowth, inactive boundaries, independent caloric reference shifts,
invalid charts/inventories, cancellation and bounded nonlinear work.

### Thermally resistive supported deposits

`sublimation::coating::resistance` stores the existing deposit and nonvolatile
substrate at separate temperatures in their original caloric charts.
`heat_resistive_surface` distributes independently supplied heat at fixed
positive masses. `advance_resistive_surface_with_heat` additionally grows or
sublimes the deposit using the original molecular, gas-film and finite-vapor
laws at the deposit temperature. Explicit conductance gives
`q = dt*G*(T_substrate,1 - T_deposit,1)`, positive toward the deposit. The
balances are `M_s*dh_s = Q_s-q` and `m_d0*dh_d + dm*L = Q_d+q`.
The interface heat cancels globally; carried vapor enthalpy includes latent
heat once. A finite nonisothermal gas owner can use the heat-only kernel with
its accepted remaining deposit mass and original stream-energy correction.

Total conservation eliminates substrate enthalpy, leaving one allocation-free
bounded scalar solve. Both local balances and represented total storage are
admitted. Stiff contacts are polished against actual stored/transported heat,
not a large predictor heat; an unrepresentable temperature jump refuses.
Insulated or zero-duration exchange preserves an unheated support exactly.
Original chart refusals, one shared evaluation allowance and cancellation
remain transactional. Bare inputs and complete coating extinction explicitly
require an interface transition; no seed or common-temperature birth is
invented. Conductance and exposed area are caller-owned and frozen over the
step. This two-node, first-order model does not resolve through-layer profiles,
porosity, geometry growth or contact evolution. Independent two-body algebra,
exponential time refinement, nonlinear calorics, latent/finite-vapor balances
and vapor-pressure feedback test the numerical and conservation claims.

### Original nonlinear heat laws on a thermally resistive coating

`sublimation::coating::resistance::ResistiveInterfaceLaw` supplies interfacial
power from the two original material temperatures, positive from support to
deposit. `heat_resistive_surface_with_interface` and
`advance_resistive_surface_with_interface` evaluate that law at the accepted
endpoint in the same separate energy balances as the fixed-conductance APIs.
The phase path retains the original molecular/finite-vapor/gas-film transfer,
complete carried vapor enthalpy and latent heat account.

The interface law must be continuous, nondecreasing in support temperature,
nonincreasing in deposit temperature and passive on its declared domain.
Every queried endpoint enforces finite heat rate, the down-temperature sign,
and exactly zero heat at equal temperatures. Violations abort the shared
bracketer rather than being treated as caloric-domain bounds. These sampled
checks do not prove global monotonicity or constitutive validity. The law owns
its conductivity integral, radiation/contact closure and source context;
internally nontrivial queries also own their bounded work and cancellation.

Total energy still eliminates the support enthalpy, leaving one original
bounded scalar solve with no fitted conductance, effective heat capacity or
second caloric solve. Each endpoint trial spends the same work allowance;
inactive zero-duration contact performs no interface query. Existing authored
constant-conductance entry points keep their multiplication order and exact
zero-contact branch. Both condensed masses remain positive; no bare/extinction
transition or spatially resolved temperature profile is inferred.

Six G0/G1/G3/G4 regressions compare an independent quadratic Kirchhoff solution,
manufactured nonlinear calorics with fourth-power heat transfer, analytical
finite-vapor deposition/sublimation, the original constant-conductance limit,
physical sign/equilibrium refusals, zero-time heat impulses, one trial budget
and cancellation. The complete resistive/geometry scalar lane passes 23 native
tests using the original production owners and real `fs_math`.


### Mass-owned thermal resistance of a growing supported coating

`sublimation::coating::resistance::geometry::ConductiveCoatingGeometry`
joins the original constant-density `RetreatGeometry` to a full-thickness
thermal circuit. Planar, spherical and cylindrical constructors retain a
fixed support; cylinder end caps remain inert. The original geometry owns
the actual deposit mass, shell thickness, exposed outer area, recession and
occupied volume. The support contact area is retained separately, so a
growing curved outer surface cannot silently change the area to which an
authored contact resistance is referred.

For frozen positive conductivity `k`, the layer resistances are
`t/(k*A)`, `t/(4*pi*k*ri*ro)` and `ln(1+t/ri)/(2*pi*k*L)`.
The caller's nonnegative areal contact resistance is divided by the inner
area and added in series; the reciprocal is the conductance consumed by
the original separate-temperature kernel. Spherical resistance uses the
original factored thickness. Cylindrical resistance compensates the
rounding of `1+t/ri` around the shared deterministic logarithm. A positive
film therefore retains resistance even when its thickness is below one ULP
of the substrate radius. Each query is allocation-free with bounded scalar
work; no hardware-speed claim is made.

G0/G1/G3 tests in `resistance/geometry_tests.rs` compare independent planar
and thick-shell Fourier solutions, thin-shell limits, radial subdivision,
and the original two-body caloric solve. Two original phase transactions
verify that actual deposited mass increases the next interval's thickness
and resistance while mass and energy remain with the existing phase owner.
Density mismatch, bare/extinct films, invalid transport properties and
unrepresentable circuits refuse instead of being clipped.

These are quasi-steady full-thickness Fourier resistances feeding a
two-node transient approximation. They do not resolve the deposit's spatial
temperature profile, infer porosity or densification, add thermal expansion,
or supply nucleation, extinction transitions or gas-volume work. Conductivity
and contact resistance are explicit original-law inputs; an interval owner
must refresh them and the mass-owned geometry on each accepted/probed step.
The analytical radial reductions follow MIT Unified Thermodynamics §16.5:
https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node119.html.

### Supported deposits sharing one finite vapor inventory

`sublimation::coating::chamber::advance_supported_chamber` joins bare or coated
substrates to one finite, thermostatted vapor chamber. It reuses the original
supported-surface caloric solve and the same bounded scalar pressure closure as
whole-solid chambers. Every surface sees the common endpoint partial pressure;
the final gas inventory is reconstructed from the accepted signed transfers.
Warm coatings can supply colder bare substrates within the same interval,
without a numerical seed or independent copies of the available vapor.

The original substrate and deposit caloric datums remain distinct. All deposits
must describe the same volatile and a common transported enthalpy datum; initial
and endpoint molar masses are checked. Only deposit mass can sublime. Bare
surfaces contribute no fictitious emission to the pressure predictor. Actual
mass conservation and nonlinear pressure closure remain separate residuals;
latent energy is already part of the transported vapor energy. Compensated
accounts include substrate, deposit and vapor energy. Work is linear in active
surface count per bounded pressure trial, with no dense coupled matrix.

`advance_supported_chamber_with_heat` couples one signed heat input per
supported surface to that common pressure and inventory. Heating and cooling
participate in the same latent/storage endpoint, so the pressure bracket may
expand beyond the initial equilibrium envelope. Each receipt and the aggregate
account record supplied heat once, preserving independent substrate datums and
the common deposit/vapor datum.

If the frozen-temperature pressure predictor has no admissible local caloric
endpoint, the shared solve tries the original pressure, vacuum and the initial
equilibrium envelope as bounded starting anchors before bracketing. Every
refused trial consumes the original pressure budget; cancellation propagates
and no caloric chart or inventory limit is relaxed. A narrow-temperature
manufactured case verifies cooling-driven deposition onto a bare support whose
latent heat makes the joint endpoint admissible despite the failed predictor.

`CoatingSimulation::exchange_chamber_with_heat` accepts sparse heat entries
keyed by persistent cell identity. Every active coated substrate participates,
including bare or unheated receivers; future births cannot consume vapor.
Duplicate, inactive, unknown and nonfinite heat entries refuse. The method
validates original substrate calorics and represented inventory/storage before
one joint publication of all cells, coatings and the shared gas. Per-surface
energy admission uses physical heat, sensible storage and latent scales rather
than arbitrary absolute enthalpy references. Authored sources and elapsed time
remain caller-owned. Four regressions cover an independently manufactured
hot-coating/cold-bare endpoint, datum shifts, zero-heat compatibility, identity
mapping, competing receivers, conservation, budgets and final cancellation.

G0/G1/G3 regressions cover analytical common backpressure with gas-film
resistance, warm-to-bare transfer from an empty chamber, exhaustion, independent
caloric reference shifts, live first-order timestep refinement, surface-order
symmetry, work limits, incompatible molar masses and cancellation. These are
equation and conservation checks, not experimental deposition validation.
The existing thin-coating, fixed-area and prescribed-gas-temperature assumptions
remain; nucleation, porosity, gas sensible heat and resolved through-layer
transport require their own physical owners.

### Supported coatings on the native thermal timeline

`sublimation::coating::heat_supported_surface` applies signed heat at fixed
substrate/deposit masses by solving
`M_s*(h_s1-h_s0) + m_d*(h_d(T1)-h_d(T0)) = Q`.
It reuses the bounded endpoint solver, including backtracking when the
substrate-only outer trial leaves the domain but the combined-storage solution
is admissible. Independent material datums remain intact. Unrepresentable heat
is either explicitly bounded by the declared residual tolerance or refused.

`enthalpy::sublimation::coating::CoatingSimulation` owns this storage together
with the original thermal timeline and supported phase operator. Substrate
cells retain their original material, mass, enthalpy, geometry and identity;
deposited mass and its enthalpy remain separate, including exact zero inventory.
The wrapper's mass/energy queries include BOTH condensed inventories. Each
selected surface explicitly supplies its original pure-solid h(T), compatible
phase law, kinetics and independent prescribed/thermostatted vapor reservoir.
Future births may carry an authored coating; its mass and sensible energy enter
the incoming account exactly once.

`CoatingSimulation::new_shared_chamber` instead selects one common volatile
law, one original deposit h(T), and one finite chamber for all authored
`ChamberCoatingSurface`s. `advance_chamber_to` supplies explicit surface and
pressure budgets and calls the joint supported-deposit operator on each native
thermal substep. Active cells are ordered by persistent identity; future cells
join only at their scheduled birth, with the then-current chamber inventory.
`shared_reservoir` is the single gas stock. Per-coating reservoir observations
repeat it and must not be summed. Individual reservoir replacement refuses in
this mode; `set_shared_reservoir` replaces the external bath atomically.

`CoatingSimulation::advance_coupled_to` and
`advance_coupled_chamber_to` consume each native substep's authored heat in
the same supported phase endpoint. Ordinary cells retain their native thermal
update. Coated substrates bypass the heat-only intermediate and receive their
net source, bath and contact heat exactly once, together with substrate/deposit
sensible storage and latent transfer. The shared variant couples all active
coatings through the original common-pressure solver; unheated bare receivers
remain active.

These advances reuse the existing heat/event loop, start exposure, explicit
CFL and whole-call budget. Future births join at the event endpoint and carry
their original coating inventory; endpoint coverage cannot rewrite the
preceding exchange. Accepted native endpoints pass their original substrate
validators and finite-storage gates. One private candidate contains the entire
timeline and gas inventory. `CoupledCoatingAudit` reports complete storage,
external heat, carried vapor energy and actual mass/energy residuals without
a fictitious heat-only intermediate account. Local compensated conservation
checks and represented-inventory corrections retain their original owner.

Three coupled-timeline regressions exercise strong heating and cooling whose
split intermediate leaves the caloric domain, common vapor with native contact
heat and an unheated bare receiver, and both modes' exact birth/coverage
intervals, whole-call budgets and final cancellation.

The material adapter `fs_material::thermal::PureSolidThermalLaw::
specific_enthalpy_at_temperature` supplies original solid h(T) directly from
its borrowed curve. It retains exact source knots, varying sensible slopes and
the melting-onset endpoint, with no constant-Cp surrogate or extrapolation.
Ambiguous all-solid isothermal states and unrepresentable inverse queries
refuse; mixed/liquid states never extend the solid branch's support.

Native heat links, ambient and reservoir exchange, source power, contacts,
separation and boundary events retain their original owner. Heat storage
includes the coating before the following phase update on every CFL substep.
The substrate's existing CFL bound remains conservative for continuous
nondecreasing deposit h(T). Endpoint coverage changes cannot retroactively
alter deposition, and future cells have no pre-birth exposure. Deposit
exhaustion preserves the substrate and its subsequent ability to collect vapor.

The complete candidate includes event time, substrate enthalpy, deposit mass
and vapor inventory. Refusal or cancellation publishes none of them. One
full-advance thermal budget and bounded local nonlinear solves apply. Local
energy differences and compensated sums preserve caloric reference handling;
actual mass-rounding corrections are charged to the energy residual, and
opposing local defects cannot cancel out of its absolute admission.

The split APIs retain first-order heat/phase splitting; the coupled APIs use
first-order explicit heat transport and implicit phase endpoints. Both retain
fixed reference area and conductances and a shared substrate/coating
temperature. They do not resolve thick-layer resistance, porosity, melting,
moving geometry, dynamic
observer-based sources, or evolving gas sensible heat. Finite baths supplied to
the independent constructor remain independent. Seven native integration tests cover combined
storage/conduction/source work, zero-inventory growth and regrowth, finite vapor
conservation, exact event exposure, coated births, cancellation/budgets,
restart and independent caloric datums; seven heat tests additionally cover
nonlinear storage, domain backtracking and unrepresentable heat.
Five shared-timeline regressions additionally cover competing bare substrates,
canonical replay, late births without replenishing vapor, empty pre-birth
intervals, contact shielding, combined heat storage, independent substrate
datums, shared work limits and late cancellation. Their selected-source native
lane exercises the original timeline and phase bodies; full workspace gates
remain a separate dependency requirement.

## Numerical and physical limits

The coupled endpoint is a first-order reduced model. Area remains frozen over
one phase step; recession uses endpoint solid density. Finite gas temperature
is prescribed: its thermostat owns subsequent sensible-heat equilibration.
The dilute film does not model non-dilute Stefan flow, pore resistance, plumes,
curvature, nucleation, melt competition, radiation or unresolved gradients.
Those mechanisms require explicit compatible models, not hidden coefficients.

For a whole-solid body, returned zero mass denotes extinction. Its spatial
consumer must update topology, contacts and conservative remapping or refuse
that endpoint. A supported coating can instead reach zero inventory while
retaining the substrate and its topology, as described above.
Reporting recession does not itself move a mesh, alter a level set or change a
flow/electrical geometry map. Uniform-cell use requires a justified internal
temperature approximation. Refining time cannot repair an invalid spatial
model.

The scalar flux and finite-reservoir kernels allocate no memory and perform
constant work. The implicit solver uses constant memory and at most its
declared evaluation budget. This is an algorithmic statement, not a measured
hardware throughput or cross-ISA reproducibility claim. Transcendentals use
the existing fs-math owner; current executed evidence is same-ISA.

## Focused verification

G0/G1/G3 tests cover signed deposition, equilibrium, unequal gas temperature,
analytic tangents, kinetic/film limits, finite-reservoir exact relaxation and
semigroup behavior, tiny timesteps, inventory limits, implicit cooling and
heating, mass/energy/recession consistency, caloric-datum shifts, first-order
refinement, cancellation and exhausted nonlinear budgets. Material tests
compare the actual correlations to published reference values and their
analytic derivatives to finite differences. These tests establish the
declared reduced model; they are not experimental validation of a complete
sublimating object or a general spatial Stefan solution.


## Phase exchange on the existing heat timeline

`enthalpy::sublimation` composes the surface operator with the original
`DepositionSimulation`. `sublimate_cell` atomically changes the active cell's
actual mass and specific enthalpy plus its vapor inventory, without advancing
the clock. `advance_sublimating_to` performs bounded first-order heat/phase
splitting with one total thermal work budget and one global mass/energy audit.
It honors birth and boundary/contact times exactly, uses the exposed area from
the start of each interval, and never sublimates a not-yet-born cell. The
existing caloric chart remains authoritative; rejected out-of-domain trial
points do not invalidate a successfully converged in-domain endpoint. All
failures and cancellation preserve the original timeline and reservoir.

`sublimate_cell_with_heat` supplies a signed energy input to that same atomic
cell/reservoir operation. It solves original solid storage plus complete vapor
enthalpy against the supplied heat at one endpoint, allowing heat-driven
sublimation and cooling-driven deposition without an out-of-domain
sensible-only intermediate state. The caller owns the source debit and duration;
this direct operation neither advances the clock nor integrates authored heat,
convection or contacts. The zero-heat method and split timeline retain their
existing composition. Native regressions cover both signs of latent transfer,
unchanged authored boundary data, invalid heat and final-checkpoint rollback.

`advance_coupled_sublimating_to` composes native heat transport and the
implicit solid/vapor endpoint on the same thermal substeps. Each step computes
the original source, convection, bath and contact watts once; its selected
solid consumes that signed heat in the latent solve, and every other cell
receives the original explicit update. Internal contact heat cancels between
its two owners. The original conductive CFL is recomputed after mass changes,
and one total native substep budget spans all births and boundary events.
No sensible-only phase intermediate is required. The new
`CoupledSublimationAudit` reports complete open-system balances separately
from the existing split operator's thermal subaccount.

Both timeline routes use the same compensated relative-enthalpy account.
The coupled route retains first-order explicit heat quadrature and implicit
local phase feedback; nonlinear convergence does not establish temporal or
spatial accuracy. Three native regressions cover authored heating/cooling,
combined source/bath/contact exchange, birth/coverage timing, work exhaustion
and final cancellation. Gas temperature, geometry and pure-solid admission
retain their existing ownership.

The whole-window audit uses compensated per-cell changes relative to each
cell's starting or birth enthalpy, together with vapor enthalpy on that same
relative reference. Mass conservation is checked separately. Unchanged large
stored energies therefore do not erase sublimation cooling or deposition
heating; full original-datum storage and vapor diagnostics remain available.
Maximum phase intervals use integer-indexed endpoints and restart at birth or
boundary events, avoiding a fictitious terminal phase solve from accumulated
clock roundoff. Native conductive substeps retain their original budget.
Regressions compare sublimation and deposition with an insulated high-energy
neighbor against isolated runs, and check decimal-step and birth/cover timing.

The whole cell must remain a pure solid on its original caloric branch.
Temperature alone cannot establish solid fraction on a melting plateau;
mixed-phase inventories require a separately owned partition. The supplied
finite reservoir belongs to this one surface and remains thermostated. A
joint chamber shared by several surfaces needs a joint backpressure solve.
Reference geometry/contact areas remain fixed and recession is reported.
Complete cell exhaustion refuses until a consumer owns topology removal.
There is no implied moving mesh, gas-energy solve or material calibration.

Vapor energy includes the solid's original sensible-enthalpy datum plus latent
heat; the separately reported latent account is never an additional debit.
Ten new native integration regressions cover balances, frost growth, covered
surfaces, birth/contact timing, nonlinear chart trials, transactional refusal,
whole-advance budgets and first-order refinement. The combined exact-source
thermal/material suite passed 95 tests with one existing performance test
ignored. `cargo run -p fs-material --example sublimating_ice` demonstrates the
real ice law with explicit illustrative calorimetry in a maintained vacuum and
a finite chamber; both native example cases passed mass/energy assertions.


## Several surfaces sharing one vapor chamber

`sublimation::chamber::advance_chamber` solves one shared endpoint vapor
pressure together with the original conservative surface enthalpy operators.
Its backward-Euler mass balance admits warm-to-cold redistribution, including
condensation of vapor emitted during the same interval into an initially empty
chamber. Depositing surfaces never receive independent copies of the original
vapor inventory. One volatile species with constant molar mass, pure-solid
caloric branches and compatible enthalpy references are required. Gas volume
and temperature stay prescribed.

`advance_chamber_with_heat` accepts one finite signed heat amount per input
surface. Each original caloric balance includes its own heat once while all
surfaces retain the same solved vapor pressure. The joint receipt exposes the
compensated total supplied heat; local residuals already subtract their heat.
The pressure bracket can expand above the initial-temperature equilibrium
envelope when heating raises the final saturation pressure.

`DepositionSimulation::sublimate_chamber_with_heat` applies that joint endpoint
to actual cells and one vapor inventory. Heat follows caller cell identity
through canonical sorting, and unheated/actively cooled surfaces continue to
participate in the shared balance. Every candidate material and the final
checkpoint are admitted before any cell or reservoir write. The direct
operation leaves the clock and authored heat sources unchanged. Six focused
regressions cover hot-to-cold redistribution, pressure-envelope expansion,
caloric reference shifts, input validation, cell-order mapping and rollback.

A frozen-temperature predictor and safeguarded scalar pressure interpolation
avoid a dense coupled matrix; each trial costs O(surfaces), with explicit
surface, local-iteration and total pressure-trial ceilings. Unsupported outer
trial bounds are backed away from. The same accepted transfers reconstruct
final vapor mass and drive every solid mass, latent/sensible-energy and
recession account. The nonlinear pressure/mass closure residual is reported
separately from actual published mass conservation. Refining time controls the
first-order integration error; the residual does not certify model accuracy.

Eight new native regressions check independent simultaneous algebra, transient
refinement, stationary-pressure circulation, initially empty chambers, shared
inventory, area partition/order, refusal and finite-capacity energy closure.
The combined thermal/material suite passed 113 tests with one existing benchmark
ignored. Thermally evolving gas, molecular heat accommodation, mixed phases and
mesh/topology changes remain separate physical owners.

### Relative entropy on the shared thermal timeline

ThermalLaw may provide its exact relative integral of dh/T at the original
fixed-pressure/composition convention. Missing entropy remains unavailable;
no substitute heat capacity or entropy reference is invented. The opt-in
advance_to_with_entropy and BackwardEulerStep::evaluate_with_entropy use the
original accepted thermal states, actual interface coefficients and the same
start/endpoint flux quadrature. EntropyAudit separates exact material storage,
actual ambient/support reservoir entropy, conductive production, unclassified
source energy divided by cell temperature, and the numerical balance defect.
The source flag survives cancellation of signed sources. Birth entropy is
excluded on both sides; subsequent storage starts from actual birth enthalpy.
A trial audit is accepted only with its trial. Ordinary energy-only APIs do no
entropy integration. G1/G3 tests cover isolated exchange, refinement, reservoirs,
latent storage, equilibrium, deferred birth, resolved implicit sources, missing
material entropy and cancellation without publication.

`advance_advected_with_entropy` composes the original explicit `Advected`
entropy owner with this thermal account. Net stream entropy, upwind mixing and
advective balance residue populate `EntropyAudit::advection`; they are excluded
from unclassified source energy/T. Thermal observation sits inside the flow
wrapper, so actual nonadvective source flags survive signed cancellation without
subtracting large transport accounts. Both owners use the same accepted states,
CFL/event grid and original source/history hooks. The returned observer retains
its inner history and enabled flow entropy account, with the usual per-advance
`begin()` reset semantics. Balanced same-law bulk, admitted liquid-carrier flow,
and stationary mixtures with original source entropy anchors are covered.
Selective streams require their own original outgoing-law bindings; their
signed entropy contribution does not supply separator work or phase equilibrium.
G1/G3/G4 regressions cover actual bath temperatures,
datum shifts, latent storage, unchanged heat trajectories, retained histories,
events, unavailable entropy and cancellation without publication.

Stationary mixture entropy uses `ThermalLaw::advected_entropy_anchor` and the
same law's exact `dh/T` integral. `CaloricEntropyAnchor` supplies an original
enthalpy and reduced specific entropy at that state, plus a shared source
convention identity. Reduced entropy retains chemical mixing, pressure and phase
terms but analytically omits one common set of linear species reference
constants. Frozen storage and admitted local species balance cancel those
constants; subtracting large rounded absolute source entropies is not an
admissible replacement. Missing anchors, incompatible conventions and nonfinite
source values refuse before publishing a state. The source coordinate is
evaluated once per participating cell or inlet, then reused across its faces.

Unlike-composition/pressure mixing terms remain signed accounting contributions;
same-law per-face nonnegativity does not generalize without physical pressure
and work assumptions. Formation-enthalpy heat contributions are compensated
per receiving cell before division by its temperature, separately from the
source entropy sum. This preserves small isothermal mixing entropy beside large
positive or negative formation energies. The native transported-power account
remains independent, leaving its finite continuity/arithmetic difference in the
reported residue. Analytic ideal-mixture regressions cover nonzero mixing at
matched temperatures, unequal-temperature mixing, recirculation and reversed
flow orientation, signed pressure terms, source ambiguity and cancellation.


## Shared chambers on actual multi-cell heat histories

`enthalpy::sublimation::chamber` applies one joint phase exchange to canonical
cell identities using each cell's original caloric law and exposed area.
`sublimate_chamber` validates every candidate before committing any solid or
vapor state. `advance_chamber_sublimating_to` composes original heat transport
with joint exchange, honors start-of-interval exposure and exact scheduled
birth/contact boundaries, and retains one whole-advance thermal work budget.
The report checks combined solid/vapor mass and the original total energy
account, while keeping each surface's recession separate.

Whole cells must remain pure solid; gas temperature and reference geometry
remain fixed. Complete extinction requires a separate topology owner and
refuses atomically. Ten native integration regressions passed for inventory,
redistribution, original caloric charts and heat graph, canonical order,
event timing, work ceilings and late cancellation. The runnable fs-material
example `ice_cold_trap` uses the real ice law to compute frost growth on a cold
film supplied by a warmer film through an initially empty shared chamber.


## Molecular counterflux for nonisothermal interfaces

`sublimation::molecular` returns both nonnegative half-space mass fluxes and
one signed net flux from the original bounded constitutive and kinetic owners.
The equal-temperature path retains stable pressure-difference evaluation.
Different surface and local vapor temperatures remain distinct; volatile
partial pressure is never replaced by total carrier pressure. Scalar work is
constant and allocation-free. Four native regressions cover equilibrium,
unequal temperatures, pressure semantics and refusal boundaries.

An energy consumer must retain the opposing streams: zero net mass does not
imply zero heat exchange. For Maxwellian ideal vapor with equilibrium internal
modes, each stream carries h(T)-R_specific*T/2 in an explicitly aligned energy
reference. The material gas adapter owns that caloric evaluation and both
endpoint energy accounts. This interface assumes energy/angle-independent
equal sticking and specular nonsticking reflection; it supplies no gas-film,
diffuse thermal accommodation, Knudsen-layer or moving-geometry model.


### Integrated turnover through supported-film birth and extinction

`sublimation::molecular::turnover::integrate_molecular_turnover` integrates the
original kinetic coefficients at frozen surface/gas temperatures and exposed
area. The finite ideal vapor follows `dm_v/dt = c - lambda*m_v`; both gross
streams are integrated analytically. Molecules may return during the same
interval, so neither stream is limited to a fraction of its original inventory.
A prescribed bath uses the same original coefficients and sign convention.

The caller explicitly admits barrier-free heterogeneous deposition and a
zero-residence bare-interface closure. At zero coating mass, outgoing sticking
molecules are limited to the smaller of the saturation capacity and the
incoming stream. An undersaturated bare support immediately reemits incoming
sticking molecules at its own temperature; a supersaturated support begins
accumulating coating. The substrate never supplies volatile mass. Complete
film loss is located analytically, and equal incoming/reemitted streams continue
for the remaining bare interval. Those streams can exchange heat even at zero
net mass. No arbitrary seed or fractional-inventory termination threshold is used.

`MolecularTurnover` retains the signed transfer independently of subtracting
large gross streams, the actual endpoint inventories, first extinction time,
and separate compensated stream/deposit/vapor mass residuals. Tiny return
streams use a cancellation-free exponential ratio; independently representable
late tails survive exponential underflow. Consumers must admit the residuals
against their physical mass budgets and retain or refuse topology differences
caused by inventory encoding. Zero-time, sealed and exact-equilibrium states
preserve inventories exactly. Scalar work is allocation-free and bounded; there
is no nonlinear iteration or hardware-performance claim.

Eight G0/G1/G3 native tests cover the affine solution, gross recycling, birth,
extinction and subdivision across the bare boundary, the prescribed film limit,
tiny return streams/events, exact stationary/no-op states, late tails and invalid
inputs. This operation supplies no caloric evolution: its consumer must close
both original energy accounts using the distinct stream temperatures. Fixed
temperature and area, bulk-solid saturation, complete sticking-stream thermalization
and the original dilute film remain approximations. No nucleation barrier,
submonolayer coverage, residence kinetics, morphology or resolved film is claimed.

### Molecular normal momentum and explicit Maxwell reflection

`sublimation::molecular::momentum` supplies allocation-free half-space momentum
moments for the SAME admitted gross mass streams. The crossing-molecule mean
normal speed is `sqrt(pi*R_specific*T/2)`; emission and sticking loads add even
when net phase transfer is zero. Reflected nonsticking molecules are a distinct
account, including carrier species. Explicit Maxwell diffuse fractions retain
incident mass and reemit that fraction at the surface temperature; the rest
reflects specularly. A thermal accommodation measurement alone does not select
an angular scattering model. At equal temperatures the complete load recovers
the original gas pressure independently of the sticking probability.

Finite input, probability, overflow and nonzero-underflow admission applies.
These moments add neither heat nor a second ambient-pressure load. Their
stationary surface and zero-drift local Maxwellian assumptions do not supply
moving-wall, continuum-stress or gas bulk-momentum evolution. Six regressions
cover independent Gaussian quadrature, equilibrium and sealed limits, vacuum
recoil, zero-net counterflow, Maxwell reflection and extreme-scale refusal.


## Original thermal storage for external solid/gas coupling

`enthalpy::sublimation::external` exposes the original pending events, work
budget and cell caloric evaluation, and transactionally applies prepared solid
mass/enthalpy endpoints to the existing thermal storage. Every full expected
cell, positive remaining inventory, caller-admitted pure-solid h interval,
original caloric endpoint and original-geometry exposure is checked before a
single write. Endpoint contact events may change current exposure; they cannot
retroactively alter the preceding interval's declared area. A surrounding
coupled owner mutates its private trial without another graph clone per phase
step. Joint gas/solid conservation and physical phase admission remain with
that owner; this boundary does not invent calorimetry or geometry evolution.

The material gas timeline now composes this boundary with the actual original
heat operator and one jointly evolving gas. Four native owner tests cover
state preservation, stale/bad batches, endpoint covering and late cancellation.
Both existing isothermal sublimation timelines also jump directly to the next
event/target while the WHOLE graph is empty, avoiding unbounded tiny pre-birth
intervals. Two regressions prove a long pre-birth gap with zero thermal work
allowance. The native sublimation suite passed all 51 selected tests.

## Heat-driven condensed-fuel pyrolysis

`enthalpy::pyrolysis::PyrolysisSimulation` extends the original activated
thermal cells with separate virgin and char inventories. It retains their
identity, geometry, contacts, boundary conditions and extensive storage, and
uses the same conductive/boundary heat operator through a private caloric
callback. No second conduction graph or alternate transport algebra is used.
`CondensedCaloric` adds source-owned h(T) to the original `ThermalLaw::T(h)`;
admission checks that the charts agree. Virgin and char may have distinct
nonlinear sensible calorics. `VolatileEnthalpy` supplies the escaping gas's
physical enthalpy on a compatible formation reference.

`new_with_thermal_materials` additionally admits explicitly named heat-only
materials in the same activated graph. Finite molds, walls and phase-change
heat sinks retain their original mass, geometry, `ThermalLaw::T(e)` and cached
thermal slope bound. They receive no virgin/char inventory and release no
volatile material. The ordinary `new` constructor still requires reactive
models for every active material; missing kinetics never silently disables a
reaction. Duplicate, unknown and simultaneously reactive/heat-only declarations
refuse. Reacting virgin and char branches must declare enthalpy storage for the
open volatile balance. Closed-mass heat-only bodies may retain internal energy
or their original declared sensible-storage convention.

Both explicit and implicit transport retain the same paired heat, boundary and
stored-energy accounts across these materials. The implicit local query uses
the heat-only body's original energy inverse directly, without a fictitious
reaction or an apparent heat capacity. Zero thermal slope is retained on a
sharp latent plateau; the nonlinear residual still decides the accepted heat.
Focused G1 tests independently solve two-body explicit and backward-Euler
calorimetry, including finite-capacity and isothermal latent sinks, and reconcile
remaining reactive storage plus actual outgoing volatile enthalpy. Domain
failures and cancellation preserve both original histories. The homogeneous
receding-slab owner refuses graphs containing heat-only bodies before indexing
reactive inventory; moving mixed-material geometry is not inferred.

`with_conductivity` binds selected original links to their existing
`ConductiveLink` half-distances/contact resistance and `ConductivityLaw`
sources. Heat-only bodies query their original stored-energy coordinate,
including k(h) variation while temperature is on a latent plateau. Reactive
cells query each present constituent's ORIGINAL h(T) at the common solved
temperature, then apply the explicitly selected Series or Parallel arrangement
using actual virgin/char volumes and caller-sourced constant constituent
densities. Absent components make no conductivity query; a present zero-k
component insulates a series arrangement. The slab and fixed-geometry paths
share the same mixing arithmetic.

Conductivities resolve once per participating cell at each explicit substep or
implicit candidate. Original series/contact arithmetic supplies G, and the
existing heat owner applies the same equal-and-opposite pair transfer. Newton
uses the existing one-cell reactive probes for both temperature and conductivity
response, including the dG/d(absorbed-energy) contribution to the paired flux.
It does not rerun a complete graph/reaction solve for every Jacobian column.
The original no-attachment constant-conductance path is retained. Adaptive
implicit time control automatically includes the attached nonlinear transport.

This attachment owns fixed reference face metrics only. Constituent volume
fractions do not infer pore-gas transport, porosity, swelling or recession from
lost volatile mass; an effective porous-char conductivity must be admitted by
its original source. Bound links reject manual coefficient overwrites, and
moving-slab adoption refuses the fixed attachment. Original source domains,
heat/mass accounts, cancellation and refusal remain authoritative. Focused
G1 tests cover independently solved nonlinear and char-insulating endpoints,
continuous char-mixture and latent-buffer trajectories, explicit/implicit
refinement, absent-source/zero-k behavior and transactional workspace reuse.

One first-order non-oxidative reaction converts virgin material into an explicit
char mass yield and outgoing volatiles. Its accepted endpoint temperature
determines the backward-Euler inventory, retained char and outgoing enthalpy:

```text
k(T) = A exp(-E/(R T))
m_virgin_new = m_virgin_old / (1 + dt k(T_new))
extent = m_virgin_old - m_virgin_new
mass_out = (1 - char_yield) extent
energy_out = mass_out h_volatile(T_new)
```

The local energy equation balances the change in remaining solid storage plus
that exact outgoing enthalpy against supplied heat. Reaction enthalpy is already
in the compatible material references; it is never a second source or debit.
The independent virgin/char caloric changes and reaction enthalpy differences
are evaluated before gross formation energies are combined. Per-cell releases
retain their integrated mass, enthalpy and consumed virgin mass for downstream
gas coupling, even when emission temperatures change within one advance.

Tiny rate-duration products use stable loss fractions. Published consumption
is reconciled to actual retained inventory, so an unrepresentable decrement
cannot emit fictitious mass or heat. The absolute kinetic mass-roundoff and
reaction-energy equivalents are reported separately and checked against the
caller budgets; they are numerical discrepancies, not physical transfer streams.
Local, cumulative absolute and whole-advance mass/energy residuals cannot be
hidden by cancellation between cells. Errors and cancellation publish no
candidate. Boundary and conductance changes use original transactional admission.

This initial lane requires endothermic single-step kinetics, monotone sensible
caloric charts, explicit temperature/source domains and fully activated fixed
topology. It refuses pending geometry events, exhausted work budgets and
unrepresentable depletion. It does not seed residual solid mass. Conduction
retains the original explicit stability bound using the current condensed mass
and mixture capacity; local reaction is implicit in its final temperature.
Near depletion this fixed conductive mesh can become stiff. Moving surfaces,
porosity/permeability, pore-gas transport, oxidation, latent-phase state,
multi-step wood chemistry and fire-spread closure remain separate capabilities.
The coupled method is first order, with no unconditional stability or entropy
theorem. Conductivity changes require an explicit caller-supplied law.

Eleven focused regressions cover independent isothermal and adiabatic
solutions, distinct char calorics, first-order temperature/mass refinement,
equivalence to original heat transport with reaction disabled, spatial heat-
driven mass release, caloric-datum shifts, small-rate representation, and
domain/cancellation/work-budget refusal. The native full thermal library suite
passes 164 tests with one existing performance test ignored. This is focused
native Cargo evidence, not a repository DSR gate or experimental validation.

`heated_pyrolysis_slab` exercises the actual 1D graph with a resolved front
half-cell heat balance, irradiation, emission and convection. It uses retained
NIST/MaCFP scalar PMMA inputs in an authored heating scenario; its optional
20% char sensitivity is explicitly not a wood material fit. `macfp_tga` uses
the same kinetic update against the retained MaCFP FDS material-temperature
history. Source data, units, calibration/verification distinctions and commands
are in `data/reference/compartment_fire/pyrolysis/README.md`.

### Implicit conduction and finite-rate decomposition

`pyrolysis::implicit` removes the explicit conduction stability limit while
retaining the original physical model. Its nonlinear unknown is absorbed cell
energy. Every trial calls the existing `solve_cell` for temperature, virgin/char
inventory and escaping volatile enthalpy, then the original
`heat_rates_with_calorics` evaluates endpoint conductive/Robin/source powers.
The actual residual is absorbed energy minus duration times endpoint power.
No alternate heat graph, heat capacity, kinetic law or heat-of-reaction source
is introduced.

`implicit_workspace` prepares original dense LU and numerical buffers for at
most 128 cells. The Jacobian combines the original link/boundary conductances
with local caloric-response probes, avoiding a complete graph reevaluation for
every column. Original bounded heat-rate maps and candidate/result objects
still allocate; this is not an allocation-free whole-advance claim.
`advance_implicit_to` returns the unchanged original release/audit format plus
actual nonlinear, local-solve, linear-solve and backtracking work counts. Local
and whole-advance conservation gates include original kinetic roundoff and the
separately reported graph heat-transfer summation discrepancy. Cancellation or
nonlinear/domain refusal publishes no candidate.

Seven focused native tests cover an independent stiff two-capacity endpoint,
simultaneous heating/conduction/reaction, char and caloric-datum invariance,
actual Robin coverage, temporal refinement, activated Arrhenius comparison
with the original fine explicit method, and refusal/cancellation/reuse. The
stiff case completes in one implicit interval while the explicit owner exhausts
its four-substep budget. This remains first-order backward Euler with caller
time-resolution controls; removing the conduction CFL is not an unconditional
nonlinear-convergence, temporal-accuracy or entropy theorem. Geometry and
authored coefficients stay fixed within each original advance.

`advance_adaptive_implicit_to` adds bounded local time-error control to that
same fixed-geometry model. Each comparison solves one native full interval and
two native half intervals with the reusable workspace. Only the two original
fine states and their mass/heat/volatile receipts are accepted; no composition
or enthalpy extrapolation can create negative inventory or skip a latent
plateau. Accuracy rejection halves the unpublished comparison. Nonlinear,
material-domain, conservation and cancellation failures propagate immediately.

The error estimate includes original specific storage, temperature, separate
virgin and char inventories, and each cell's outgoing volatile mass and
enthalpy. Storage participates even when a finite latent buffer stays at one
temperature. Relative tolerances scale changes during the interval rather
than gross formation enthalpies. Outgoing-energy comparisons use actual native
release masses and the original volatile law at each accepted endpoint,
subtracting a common cell reference before multiplication. This avoids
subtracting rounded extensive formation energies; accuracy remains limited by
the precision of the original source charts. These are local temporal
discrepancies, not certified whole-trajectory errors or an entropy theorem.

Whole-call conservation tolerances and the original accepted-substep ceiling
are never renewed between comparisons. `maximum_comparisons` bounds all three
native solves, including rejected work; each original nonlinear solve retains
its iteration/backtracking/scalar bounds. Reports separate accepted physical
accounts from all attempted numerical work and return the next suggested span.
Tests compare independent exponential virgin/char decay and an activated
Arrhenius reaction exchanging reciprocal heat with a latent buffer, verify
exact equivalence to native accepted fine steps after rejection, and exercise
datum shifts, budgets, source refusal, late cancellation and workspace reuse.

### Conservative noncharring slab recession

`pyrolysis::slab::RecedingSlab::new` adopts the original homogeneous x-aligned cell
chain with caller-sourced constant density, conductivity and cross-section.
Actual retained mass sets thickness through `L = m / (rho A)` with a fixed
back face. This is the single-material density/length relation in the pinned
FDS Solid Chapter's shrinking/swelling section; fixed-back anchoring and
conservative rezoning are explicit choices of this implementation.

Each bounded interval uses the original implicit conduction/reaction owner.
Excessive actual per-cell mass loss retries the complete unpublished interval,
without accumulating rejected volatile releases. Accepted specific enthalpy is
reconstructed with limited linear slopes in retained-mass coordinate and
integrated over overlaps onto equal-mass cells. The original caloric inverse
recovers temperature; original cell geometry and series conductances are rebuilt.
Mass and extensive enthalpy are transferred, not averaged temperatures. Cell IDs
are spatial slots after rezoning rather than tracked material particles.

Original volatile mass/enthalpy integrals are unchanged and retain their
production slots and homogeneous source binding. The original aggregate audit
includes actual remap storage/discrepancy once; separate mass/energy residuals,
rounded-geometry mass equivalent and enthalpy redistribution are returned.
The redistribution diagnostic is not lost heat or a PDE error bound. No seed
mass, heat-capacity floor or release clipping extends the physical domain.
Geometry and coefficient evolution remain first-order split from each implicit
interval. The charring extension below supplies separate residue geometry.
Pores, mechanical work, thermal expansion, gas-domain volume displacement and
complete burnout require separate models.

Seven native groups verify independent isothermal inventory/recession,
nonlinear-caloric overlap conservation, monotone reconstruction/refinement,
formation-datum shifts, reaction-disabled identity, complete retries and
cancellation/replay/source/work refusals. The combined native thermal suite
passes 178 tests with one existing performance test ignored; DSR is unavailable
in this execution environment, so this is not a repository-level gate result.

`heated_pyrolysis_slab` retains its original five arguments and fixed default;
an optional sixth `receding` argument enables the new noncharring path with the
actual front half-cell resistance. Original NIST/MaCFP scalar inputs and the
effective volatile enthalpy remain unchanged. The 120-second 24/48-cell runs
complete, with emitted-mass spatial sensitivity of about 1.88%; their actual
accounts are retained in
`data/reference/compartment_fire/pyrolysis/receding_slab_example_results.json`.
These authored heating cases do not establish experimental PMMA validation.

### Conservative charring slab shrinkage and swelling

`RecedingSlab::with_char` adds retained residue to the same finite-rate
virgin-to-char-and-volatile material, implicit heat solver and release accounts.
`CharProperties` supplies independent char density, conductivity, source and
a declared `ConductivityMixing::Parallel` or `Series` arrangement. The
existing `new` constructor remains the noncharring path.

The declared condensed volume is `V = m_virgin/rho_virgin + m_char/rho_char`.
The back remains fixed; front motion may be positive recession or negative
recession (swelling). This additive-volume choice is explicit, not the general
conditional FDS shrinking/swelling algorithm. Parallel conductivity is the
volume-fraction arithmetic mixture used by the pinned FDS Solid Chapter;
series conductivity is its harmonic counterpart for the separately declared
unresolved arrangement. All rebuilt internal links and the public
`cell_conductivity_w_m_k` and `front_half_cell_conductance_w_k` queries use the
same current component mixture.

Equal-retained-mass rezoning conserves total mass, char mass and physical
enthalpy, with virgin mass recovered as the complementary component. Char
fraction and sensible enthalpy use a paired limited reconstruction that
preserves donor means, nonnegative component inventories and the common
caloric-domain bounds. Formation references travel with their component
masses. The original scalar energy owner recovers mixed temperature with
reaction disabled during inversion; rezoning creates no volatile stream or
additional reaction heat. Separate virgin/char remap residuals prevent total
mass cancellation from hiding composition drift.

Frozen-geometry intervals are bounded by both the existing mass-loss limit
and `maximum_cell_volume_change_fraction` (default 0.25). This also controls
swelling with little mass loss. An excessive change retries the entire
unpublished interval, retaining no rejected release. Reports include maximum
volume change, component and total mass discrepancy, physical-energy
discrepancy, geometry roundoff and original caloric evaluation work.
Cancellation or failed budgets return no candidate.

Five additional native groups establish independent component inventory,
volume, fixed-back motion, outgoing enthalpy and both conductivity limits;
nonlinear-caloric remap and series-link consistency; better-than-threefold
remap refinement with independent phase energy-datum shifts; swelling retries;
and source refusal, cancellation and deterministic reuse. All 12 slab tests
pass, and the focused native thermal library passes 183 tests with one existing
performance test ignored. This is not a full-workspace DSR result.

Geometry/coefficient evolution remains split and first order in time, and
remap refinement does not certify global PDE accuracy or an entropy theorem.
No gas-volume displacement or pressure work is inferred. Pore transport,
permeability/porosity evolution, char oxidation, cracking, multi-step wood
chemistry, latent phases, complete burnout transitions and an experimentally
validated wood material remain outside this closure.


## Geometric conduction with evolving material conductivity

`enthalpy::conduction::ConductionObserver` prepares existing thermal links with
explicit centre-to-face distances, contact resistances and immutable scalar
`ConductivityLaw` bindings. It queries each participating cell once per actual
explicit substep and shares that value across incident faces. Internal links
reuse the original series resistance owner. A surface reservoir adds the
half-cell material resistance and an explicit areal film/contact resistance;
zero external resistance prescribes temperature at the face, not at the cell
centre. Zero conductivity explicitly insulates. Geometry, reference mass by default,
material identities and face bindings remain fixed and are checked on reuse.
Material selection and mixed-phase microstructure remain the constitutive
owner's responsibility.

The original observer protocol now accepts step-local reservoir conductances.
Both pairwise links and reservoirs enter the same original heat-rate and CFL
calculation. Their authored values are restored before accepted-step hooks and
events. No reservoir heat is reclassified as an applied source. Composed
observers reject conflicting reservoir ownership, and the existing advection
wrapper forwards the hook. `advance_observed_with_entropy` forwards each
physical hook once and retains the exact applied sources, links and reservoirs
for the original entropy account, including actual bath temperatures. General
internal transfers and applied sources retain explicitly unclassified entropy
contributions; temporal entropy defects remain distinct from physical
conduction production. Failure or cancellation publishes no partial physical
state.

The prepared graph uses O(cells + faces) storage and
O((cells + faces) log(cells)) work per substep with ordered maps. This is an
algorithmic bound, not a measured throughput claim. The existing explicit
coefficient-freezing and caloric CFL applies; it is not an implicit solver or
a time/spatial error certificate. Focused G0/G1/G3 tests cover series/Robin
resistance, dynamic boundary CFL and entropy, exact event splitting, independent
nonlinear transient and Stefan front/heat refinement, composition, insulation,
domain/work limits and transactional cancellation.

### Mass-dependent resistance of planar phase layers

`ConductionObserver::with_planar_phase_layers` explicitly binds selected cells
to constant bulk density and a fixed planar footprint. Current conserved mass
sets each normal half-thickness to `d = (m / rho) / (2 A)`. The original
Fourier/contact resistance owner uses this distance for both link sides,
reservoirs and ambient faces, while the original caloric and CFL owner uses
the same current mass. Growing frost therefore changes both storage and
resistance; its conductivity still comes from one original query per cell.
Unselected cells retain the fixed-mass refusal.

Admission checks the reference mass against density times reference volume,
the total two-face area against `2 A`, authored normal half-distances and
per-side contact measures. Source identity, reference mesh, contact coverage
and footprint remain fixed. A nonzero positive mass/thickness must remain
representable. All selected channels must be normal to those planar faces;
this is an explicit caller assertion, not an inferred geometric orientation.
The implicit chain/graph can reuse a selected layer at fixed mass; any later
mass change still refuses and requires rebinding the fixed-mass solve.

This closure supplies no mesh motion, sidewall growth, porosity, nucleation,
complete extinction, gas displacement work or density evolution. In particular,
mixed-phase calorics do not justify constant mixture density automatically.
The existing first-order thermal/phase split requires timestep refinement.
The four G0/G1/G3/G4 groups in `conduction/phase_layer_tests.rs` exercise growth
and retreat resistance, both link sides, finite heat storage/conservation,
temperature-dependent conductivity shared by both baths, per-side coverage,
unchanged default admission, cancellation and missing implicit bindings.

### Original material density in planar heat transport

`ConductionObserver::with_material_density_layers` selects the same planar
geometry with the original thermal law's optional `bulk_density_kg_m3(h)`.
The current conserved mass and density at the actual stored enthalpy set
`d = m / (2 rho(h) A)` for every normal link, reservoir and ambient face.
Density is evaluated once per selected cell during coefficient preparation;
the half-distance is reused across faces. Final-state density and distance
are admitted before an advance publishes, including after its last substep.
Unavailable density, unsupported enthalpy or nonpositive/unrepresentable
response refuses transactionally. The constant-density path is unchanged.

The original enthalpy law owns the pressure constraint and phase mixture.
Heating never reconstructs mass from density times the frozen reference
volume, and latent density changes remain visible when temperature is constant.
This mode requires enthalpy storage; an internal-energy law with moving volume
would need its own pressure-work equation. No extra p*dV heat is added.
The declared fixed-footprint, normal-expansion reduction supplies no lateral
strain, stress, mesh motion or gas displacement. Temporal and spatial
refinement remain necessary.

The implicit chain/graph admits these layers at fixed mass when the original
`ImplicitConductionLaw` additionally supplies `d(1/rho)/dh`. Each endpoint
queries the original density and tangent once per selected cell and reuses
its half-distance across normal links, reservoir and ambient faces. The exact
Newton derivative includes `d'/d = rho * d(1/rho)/dh` alongside `k'/k` in the
series-resistance derivative, including on isothermal latent plateaus. Native
state-dependent contact resistance composes with this same derivative. No
mass is reconstructed from the reference mesh and no second pressure-work
term is added. The same endpoint serves backward Euler, SDIRK2 and their graph
owners. Fixed-density layers need no extra tangent. Missing density tangents
still refuse.

Custom link, boundary and two-bath owners may implement their optional
`state_at_specific_enthalpies_with_geometry_tangent` or
`state_at_specific_enthalpy_with_geometry_tangent` query. It receives the
CURRENT normal half-distances and each local `(dd/dh)/d` in kg/J, and returns
the complete original flux tangent including that chain rule. The existing
fixed-geometry query remains unchanged; `None` explicitly refuses density
coupling even where a local slope happens to vanish. Unsupported prepared
source-integrated geometry is not silently substituted or differentiated.

`GrayAmbientBoundary`, `GrayTwoBathBoundary` and `GrayRadiativeLink` provide
this geometry-aware response through their SAME bounded surface solve. Their
half-cell resistance derivative uses `k'/k-d'/d`, including on latent plateaus,
and composes with original surface-coefficient and gap-resistance derivatives.
The geometric chain rule uses the existing source samples and adds no heat
account or pressure-work term.
Separate bath powers/temperatures and equal/opposite finite-body transfers
retain the native energy and entropy owners. Area, coverage, coatings, view
and independently declared gap geometry stay fixed; no gap displacement or
mechanical motion is inferred from a changing half-thickness.

G0/G1/G3 tests in `conduction/density_layer_tests.rs` cover a heating trajectory
with evolving resistance, both external faces, one density query per cell,
constant-density equivalence, missing data and final-state refusal. Material
integration in `fs-material/tests/density_conduction.rs` exercises the original
solid/liquid density law through a latent plateau and conserved heat transport.
`implicit/density_tests.rs` additionally checks analytic expanding/contracting
backward-Euler endpoints, the complete energy-coordinate link and two-bath
Jacobian, fixed-density reduction, closed chain/graph SDIRK conservation and
unsupported-law refusal. `fs-material/tests/phase_conductivity/density.rs`
drives a latent-density endpoint through the original source adapter and
checks its independently eliminated energy equation.
`implicit/density_tests/gray.rs` checks an analytic latent convection tangent,
nonlinear gray surface/gap finite differences, single original surface/gap queries,
both boundary channels against an implicit energy equation, distinct-bath
entropy, closed graph SDIRK conservation and cancellation.


## Implicit enthalpy transport on geometric conduction chains

`BackwardEulerStep::new_with_duration` accepts an independently owned physical
duration and its rounded absolute endpoint. They must agree within clock
roundoff; the endpoint still advances and original thermal event admission
still applies. Storage, source heat, transport and entropy quadrature all use
that same supplied duration. This lets integer native fluid clocks retain their
actual elapsed time at large absolute timestamps without replacing it by the
subtraction of rounded clocks. It adds no new heat-rate or integration owner.
A focused native regression defines three seconds of source heating where the
representable endpoints differ by four, plus inconsistent-clock, invalid-duration
and cancellation refusals. Execution is unverified on the current build host.

`enthalpy::backward_euler::BackwardEulerStep` additionally accepts optional
endpoint reservoir-conductance replacements through
`evaluate_with_reservoirs` and `evaluate_with_entropy_and_reservoirs`.
Original bath temperatures, areas, heat rates and entropy remain authoritative;
authored coefficients remain unchanged through evaluation and acceptance. Existing
coefficient callbacks retain their original behavior. Unknown cells, absent
baths, negative/nonfinite conductances and cancellation refuse transactionally.

`enthalpy::conduction::implicit::ImplicitConductionChain` solves fixed-mass,
nearest-neighbor chains using that same residual, not a second heat-rate owner.
Every caloric and conductivity law is bound to the original borrowed object.
A fused constitutive query supplies T(h), k(h), their analytic local derivatives
and an explicit admitted enthalpy interval. Physical latent plateaus retain
zero dT/dh while their conductivity may continue changing with phase inventory.
The exact endpoint Jacobian differentiates both the temperature difference and
the half-cell/contact/surface conductance, in cell ENERGY coordinates.

An initially insulating half-cell may have a source-declared one-sided
conductivity onset. Against a conducting neighbor or a finite boundary film,
`G=A*k/(d+k*R_other)` gives `dG/dh=A*(dk/dh)/d` at zero conductivity. The native
chain/graph Jacobian retains that finite derivative for links and both bath
channels while the current heat remains exactly zero. If the opposite half-cell
is also insulating, varying only one side cannot conduct. Equal-temperature
double onsets have zero complete-flux derivative and may evolve under independent
sources. Two active insulating onsets at unequal temperatures have a
direction-dependent harmonic-mean derivative; native Newton refuses that case
unless a complete link law supplies its own admitted tangent. Negative or
nonfinite conductivity/resistance and invalid source trials still refuse; no
coefficient clipping, artificial conductivity or extra heat storage is introduced.
G1 tests in `conduction/implicit/tests.rs` check analytic onset Jacobians, a
driven two-body cubic balance through both chain and graph, both bath channels,
and the exact quadratic two-body response from initially perfect contact.
Those three native regression definitions are syntax checked but unexecuted in
this environment because the complete dependency build is unavailable. Three
isolated tests of the unchanged production conductance/tangent functions passed
against independent resistance formulas and one-sided differences; this does
not establish a whole-solver or continuum-validation result.

Only three matrix diagonals are assembled. The reusable `fs-la` adjacent-pivot
workspace keeps the linear solve O(cells) in work and storage; original
map-based graph/constitutive bookkeeping remains. Newton candidates whose
nonzero increment cannot change the stored enthalpy are represented by exactly
zero increment before evaluation. This does not bypass or relax the original
underflow policy, heat rates, per-cell residual or total-energy acceptance.
Every resulting candidate is reevaluated by the original owner.

Declared controls bound timestep, accepted intervals, Newton iterations,
backtracks, interval halvings and total residual evaluations, including
rejected trials and final entropy queries. Bounded backtracking respects the
material domains; nonconvergence or a singular linear system may trigger a
whole-interval retry. Scheduled events split intervals at their exact times.
Only fully accepted intervals contribute heat or entropy. A complete call
also checks independently summed original stored energy against all supplied
heat. Failure may alter scratch but publishes no physical history. The report
retains original heat/entropy, maximum accepted local residual and work counts.
Every admitted caloric law must provide the original entropy integral.

This path removes the explicit diffusion CFL restriction for an accepted
implicit interval. Backward Euler remains first order; interval halving here
is nonlinear globalization, not a time-error estimator. Global uniqueness,
unconditional nonlinear convergence, general mesh motion and constitutive
accuracy beyond the supplied laws are not claimed. Fixed-footprint planar
density evolution is the explicitly selected geometric reduction above;
general network connectivity is supplied by the separate graph owner below.

Five focused chain regressions cover stiff exact solutions, latent traversal,
original-residual Jacobian differences, nonlinear time refinement, event
splitting, cancellation/budget rollback, singular falling-conductivity retries,
and a 96-cell diffusion tail with sub-ULP energy updates. Two added reservoir
regressions check original heat and actual-bath entropy plus refusal/rollback.
The full pinned-main native thermal source passes 191 tests, with one existing
ignored test; required workspace/DSR gates remain a separate unavailable check.

### Second-order implicit network convection and diffusion

`ImplicitConductionChain::advance_sdirk2_to` and the graph equivalent add
two-stage L-stable SDIRK2 with `gamma = 1 - 1/sqrt(2)`. Both stages solve the
original cell-joule residual using the original caloric laws, conductivity and
interface tangents, domain-respecting Newton search and direct/sparse linear
owner. The second stage retains the actual physical initial enthalpies; its
known first-stage contribution is the native integrated heat-rate vector,
not a synthetic initial temperature, an inferred storage difference, or an
external heat source. The accepted endpoint is the second stage.

First-stage heat is retained only for its final entropy-audited candidate,
with O(cells) extra storage and no extra per-Newton rate-vector allocation.
The ordinary backward-Euler entry points keep their original quadrature and
arithmetic. Both stages and every rejected interval share the existing work
budgets. Each fully accepted physical interval counts once; a rejected second
stage publishes neither stage. Independent local, interval-total and whole-call
joule gates remain active. Bath heat/entropy and conductive generation use the
two positive stage weights; exact caloric entropy storage is measured from the
original initial state and the temporal entropy defect remains explicit.

This entry point admits fixed-mass chains/graphs, including their existing
nonlinear bulk, contact and integrated interface laws, with piecewise-constant
authored cell source/sink powers, ambient/reservoir temperatures and native
Robin htc on unbound ambient channels. Sources enter both native stages with
the same positive weights as conduction; each joule is counted once in the
external energy account. Their entropy remains explicitly
unclassified at-cell heat, not an inferred reservoir entropy or dissipation law.
Admitted `Event::Boundary` changes split each interval at the exact event
time. Both stages use the left-limit boundary values before processing the event
through the original timeline/observer. Simultaneous changes retain authored
order; target-time events change the returned boundary without transferring
heat under its new values. Each interval's heat and reservoir entropy use its
actual bath temperatures; a bath jump adds no instantaneous cell heat or
entropy. No storage or cost of changing an external reservoir is inferred.
Prepared material, geometry and boundary bindings remain unchanged.
One whole-window preflight, with bounded cancellation polling, rejects changes
to reservoir presence/area or authored conductance, bound ambient htc, unknown
cells and all topology/birth events before stage work. Geometric and complete
boundary laws own their conductance; accepting changes to overridden authored
coefficients would silently ignore the requested physics. Native Robin htc
may change, including to zero for insulation. Complete boundary laws receive
the current bath temperature, including the gray exterior convection/radiation
law; their emissivity, film coefficients and surface resistances remain fixed.
The same whole-call step and nonlinear work budgets cover every event interval;
events never renew them. The nonadaptive constant-boundary path avoids event
processing/cloning at steps where no event is due. External source observers
and moving mass remain outside this entry point. Complete-interval halving in
`advance_sdirk2_to` handles nonlinear failure and the passive entropy admission
below; it does not estimate temporal accuracy. Second-order time accuracy
requires a smooth admitted solution and independent timestep refinement;
constitutive/phase/limiter corners may lower it. L-stability is not unconditional
positivity, monotonicity or a maximum-principle guarantee: every stage retains
the original material-domain checks, without temperature/heat clipping.

For passive thermal networks, a converged positive-temperature endpoint must
also retain nonnegative **cell-plus-bath** entropy change within the existing
algebraic/roundoff allowance. The check subtracts the actual stage-weighted bath
entropy input from the original exact caloric storage change; ordinary cooling
may decrease cell entropy while increasing total entropy. Nonnegative-conductance
links, contacts and single-bath boundary laws retain their original ownership.
Any attached advection, heat-source law, unclassified source power or generalized
two-bath exchange callback excludes this check: those callbacks do not assert
passivity, and no missing source/work/stream entropy is inferred. Passive
radiative laws keep their existing bulk/bath convention, without a claim about
photon entropy or an additional radiative field closure.

The allowance converts the two stages' absolute local joule acceptance envelope
`n * cell_tolerance_j / gamma` through the minimum actual initial/stage/bath
temperature. A compensated absolute-exchanged-heat reduction supplies a small
floating-point allowance, counting both entropy endpoints of each transfer.
Neither an absolute enthalpy nor an entropy reference datum sets this scale.
Only a negative candidate needs these extra coefficient/coverage reductions;
already evaluated stage laws are reused without another caloric integral or
nonlinear constitutive solve. `EntropyDecrease` retries the original complete
interval through both fixed and adaptive controllers, and exhausted retries
publish neither a state nor its accounts. Stored heat, temperatures, reported
entropy and the signed temporal defect are never repaired. This finite-step
admission is separate from temporal accuracy and maximum-principle guarantees.

An attached `with_implicit_advection` binding now participates in both stages,
including `with_implicit_advection_reconstruction` on the graph. Its prescribed
mass rates, inlet enthalpies, original caloric references and fixed inventory
remain unchanged. Every Newton trial queries the original transport owner at
that stage's actual enthalpy; a limited face and its generalized tangent are
never frozen at the initial state or first stage. The native stage heat vector
retains equal opposite internal-transfer rates, and the external energy account
integrates inlet-minus-outlet enthalpy with the same two positive weights as
authored power and bath heat. Internal circulation is not external heat.
Signed transport Q/T stays unclassified at-cell entropy with those same weights
unless the original flow entropy opt-in is selected. Bulk/liquid carriage and
source-anchored stationary mixtures, including bound selective streams, use
their own supported entropy closures. Reconstructed faces and attached explicit
flow histories retain their separate endpoint restrictions.
No explicit flow CFL is imposed; finite-volume spatial error still requires
independent refinement and the supplied flow observer's history is not advanced.

Focused regressions compare unequal-capacity equilibration and bath cooling
against independent analytic stage solutions, including physical weighted
entropy. Periodic graph Fourier decay and nonlinear `k(T)` cooling demonstrate
second-order temporal refinement. Stiff modal damping, latent-plateau storage,
second-stage domain/work rollback and unsupported event refusal exercise the
stated admission boundary. Constant heating/cooling regressions cover unequal
capacities with independent two-stage temperatures and entropy, driven bath
relaxation with a three-rung second-order ladder and adaptive accuracy, and
latent-plateau storage with a large caloric reference datum.
Scheduled-source regressions independently integrate heater/sink joules across
non-grid-aligned and simultaneous events, target-time changes and continuation
on all chain/graph and fixed/adaptive entry points. A piecewise exponential
bath solution retains second-order temporal refinement. Late unsupported
events, whole-call step budgets and cancellation during a long source-schedule
preflight refuse without publishing a prefix.
Boundary-schedule regressions independently compose exponential bath solutions
on all fixed/adaptive chain/graph entry points, retain second-order refinement,
and reconstruct stage heat and entropy across bath changes. Native Robin
switching covers an insulated interval, exact target-time changes and
continuation. A gray exterior law follows the same event-aligned trajectory
and accounts as manually changing its bath between separately prepared solves.
Changes to ignored coefficients refuse before any stage work.
Convection regressions independently check both-stage paired-transfer entropy,
directed Fourier convection-diffusion phase/decay with second-order time
refinement, and open-stream/source/bath joules and entropy under a large enthalpy
datum shift. A limited latent channel has independent scalar stage solutions
and an exact semi-discrete transient: its enthalpy converges at second order
while all temperatures remain fixed. Tightening adaptive enthalpy tolerance
improves that solution. Source-event timing, late cancellation and unsupported
boundary-geometry event preflight are also exercised with attached flow.
These numerical checks do not assert spatial second order at limiter switches, unconditional
TVD/positivity, a stream-entropy closure or experimental validation.

#### Adaptive second-order network convection and diffusion

`advance_sdirk2_adaptive_to` on the same chain and graph owners adds local
h/T/phase error control through the existing `AdaptiveConductionStepper`.
Each comparison uses one full and two half native SDIRK2 intervals, sharing
all six stages' residual/linear work budgets. Only the conservative fine pair
contributes accepted heat, entropy and physical step counts; no Richardson
state or rejected/coarse history is published. The original independent local,
interval-total and whole-call energy gates remain active.

`ConductionTimeOrder` defaults to first order, preserving the backward-Euler
factor and IEEE square-root controller for existing adapters. SDIRK2 selects
second order: actual representable half spans `a,b` use fine-error factor
`(a³+b³)/(3ab(a+b))`, equal to `1/3` for equal halves, and cubic local-error
scaling with the existing safety/growth/reduction bounds. The factor uses a
dimensionless ratio to avoid cubing seconds. Bounded deterministic bisection
within the existing controller limits avoids platform cube-root/transcendental
arithmetic and overflow from multiplying a large normalized error.

The complete requested window is checked before the generic controller can
process an event or compute a comparison. The generic controller already clips
both coarse and fine trajectories at admitted boundary events and processes
them only after accepting the preceding interval. Attached prescribed implicit
advection and optional limited reconstruction participate in every stage of both
trajectories. All other events through the target remain outside SDIRK2's
admitted scope even if an initial comparison would end sooner. The adapter shares
the original material observables and change-based, caloric-datum-invariant
error scales. It adds no phase proxy, source splitting or global accuracy bound;
constitutive corners and spatial error still require independent refinement.

Focused native regressions cover accuracy/tolerance refinement, latent inventory
at constant temperature, caloric datum shifts, sparse graph conservation,
unequal representable half spans, shared six-stage work ceilings, late rollback
and complete-window refusal. The production order arithmetic also has directly
executable standalone tests for first-order compatibility, second-order scaling,
extreme finite spans and invalid controller inputs.

## Additional material conductances on the native explicit timeline

`DepositionSimulation::limit_timestep_for_additional_conductances` composes
nonnegative extra W/K with an already admitted explicit step using the original
cell masses and registered global caloric dT/de bounds. It adds the largest
extra decay rate harmonically to the current rate allowance, rather than taking
two independent CFL minima. An empty or zero map leaves the proposed step
bit-for-bit unchanged; zero caloric slope adds no restriction. Unknown cells,
invalid rates, unrepresentable steps and checkpoint failures refuse.

This helper applies no heat and changes no thermal state or authored boundary.
Call it before event clipping through the existing observer timestep hook.
It can be conservative when different cells or a caller maximum step dominate.
It bounds frozen heat-transfer coefficients, not their nonlinear derivatives
or time-discretization error. Focused native regressions cover combined
conductive/cooling limits, sequential composition, zero slopes and refusal.
The airflow material observer uses the exact hA from its start-state source
evaluation, so changing fluid data changes both cooling and native subcycling.

## Shared endpoint ambient and reservoir conductances

`backward_euler::{BoundaryConductances, BoundaryStepCoefficients}` pass endpoint
total ambient/reservoir G directly into the original heat-rate and entropy
owners. The boundary-aware residual and external-power APIs retain prescribed
temperatures, reference coverage, authored h/G and cell watts. Replacement G
is applied once; no boundary input is mutated or reclassified as generated heat.
The same G participates in signed physical watts, outgoing-conductance/CFL
accounting and boundary entropy Q/T. Existing ordinary callbacks use empty maps;
the reservoir-only APIs delegate their maps through this same general ingress.
Finite thermal stores remain original internal links, not prescribed baths.

Unknown cells, absent reservoirs, negative/nonfinite G and nonzero Robin G
without exposed reference coverage refuse. Three added native regressions cover
independent finite-cell energy/entropy/CFL, unscaled applied watts and empty-map
parity, invalid coverage and final cancellation/replay. These tests are compile
checked in this change; their Cargo runtime build was killed while compiling
the existing asupersync dependency. The earlier native test counts above are
prior validation, not runtime evidence for this extension. Geometry-model
validity and nonlinear/time accuracy remain the composing caller's responsibility.

## Declared thermal storage and closed-volume admission

`ThermalLaw::storage_basis` declares whether the original caloric law stores
specific enthalpy, fixed-volume specific internal energy, or declared sensible
heat storage whose capacity does not identify either potential. Existing
enthalpy laws retain their default; generic measured or assembled capacities
use `DeclaredSensibleStorage`. The inherited `specific_h_j_kg` field carries the selected
coordinate; the timeline applies the same conservative heat exchange and
`integral de/T` entropy account without converting one thermodynamic potential
into another. Each law can admit its own mass/geometry relationship through
`validate_cell`, called before a birth or phase exchange changes inventory.

`material_storage_basis` exposes the original registered law's declaration.
Nonzero prescribed mass flow and inlet calorics require enthalpy storage. The
single-surface, joint-chamber and external-endpoint sublimation paths also require
enthalpy before constructing or publishing open-mass exchange. Batch admission
finishes before any endpoint mutation. Zero-flow connections remain available
for heat-only internal-energy and declared sensible storage. A reversible
`integral de/T` opt-in does not grant either source an enthalpy interpretation.

Direct liquid evaporation/condensation, direct and shared-chamber sublimation,
and externally solved solid/condensed endpoints validate the complete candidate
cell with its original caloric law before committing mass or enthalpy. This is
the same admission used by external liquid exchange. A caloric law bound to a
fixed component inventory or reference mass/volume cannot silently acquire a
different heat capacity through mass loss. All cells and the shared vapor
reservoir remain unchanged if any candidate is refused; zero-transfer endpoints
remain admissible when they preserve the original binding.

Focused regressions exercise two rigid bodies exchanging heat with conserved
internal energy and increasing thermal entropy; invalid mass/volume births;
forward/reverse circulation and boundary-flow refusal; and all three sublimation
entry paths. Mixed enthalpy/internal-energy batches preserve the original solids,
vapor, mass, stored energy and clock on refusal. This adds no open-mass internal-
energy model, pressure-work conversion or equation of state. Such conversions
remain the explicit responsibility of the physical material/equipment owner.

## Adaptive implicit thermal accuracy

`enthalpy::conduction::implicit::adaptive` adds first-order step-doubling
control around the original `ImplicitConductionChain`. One full backward-Euler
interval is compared with two actual half intervals; only the conservative
half-step trajectory is accepted. Comparison trials never contribute physical
heat, entropy or accepted time, and no Richardson-extrapolated state is formed.
Unequal representable half spans receive their first-order error correction.

`enthalpy::conduction::implicit::adaptive::AdaptiveConductionStepper` supplies
a typed endpoint seam for the existing adaptive controller. The original
conduction chain and material-integrated layered endpoints share the same
step-doubling comparisons, change-based storage/temperature/phase error scales,
exact event grid, original energy/entropy accounts and whole-call work ceilings.
An adapter must retain original cell/material identities, charge every residual
evaluation and completed nonlinear iteration including refused trials, and
return accepted original backward-Euler endpoints without internal physical
subdivision. This interface adds no conductivity, caloric law, extrapolated
state, global-error certificate or trajectory adjoint.

`ImplicitConductionLaw::time_observables` supplies original temperature and an
optional original phase mass fraction. The per-cell estimator checks stored
specific energy, temperature and, when requested, phase inventory. Relative
scales use changes from the accepted state, rather than absolute enthalpy or
Kelvin values. Phase accuracy refuses when an original observable is unavailable
rather than inferring phase from temperature. The h-labelled coordinate retains
the declared thermal storage basis.

Controls retain the original local/whole-call joule gates and bound comparison
spans, comparison count, accepted half steps, consecutive reductions, Newton
work and the shared residual-evaluation budget. Scheduled events are hit exactly.
Source refusal, exhausted work, cancellation and unresolved accuracy return no
partial trajectory. Reports separate local time-error estimates from algebraic
energy residuals, count all attempted nonlinear work, and expose the next proposed
comparison span for continuation.

Four native regressions cover analytic cooling and caloric-datum invariance;
phase-driven refinement on an exactly isothermal latent plateau; non-dyadic bath
events and fine-only heat/entropy; and atomic refusal for missing observables,
minimum-step limits, trial/residual budgets and late cancellation. A local
step-doubling indicator is not a certified global time error or spatial error
bound. The existing fixed-chain geometry and material-domain restrictions remain.

## Implicit conduction on branched thermal graphs

`enthalpy::conduction::implicit::ImplicitConductionGraph` extends the same
fixed-mass enthalpy solve to branched, cyclic and disconnected graphs, including
parallel interfaces and asymmetric covered contact areas. It delegates to the
original backward-Euler residual, nonlinear line search, reservoir/event owner,
local and global energy gates, entropy audit and adaptive step-doubling logic.
Cells, caloric laws, conductivity laws and geometric bindings remain original.

The graph Jacobian retains the analytic derivatives of both T(h) and k(h),
including zero dT/dh on latent plateaus and energy-coordinate mass factors.
It stores one diagonal and two off-diagonal contributions per physical link.
The existing `fs-solver` FGMRES owner applies this generally nonsymmetric
operator without assembling a dense global matrix. Each current Jacobian is
checked for acyclic off-diagonal support after combining parallel contributions.
Forests use exact leaf elimination, retaining unequal directional derivatives,
one-sided couplings and disconnected cells. Forest assembly and factorization take
O(cells + Jacobian contributions), and each forest inverse action takes
O(cells) without scratch allocation. Cyclic graphs use the existing `fs-sparse`
zero-fill ILU owner on the same nonsymmetric Jacobian, with checkpointed setup
and substitution. Canonical duplicate assembly takes O(entries log entries);
factor storage and each inverse action remain O(cells + entries), while setup
work depends on sparse row intersections. Neither path shifts the physical
matrix or assumes positive definiteness. Unsafe factors or inverse arithmetic
fall back to diagonal scaling. The Krylov solver still
verifies the true residual of the original, uncombined contributions.
Each operator action is O(cells + links), and bounded restart storage is
O(cells * restart + links). Original map-based physical bookkeeping remains.
The direct `ImplicitConductionChain` path retains its O(cells) tridiagonal solve.

`ImplicitConductionLinearControls` bounds restart dimension, restart cycles
per nonlinear direction and total Krylov iterations for the whole advance.
Coarse, fine, rejected nonlinear and rejected adaptive trials share that
iteration ceiling; `linear_iterations` reports completed work. Linear
nonconvergence may use the same interval-halving path as nonlinear refusal.
Cancellation and exhausted budgets publish no physical candidate. Linear
residual convergence never replaces the original joule acceptance gates.

Regression coverage includes stiff star equilibration with unequal capacities,
independent endpoint elimination, a manufactured 49-cell branched endpoint with
parallel contacts and a one-iteration linear budget, a stiff cyclic loop under
an eight-iteration budget, sparse-factor fallback, energy and entropy accounts,
latent branches,
active cycles and parallel links, an isolated cell, asymmetric covered areas,
the original-residual Jacobian and transpose, adaptive analytical diffusion,
shared work budgets and cancellation. Geometry and reference masses remain
fixed; backward Euler remains first order. No unconditional nonlinear
convergence, spatial accuracy, hardware throughput or experimental validation
beyond the supplied material laws is claimed.

## Evaporation on the original thermal timeline

The shared `evaporation::EvaporationLaw` selects an explicit immutable
`InterfaceKinetics`. Its default `Equilibrium` retains the original Stefan-film
arithmetic. `HertzKnudsen` and `Schrage` require a caller-supplied constant
accommodation probability in `[0,1]`; zero seals the interface. The borrowed
`KineticEvaporation` adapter delegates the complete original evaporation and
thermal source, including caloric identity, storage basis, density, cell
admission, entropy and species-advection bindings. No material property or
accommodation probability is inferred from a name.

Finite kinetics solve the nonlinear series balance
`j = K*(p_sat-p_i) = C*ln[(p-p_v)/(p-p_i)]`, with
`C = h_m*M*p/(R*T_f)` and `K = a*sqrt(M/(2*pi*R*T_s))`.
Hertz--Knudsen uses `a=alpha`; Schrage uses `a=2*alpha/(2-alpha)`.
The local gas and liquid share `T_s`; the representative film temperature
remains the caller's `T_f`. The bounded, allocation-free scalar solve retains
pressure differences directly and solves the smaller resistance's pressure
drop, preserving evaporation, condensation and both limiting regimes. The
analytic surface-temperature tangent includes `dK/dT_s=-K/(2*T_s)`.
At zero net pressure drive it retains the combined conductance; a physically
sealed film or interface has both zero rate and zero tangent.

Scalar implicit latent steps, all original timeline variants, and the shared
implicit network's initial and endpoint rates use this same selected closure.
Mass, recession and complete vapor-stream enthalpy consequently derive from
one transfer without a second latent source. G0/G1 regressions use manufactured
interfacial pressures to check both nonlinear equations, finite-difference
tangents, film/kinetic limits, pressure gaps below the total-pressure ULP,
bitwise default equivalence, source delegation, scalar energy closure and
native heat-network/inventory-retry integration. These establish implementation
behavior when executed, not measured accommodation data or experimental accuracy.
The equations' stated small-temperature-jump scope follows Lu et al., Nature
Communications 10, 2368 (2019), equations (7)--(8),
<https://doi.org/10.1038/s41467-019-10209-w>.

The finite interface does not extend the gas film to pure vapor, boiling,
rarefied transport or a resolved Knudsen layer. It supplies no nonisothermal
molecular heat accommodation, curvature correction, nucleation or spray model.
The original material, gas-film, geometry and pure-liquid validity checks remain.

`ParallelFilmEvaporation` borrows one original liquid law and an immutable
slice of `EvaporationPatch { area_m2, mass_transfer_coefficient_m_s }`. Its one
inventory uses the returned total area. All patches share the actual liquid
temperature, gas pressures, representative film temperature and selected
interface kinetics; each retains its own molecular coefficient. The result is
the compensated area-mean of independently closed nonlinear patch fluxes and
temperature tangents. Coefficients are never averaged before finite kinetics.
Latent heat and recession derive from the same represented mean mass and the
original common latent heat/density. Zero-area patches contribute no evaluation;
zero coefficients seal individual patches. Finite nonnegative inputs and
positive represented total area are required. Uniform coefficients preserve
the original single-film arithmetic and require only one flux evaluation.

The object-safe `EvaporationLaw::evaluate_film` hook receives already-admitted
source properties, so adding patches does not multiply original material
lookups. Public flux queries, scalar implicit heat/phase steps and both initial
and endpoint network rates dispatch through it. `KineticEvaporation` delegates
this hook with the selected kinetics, preserving patch physics in either wrapper
order. Both adapters retain all original thermal capabilities. Construction and
evaluation allocate no memory; heterogeneous evaluation is bounded by the patch
slice. Callers retain geometry/connectivity, distinct-face admission, work bounds
and the one original mass/energy account. No local gas field or second inventory
is introduced. Added regressions cover independent per-patch Stefan/kinetic
equations in both directions, changing surface temperature, tangent and wrapper
composition, one material query per state, uniform fast-path equality, original
thermal delegation, scalar/network heat balance and cached inventory retries.

`evaporation::droplet` adds a spherical pure-liquid inventory whose only retained
coordinates are mass and original specific enthalpy. Radius, area and volume
are derived from mass and the original liquid density at the actual temperature;
accepted geometry is rebuilt from both final mass and final density. The
stagnant continuum sphere uses the diameter-based `Sh=2` result, `k_m=D_AB/r`,
with explicitly supplied binary molecular diffusivity. The original nonlinear
Stefan film, selected finite interface kinetics, heat balance, latent energy,
inventory refusals and cancellation remain the sole phase/thermal owners.
Geometry construction and the adapter allocate no memory.

`advance_droplet_with_heat` freezes initial area, radius and molecular diffusivity
within one interval and refreshes them after acceptance. This geometric split
is first order; refine the whole timestep. The receipt retains both spherical
geometries and the original phase/heat account. Its planar recession coordinate
is not the spherical radius change. Constant-property, constant-driving-force
steps converge toward the classical diameter-squared law; finite humidity,
thermal cooling, density variation and finite kinetics need not follow that law.
Added regressions cover independent sphere geometry and `D/r`, timestep
convergence with the original caloric inverse, heat-driven density/size changes
at zero mass transfer, finite-kinetic sealing, input and inventory refusals, and
cancellation. These tests were added but could not be executed in this session:
the pinned build was killed in the unchanged upstream dependency `asupersync`.
The isolated, uniform-temperature sphere assumes quasi-steady continuum radial
transport and negligible forced/buoyant convection; it adds no moving spray,
Kelvin correction, rarefaction, nucleation or complete-dryout event.

`enthalpy::evaporation` connects the continuum liquid/gas-film owner to the
actual `DepositionSimulation` cell mass and original caloric enthalpy.
`evaporate_cell` changes these two coordinates atomically and leaves the clock,
identities, reference geometry, contacts and authored sources intact.
`advance_evaporating_to` composes original heat transport followed by implicit
phase transfer on event-clipped intervals. The interval's initial exposed area
survives endpoint covering, and no phase exchange occurs before cell birth.
The whole request shares the original thermal work limit.

The split timeline uses the same birth-relative, compensated open-mass account
as the coupled and implicit evaporation owners. Its final gate reconstructs
each cell's enthalpy change and the vapor stream relative to that cell's
original or birth datum before reducing; an unchanged high-energy neighbor
cannot hide or reject a smaller latent transfer. Actual mass changes have a
separate balance. Full vapor enthalpy and storage changes remain reported on
their original references, and latent heat is not counted again. The split
regression compares both evaporation and condensation, including delayed
birth, against the identical isolated specimen.

`evaporation::advance_surface_with_heat` and the native
`DepositionSimulation::evaporate_cell_with_heat` solve signed supplied heat
and phase transfer at one endpoint. The original liquid chart closes
`m0*(h1-h0)+dm*L(T1)=Q`; the vapor still carries `dm*(h1+L(T1))`, so latent
heat is counted once. Supplied heat can feed evaporation directly without an
artificial sensible-only intermediate temperature outside the liquid chart.
The returned surface receipt retains `supplied_heat_j` and the complete
energy residual. A zero-duration heat impulse changes caloric storage with
zero phase transfer, subject to the same liquid and nonboiling admission.

The caller supplies total signed energy and owns its reservoir account.
Temperature-dependent convection, radiation or contact exchange still needs
its own coupled owner. The direct cell method preserves the original clock,
geometry and authored sources; it does not automatically apply their heat.
The existing `advance_evaporating_to` timeline remains first-order heat/phase
splitting. Its original calls use zero additional supplied heat.

`advance_coupled_evaporating_to` feeds the original start-state source,
convection, reservoir and contact heat directly into that phase endpoint on
each native interval. Other cells receive their original explicit heat update,
including the equal-and-opposite contact exchange. The liquid never visits a
sensible-only intermediate state. Original conductive CFL bounds are computed
again after every actual mass change; one substep budget spans the whole
request. Birth, contact and boundary events clip both heat and mass exchange,
with exposure changes applied only after the preceding interval completes.

`CoupledEvaporationAudit` reports actual mass/storage changes, incoming birth
mass/energy, external heat and full carried vapor energy. Latent energy is a
subaccount of that vapor stream. Compensated balances cancel each cell's
original or birth enthalpy datum before summing, avoiding subtraction of large
absolute energies. Failure or final cancellation
publishes no state. This method is first-order IMEX: the phase endpoint is
implicit but contact and boundary heat remain explicit and require temporal
refinement. It does not remove the original thermal CFL restriction or evolve
the prescribed gas film or reference geometry.

Six native coupled-timeline regressions cover direct latent heat supply,
independently solved source/contact/bath balances, mass-dependent CFL limits,
birth and coverage events, large shifted birth-energy datums with an unchanged
high-energy neighbor, and late event/cancellation rollback.

`evaporation::implicit::EvaporatingBackwardEulerStep` assembles original
backward-Euler heat transport and signed phase transfer at one common endpoint.
Its unknowns retain `m0*(h1-h0)` in the original stable cell order. The liquid
equation includes `dm*L(T1)` once; acceptance installs the actual remaining
mass `m1=m0-dm` and reports complete vapor enthalpy `dm*(h1+L(T1))`. Contact,
source, convection, reservoir and resolved bath/radiation inputs use the same
original backward-Euler coefficient ingress. No explicit conductive CFL is
imposed by this residual adapter.

The resolver receives both the original thermal endpoint with `m0` retained
for its energy coordinates and an explicit proposed liquid state with `m1`.
Mass-dependent sources must use that proposed inventory. One original thermal
rate pass measures physical heat before latent exchange; the native residual
pass includes the latent sink. This preserves an independent physical account
without recovering small heat inputs by subtracting large mixed source terms.
The source law is queried once per endpoint evaluation. Original authored
boundaries and topology survive acceptance; events at the endpoint remain
pending until coupled owners accept the preceding interval.

`EvaporatingBackwardEulerAudit` separates physical external heat, actual
changed-mass storage and the complete surface receipt. Compensated relative
energy and mass accounts accompany the original local and whole-system
backward-Euler gates. Represented mass-rate error, remaining inventory,
original material-cell constraints and final cancellation are checked before
returning an accepted candidate. The caller supplies the nonlinear solver and
bounds its work; this API does not itself solve the nonlinear system, estimate
time error, or evolve finite gas, geometry, boiling or dryout.

Seven native implicit regression groups cover an independent stiff two-cell
solution, temperature-dependent latent heat with an enthalpy datum shift,
proposed-mass sources and resolved bath/radiation/pair-transfer inputs, both
flux reversals, original inventory/source-domain refusal, independent
local/mass/cancellation gates, and pending endpoint event exposure.

`DepositionSimulation::solve_implicit_evaporating_step` drives that endpoint
with the original `fs_solver::NewtonKrylovState` and fallible numerical
Jacobian actions. It solves the entire active heat network with one selected
liquid interface; other cells retain their original heat-only response.
`EvaporationImplicitControls::new(range)` requires an explicit liquid branch
and keeps nonlinear tolerances, physical joule/kilogram gates and bounded work
separate. No conductive CFL restriction is introduced.

A two-entry cache avoids duplicate residual reconstructions in numerical
Jacobian actions. The final physical state and coefficient payload are always
recomputed, charged to the same evaluation budget and passed through the
original local/global energy and mass gates. Numerical convergence cannot
replace physical acceptance. Source/domain/stencil errors retain their cause;
failed work or final cancellation publishes no candidate. The method uses
numerical derivatives, requires a resolvable difference scale, and remains
first order in time without automatic time-error estimation.

`DepositionSimulation::advance_implicit_evaporating_to` now drives those solves
on the original event timeline. Start-time births/contact/boundary changes are
processed in authored order; the chosen liquid must be active before a
positive-duration solve. Later births of other cells, covering and separation
retain their original semantics. Each interval ends at the next maximum-step
endpoint or event; endpoint events change the following interval's exposure.
The last returned coefficient payload describes the accepted interval before
those endpoint events, and no payload is produced for an event-only request.

`EvaporationImplicitTimelineControls` separates maximum physical timestep from
the numerical solve controls. Newton attempts, uncached physical residual
evaluations and original thermal substeps are bounded across the entire
request. Absolute interval energy and represented-mass-rate defects accumulate
without sign cancellation under unchanged whole-request tolerances. The
implicit and IMEX timelines share the same compensated physical mass/energy
account, including birth enthalpy and complete carried vapor energy. A late
source/event failure or cancellation publishes neither intermediate state nor
an interval payload.

Five native timeline regressions compare manual interval composition and
independent source/exposure balances, resolve nonlinear time/mass-dependent
heating, enforce shared work limits and final rollback, and refine toward the
analytic variable-mass heat equation. Backward Euler remains first order; the
maximum timestep still controls physical time error independently of algebraic
convergence. That maximum-step entry point does not infer automatic time-error
control or evolving gas/geometry.

`EvaporatingNetworkBackwardEulerStep` and
`DepositionSimulation::solve_implicit_evaporating_network_step` bind multiple
liquid interfaces to that same original heat graph. Each `EvaporatingCellBinding`
owns its cell id, original evaporation law, prescribed film and admitted liquid
enthalpy range. Heterogeneous laws and simultaneous evaporation/condensation
retain separate local receipts and acceptance gates; signed transfers cannot
hide another surface's energy or represented-mass-rate defect. Unselected cells
keep their original caloric response, including latent storage, and an empty
binding slice explicitly selects heat-only transport.

Every endpoint source callback sees all proposed liquid inventories in canonical
cell order. One interface query per surface and two original thermal-rate passes
resolve the whole network; there is no per-surface heat solve or sequential phase
split. Exposed areas are prepared by one traversal of the original contact graph,
preserving asymmetric coverage and reservoir areas. All masses are installed
atomically after native local/total heat gates, independent open-system balances
and every original material validation pass. The original single-cell residual
and automatic solve delegate to these shared owners, preserving their explicit
liquid range, callback and unsized-law support. Geometry and films remain
prescribed, and endpoint convergence alone does not establish time accuracy.

`DepositionSimulation::advance_adaptive_evaporating_network_to` adds local
time-error control to that shared multi-surface solve. Each comparison solves
one backward-Euler interval and two actual half intervals. The original
`ConductionTimeOrder::First` supplies the correction for unequal representable
halves and bounded step-size factors. Only the conservative fine trajectory is
accepted; mass and enthalpy are never extrapolated.

Time accuracy measures original specific stored energy and temperature for every
active heat cell and actual mass for each selected liquid independently. It can
therefore detect latent-storage error at constant temperature and mass-loss error
when two endpoints have the same temperature. Relative allowances scale with
changes from the accepted state, preserving the original caloric datum and
avoiding an allowance set by a large retained inventory. Absolute time-error
allowances in J/kg, K and kg remain distinct from nonlinear joule/kilogram gates.

Comparisons stop at original events, which retain authored order and exposure
semantics. Only accepted fine intervals contribute external heat, latent energy,
complete carried vapor enthalpy or the returned coefficient payload. Compensated
`surface_totals` preserve each cell's signed mass/latent/vapor/recession account,
including simultaneous condensation and evaporation that cancel globally.
`last_surfaces` keeps the final accepted half-interval receipts; endpoint events
still follow that interval's coefficient payload.

All coarse work and accuracy-rejected comparisons spend the same whole-request
Newton and residual-evaluation budgets. Comparison count, accepted half-interval
count and consecutive accuracy reductions are bounded separately. Absolute
accepted energy and represented-mass-rate defects accumulate without cancellation,
and original whole-trajectory mass/energy gates still govern publication. A
source, nonlinear, material, work-budget or cancellation failure returns no
partial state. A final or event-clipped span below the step floor must still pass
accuracy. Reports expose accepted error estimates, rejected work and the next
suggested comparison span. Backward Euler remains first order; the local
estimator is not a certified accumulated-error bound, and source corners require
independent refinement. Reference geometry and gas films remain prescribed.

Adaptive liquid bindings may select cells that are active now or have an original
scheduled birth. Selection is resolved once before source callbacks; unknown
unscheduled identities refuse. Heat transport continues through existing cells
before a selected liquid is born. Each selected surface activates only after its
native birth has been admitted, with its original mass, enthalpy and current
contact exposure. It receives no earlier phase transfer, mass-error entry or
surface total. A birth exactly at the requested endpoint remains absent from the
preceding interval's coefficient payload and phase receipts, while a later birth
stays pending. Actual birth-time material/branch failures and late cancellation
still discard the whole proposed trajectory.

`EvaporatingCellSource` and
`solve_implicit_evaporating_network_step_with_films` admit original liquid laws
with source-resolved pressure, vapor partial pressure, gas temperature and mass
transfer coefficients. Initial films are queried once from the actual initial
state. `EvaporatingNetworkBackwardEulerStep::evaluate_with_films` then resolves
the complete film map at each uncached endpoint BEFORE evaluating any phase rate.
The heat callback receives those exact validated films and all proposed liquid
inventories, so heat and environmental response participate in one residual.
The thermal view contains endpoint time and caloric coordinates with the
interval's initial masses; proposed new masses are supplied in the liquid map.
Every film map must match selected active cell identities exactly; missing,
extra, wrong-id or invalid films refuse before the heat callback.

Fixed and source-resolved paths share the same physical residual, Newton loop,
two-entry cache, numerical work bounds and independent final acceptance. The
final endpoint queries both sources afresh; only its payload is returned. The
residual budget counts complete endpoint evaluations, with one additional
bounded initial-film query for preparation. Callbacks retain immutable accepted
history during all trials. A film temperature supplies the mass-transfer model;
reciprocal sensible heat remains explicitly supplied thermal coefficients.
This entry supplies environmental feedback, without creating an atmosphere
inventory or a gas caloric model. A finite-gas coupling must separately close
species/energy changes and upstream-vapor enthalpy for nonisothermal condensation.

`advance_adaptive_evaporating_network_with_films_to` advances those sources on
the same adaptive event/phase timeline. Both public adaptive entry points share
one controller, activation path, physical account and whole-request work budget;
only their per-interval source evaluation differs. Future liquid births need no
invented film state. Resolved maps contain only active selected cells, including
an empty map during a heat-only prefix.

Each coarse or fine interval resolves its own initial films and nonlinear
endpoint films from its own thermal state. In particular, the second fine
interval uses the first fine interval's changed inventory and enthalpy. Only
accepted fine transfers enter per-cell totals, and only the final fine payload
is returned before endpoint events. For a successful request the additional
initial-film queries are exactly three per complete comparison; actual endpoint
queries remain bounded by the existing shared residual-evaluation limit.
Mass-sensitive time accuracy detects changing humidity even when latent power
keeps liquid temperature fixed. External source history must remain immutable;
resolving films does not itself evolve or conserve a finite atmosphere.

`evaluate_with_vapor_enthalpy` extends the same original network residual with
an explicitly owned vapor-stream enthalpy excess above each cell's `h+L`.
The complete finite map has exactly the selected cell identities and may contain
positive, negative or zero excess. Its reduced balance is
`m0*(h1-h0) + dm*(L+excess) = dt*P1`; `SurfaceStep::latent_energy_j` remains
`dm*L`, while the complete vapor receipt includes `dm*excess` exactly once.
This admits nonisothermal condensation carrying upstream gas sensible enthalpy
without misreporting it as latent heat. The gas/material owner must establish
source compatibility and debit the opposite full stream. Native exposure,
represented mass, source errors, two original heat passes and local/global
acceptance remain shared. A distinct corrected trial/audit retains every
`dm*excess` before large caloric-reference sums; it cannot be consumed by the
older prescribed-film timeline account. This event-free ingress supplies no
gas equation of state, gas reconstruction or separate time-accuracy claim.

Heating can turn initial equilibrium or condensation into evaporation. An
exhausted endpoint now uses the larger initial/endpoint frozen evaporation
rate for its finite retry estimate. The estimate remains advisory: a retry
must resupply heat for its new interval and satisfy all original source,
inventory and energy gates. Nonfinite residuals or combined energy allowances
refuse. The scalar solve retains bounded work and allocation-free evaluation.

Six scalar regression groups cover direct latent heat supply, signed flux
reversal with temperature-dependent latent heat, boiling-trial narrowing,
zero-heat equivalence, caloric impulses, inventory retry, finite controls,
source/cancellation refusal and independent time refinement. Three native
cell groups cover simultaneous heating, original mass/energy receipts,
inventory constraints and final-gate rollback.

`LiquidCellExchange` and `apply_liquid_exchanges` admit externally solved
finite-atmosphere endpoints on those same cells. Full expected-cell equality,
unique identities, original caloric/storage and cell constraints, admitted
pure-liquid enthalpy range, positive remaining inventory and interval area are
checked for the complete batch before its final cancellation gate and writes.
The external physical owner remains responsible for joint liquid/gas balances.

The caller admits the solver's specific-enthalpy bracket as a pure-liquid
branch. Original births and every thermal substep are checked against that
branch; an internal-energy law refuses open-mass evaporation. The original
caloric law supplies all temperatures, including during the phase solve.
Signed transfer permits condensation without clamping humidity. Actual mass
change, full carried vapor enthalpy and thermal storage close one whole-window
account. Latent energy is separately visible only as a subaccount of vapor
enthalpy; it is never added as a second thermal source. Numerical balances use
the unchanged absolute whole-window joule allowance.

During condensation, a trial at the upper caloric bracket may exceed the film's
boiling pressure even though the physical endpoint is subboiling. The original
surface solver now narrows that trial within the admitted liquid chart until
it obtains an actual opposite-sign residual. It neither treats boiling as a
residual nor clamps temperature. Initial boiling, source refusals, cancellation
and the shared source-evaluation budget retain their original semantics.
Native evaporation regressions compare broad and tight brackets against the
independent mass/energy equation, including roots near the boiling boundary.

G0 unit regressions cover latent cooling and condensation warming, actual mass
and energy conservation with conduction and heat sources, event exposure,
birth timing, original branch/storage refusal, cumulative work limits and late
cancellation. A G1 regression refines the complete heat/phase split, which is
first order. These checks do not establish experimental model accuracy.
The atmosphere is prescribed: its temperature, pressure and vapor fraction do
not evolve here. Finite atmosphere energy/species evolution requires its own
owner accepting the complete transfer. Geometry is fixed reference geometry;
recession is reported, with no inferred VOF, liquid-flow, dryout, boiling,
freezing or contact/topology transition.

## Caloric identity for constitutive laws on native stored energy

`ThermalLaw::caloric_identity` optionally exposes the immutable source-law
identity, including its specific-energy datum, as 32 bytes without adding a
hashing or material dependency to this thermal kernel. `None` means that the
law has not supplied this binding. `DepositionSimulation::material_caloric_identity`
queries the actual registered owner and refuses unknown material IDs.

Constitutive consumers indexed by native specific energy can require this
identity together with the existing material and storage-basis checks. The
identity supplies no new caloric model and does not authenticate empirical
data. Existing thermal laws and temperature-only consumers remain valid when
the optional identity is absent; consumers requiring it must refuse absence.

## Physical source transitions on the observed native timeline

`DepositionObserver::next_transition_time_s` optionally supplies an absolute
physical source/state transition. The explicit native owner queries it after
same-time authored events, before preparing an active-cell interval. Returned
times must be finite and strictly future, even when beyond the requested target.
The native interval ends at the earliest event, transition or target after all
conduction and observer timestep restrictions. A later harmonic decay limiter
therefore cannot turn a source edge into an asymptotic sequence of steps.

Transitions are not authored events and do not consume their budget. They retain
the existing substep budget and transactional candidate semantics. Coupled
observers validate each child's result before selecting the minimum; generic
thermal, optical, flow and mechanical wrappers forward the hook. Default None
preserves existing stepping arithmetic, and the existing empty-grid skip is
unchanged. Focused observer regressions check exact coupled boundaries, invalid
time rejection with rollback, coincident native events and split continuation.

## Observed sources and advection in implicit conduction

`ImplicitConductionChain::advance_observed_to` and the corresponding graph
entry point accept original `DepositionObserver` source owners on the same
native timeline. The existing backward-Euler residual, analytic conduction/
contact/cooling Jacobian, Newton line search, chain direct solver and graph
Krylov solver remain authoritative. There is no second thermal integrator.
Original source-free nonadaptive advancement shares this event/retry path.

Start-state source watts and final-interval source watts enter the original
residual exactly once, alongside authored boundary heat. The source sees the
accepted interval START state; its supplied watts remain fixed during Newton.
Conduction and native boundary cooling are evaluated at the endpoint. This is
an IMEX discretization for temperature/enthalpy-dependent absorption, not a
fully endpoint-coupled optical solve. Optical temporal refinement remains
necessary even when implicit diffusion permits a much larger thermal step.

Native events and source transitions cut intervals exactly. Observer timestep
restrictions remain active. A nonlinear interval retry re-evaluates the source
for its actual shorter duration, so a linear ramp retains the correct mean.
Source `begin` runs once per advance; work used by rejected preparations stays
charged. Preparation may replace pending sources/caches but must not publish
physical history. The source's acceptance hook runs only after the original
thermal endpoint passes its cell and total joule gates. The entire thermal and
optical candidate remains transactional if a later source check refuses.

Conservative observer transfers use the original paired heat-rate owner inside
that same residual. `Advected` therefore composes balanced mass-carried enthalpy
with endpoint conductivity, contact and boundary cooling. Source, transfer and
timestep hooks retain their original order; one prepared flow epoch survives
nonlinear retries and is accounted only for the accepted duration. Advection is
explicit and retains its original timestep, caloric-datum, storage-basis and
fixed-mass admission. Internal transfers cancel from external heat, including
on latent plateaus; inlet/outlet enthalpy retains its original external account.
This introduces no moving geometry, momentum equation or pressure-work model.

Nonempty observer link, ambient or reservoir conductances still refuse because
the implicit solver owns those coefficients. Authored boundaries are restored
before native solving and acceptance. Native external heat already contains
observed sources, so component accounts are parallel descriptions of the same
energy, not additional credits. Generated-source and transfer entropy retain
their original unclassified status; an opted-in advection entropy account keeps
its own flow/mixing semantics. The adaptive step-doubling path does not accept
observed source owners; original endpoint laws below remain supported.

`advance_sdirk2_observed_to` on the chain and graph connects those same physical
source owners to the original two-stage SDIRK solver. The complete fixed-cell
boundary-event window is admitted before observer preparation. Each attempted
interval supplies its accepted-start sources, conservative transfers and
whole-interval mean heat to both native stages. Their positive quadrature
weights account for that physical heat once. The intermediate stage never
publishes observer history; the accepted endpoint invokes `step` once.
Nonlinear retries reprepare the actual shorter interval without restarting
the observer's work budget. Source CFL restrictions and exact transitions
remain active alongside the native energy, entropy and cancellation gates.
Holding state-dependent source or transfer feedback at the accepted start can
reduce the coupled method to first order; this entry does not infer a source
Jacobian, strong coupling or adaptive observer-history control.

`ImplicitConductionChain::with_heat_sources` and the graph attachment admit
immutable local `ImplicitHeatSourceLaw` objects bound to the original caloric
owner of each selected cell. Each law supplies signed total cell power [W] and
its analytic derivative with respect to original specific stored energy.
Temperature-dependent electrical or reaction-heating closures can therefore
join Fourier transport at the same nonlinear endpoint without freezing their
feedback at the interval start. Phase-dependent generation may vary with
enthalpy even when the original caloric temperature is on a latent plateau.
The caller owns source units, calibration, state bounds and the chain rule;
this seam infers no circuit, reaction inventory or cross-cell feedback.

These powers and tangents enter the existing joule residual and diagonal
Jacobian once, alongside independent authored and observed watts. The same
evaluation serves backward Euler, both native SDIRK2 stages and adaptive
step doubling. Original generated-source entropy remains unclassified. No
source evaluation advances physical history; invalid identities, nonfinite
powers/tangents, out-of-domain trials, cancellation or failed nonlinear gates
publish no candidate. Positive feedback still needs a physically admissible
endpoint and bounded nonlinear convergence. Independent two-cell solutions,
nonlinear Fourier/source finite differences, latent-energy feedback and
second-order time refinement exercise these paths.


Prescribed source histories use
`ImplicitHeatSourceLaw::state_at_specific_enthalpy_and_time` at the actual
native endpoint or SDIRK stage time. Its complete enthalpy tangent holds that
time fixed; existing autonomous laws retain their original response.
`HeatSourceSchedule` supplies bounded, continuous piecewise-linear nonnegative
power multipliers, and `ScheduledHeatSource` applies the same multiplier to
watts and the full source tangent. A schedule can be borrowed by several cells.
Its at least two strictly ordered samples declare a closed time domain:
out-of-range queries refuse, and constant tails must be authored explicitly.
An off interval still admits the original caloric/source domain. The time-free
query on a scheduled source refuses rather than inventing a physical time.

All fixed and adaptive BE/SDIRK paths clip comparisons to the next source knot,
alongside existing flow transitions and native events. Discontinuous powers
remain owned by native boundary events or interval-mean observers. The source
does not change storage, add heat credits, infer source entropy, or advance
history. Two-body triangular-pulse regressions compare conduction with the
independent modes `S'=P(t)`, `D'=P(t)-2D`, exact SDIRK pulse energy, nonlinear
source/Fourier Jacobians and temporal refinement. Backward Euler retains its
own first-order endpoint quadrature. The three focused production schedule
leaf checks and Rust syntax parsing ran in this implementation session; the
full native integration cases were added but not executed in this environment.


G1/G3/G4 checks cover independent stiff chain/graph forced solutions, implicit
Robin cooling with authored/start/interval sources, nonlinear halving of ramp
power, cumulative rejected-work limits, exact event/source cuts, split restart,
source timestep bounds, unsupported-hook refusal and late cancellation/history
rollback. Cross-crate tests additionally use actual conductor surface pulses,
state-dependent n/k response and Beer-Lambert volume deposition. Those checks
retain independent native and optical energy accounts; they do not certify
spatial, constitutive or nonlinear-optical accuracy from a small energy residual.

Focused advection regressions compare chain and graph endpoints with an
independent stiff circulating-pair solution and a large caloric-datum shift,
resolve an open latent channel with nonlinear conductivity/contact resistance,
and retain one flow preparation per accepted start across nonlinear retries.
Unbalanced flow and late source refusal preserve both original histories.

## Limited axial enthalpy reconstruction on native transport

`enthalpy::advection::reconstruction::LimitedEnthalpyReconstruction` prepares
bounded three-cell stencils on contiguous original centreline geometry. Minmod
limits physical one-sided enthalpy derivatives using actual adjacent lengths;
the signed native mass flow selects the donor's appropriate face. Each face
trace lies between its donor centre and neighboring sample, and an extremum
has zero slope. Missing stencils preserve the original donor-cell trace.

`Advected<()>::endpoint_transport_with_reconstruction` first runs the original
endpoint admission and source calculation, then corrects each selected internal
transfer with the same power on both cells. Original component-datum handling,
fixed mass continuity, caloric laws, external stream enthalpy and backward-Euler
energy accounting remain authoritative. Prepared mass, geometry and material
bindings are checked on every candidate. Topological junction gradients and
external streams on stencil centres refuse. Stationary species, selective
streams and explicit-trajectory entropy histories require separate closures;
the bulk thermal entropy limiter below supplies a closure for stateless and
native implicit graph queries.

The plan bounds authored stencils and retains only native cell identities and
prepared geometric factors; no new physical timeline or accepted-state clamp is
introduced. The caller solves the nonlinear endpoint residual and retains the
original cell and total joule acceptance gates. Regression sources exercise
nonuniform affine traces, smooth periodic face-error refinement, reversal and
datum invariance, an independent implicit cell balance above explicit CFL, sharp
extrema, invalid ownership/topology and cancellation. These checks concern
spatial face reconstruction and stated discrete equations, not a blanket TVD,
global convergence or arbitrary-source maximum-principle claim.

## Integrated material-to-bath fluxes in implicit enthalpy transport

`ImplicitConductionChain::with_boundary_laws` and the graph equivalent attach
complete original-law fluxes to existing reservoir and ambient faces. The
native backward-Euler owner still applies the signed heat once, stores the
original mass times specific energy, and accounts for entropy at the actual
bath temperature. Authored boundaries, their separate covered/exposed areas,
event times and caloric objects remain unchanged.

`ImplicitBoundaryLaw` supplies a passive secant and the analytic derivative of
the complete outward power with respect to specific enthalpy. The original
Newton operator uses its negative derivative in the heat-into-cell residual,
including the original cell-mass factor. Unselected channels preserve their
bulk-conductivity behavior. This seam admits exact source-integrated half-cell
conduction and a film jump without freezing conductivity at the bulk state or
adding another surface heat source. Material owners enforce support at the
actual material face; an external bath is not itself a material sample.

Attachment rejects duplicate or missing channels and substituted caloric
objects. Every endpoint invokes the same source law; event changes therefore
refresh the physical bath without changing source ownership. Cancellation,
invalid source states and failed convergence publish no thermal candidate.
Independent quadratic one-cell and Jacobian regressions exercise simultaneous
hot/cold baths, unequal cell mass, chain/graph agreement and native heat/entropy
closure. This remains a fixed-geometry, first-order thermal reduction; no
spatial convergence or empirical material-accuracy claim follows from its
algebraic energy residual.

### Nonlinear ambient radiation on the original enthalpy graph

`boundary::GrayAmbientBoundary` supplies a complete exterior heat flux to the
same implicit boundary seam. A bounded scalar solve resolves the exterior
surface temperature through the original cellwise `k(h)`, centre-to-face
distance and optional coating/contact resistance. Convection and gray radiation
then act in parallel at that surface. The explicitly fixed or source-bound emissivity and htc
are used once; authored htc and another radiation source must not also own the
selected channel. Air and the large black radiative surroundings share the
actual bath temperature supplied by the original boundary, including its events.

Factored fourth-power differences retain near-equilibrium exchange. The
allocation-free analytic tangent includes both `dT/dh` and `dk/dh`; phase-driven
conductivity changes can affect heat flux even at a latent plateau. The
original chain/graph solver, enthalpy storage, nonlinear convergence controls,
adaptive time owner and heat/entropy accounts remain authoritative. A surface
response exposes temperature, separate convective/radiative powers, applied
power and closure residual; these are descriptions of the same transfer, never
extra energy credits. Invalid data, exhausted surface work and cancellation
return no response or accepted thermal candidate.

This uses the native cellwise conductivity approximation, not an integrated
within-cell `k(T)` profile. It assumes an opaque gray surface with a complete
view of one isothermal surrounding, fixed area and no surface heat capacity.
Normal half-thickness may follow the explicit material-density closure above.
Separate radiative/convective bath temperatures and original-enthalpy surface
response use the closures below. Finite enclosure reflections and spectral
photon entropy require their own closures.
The existing thermal entropy account describes bulk/bath entropy exchange,
not a spectral radiation field. G1/G3 regressions independently reconstruct
Stefan-Boltzmann plus film heat and the material tangent, exercise sharp latent
storage on chain and graph paths, and refine adaptive cooling toward a separate
continuous calorimetry integral. This is numerical verification, not material
validation. The boundary relation follows [MOOSE's radiation formulation](https://mooseframework.inl.gov/docs/site/source/components/HSBoundaryRadiation.html).

### Distinct convective and radiative baths on one enthalpy face

`GrayTwoBathBoundary` resolves one common exterior surface through the original
half-cell/coating resistance. The selected native boundary supplies actual gas
temperature and its admitted events; an explicitly prescribed, fixed radiation
temperature supplies the large black surrounding. Explicitly fixed or
source-bound htc and emissivity act once at that surface. The analytic total enthalpy tangent retains both `dT/dh`
and `dk/dh`, including phase-driven conductivity on a latent plateau.

`with_boundary_exchanges` carries the two signed outward powers and their
actual bath temperatures through the original rate and entropy owners. The
new `BoundaryExchangeStepCoefficients` ingress preserves the existing public
coefficient structs. A selected channel cannot also have a legacy conductance
or boundary-law owner. Reference coverage, mass, caloric identity and authored
sources remain authoritative; no second half-cell or generated-heat source is
introduced. The heat account uses `-sum(Q_out)`, reservoir entropy uses
`-sum(Q_out/T_bath)`, and transfer generation uses
`sum(Q_out*(1/T_bath - 1/T_cell))` at each actual quadrature point.

Opposing bath powers therefore survive zero net body heat. A zero-conductivity
bulk still permits exterior exchange between its baths while its body tangent
is zero. Net surface-balance tolerances do not scale with large cancelling
gross bath powers. Backward Euler and native SDIRK retain the same energy and
entropy stage weights, including the current left-limit gas-bath event rules.

Nine new native regressions cover independently prescribed surfaces and
tangents, latent storage, zero-net-heat bath entropy, insulating bulk, separate
SDIRK stage references, gas-bath events, partially covered ambient faces with
internal conduction, and transactional ownership/refusals.
These and the five prior radiation regressions passed in a focused build of
the complete production thermal library and original numerical dependencies.
This remains fixed-area, cellwise-conductivity, gray full-view exchange
without surface storage or spectral radiation entropy. Normal half-thickness
may follow the explicit original-density closure above.

### Original enthalpy-dependent gray surface response

`ImplicitGraySurfaceLaw` binds emissivity, convective htc and their complete
original-enthalpy derivatives to the same caloric object as the body.
`GrayAmbientBoundary::with_surface_law` and
`GrayTwoBathBoundary::with_surface_law` replace zero fixed coefficients with one
source owner. Finite source bounds, emissivity in [0,1], nonnegative htc and
finite tangents are required. Each endpoint resolves the source once at its
original bulk h, before insulation/zero-area shortcuts, and holds those
coefficients fixed during the bounded common-surface solve. No extrapolation,
clipping, apparent heat capacity or additional half-cell is introduced.

This is an explicitly cell-state surface closure: the source owns its
temperature/phase chain rule at bulk h. Solved exterior temperature is not
substituted into source coordinates. A source needing distinct surface phase,
coating evolution, changing flow/pressure or another evolving coordinate needs
a separate closure. Source evaluation retains its own costs; the native face
kernel remains allocation-free.

Writing `a=1/(d/k+R)`,
`F=htc(h)*(Ts-Tgas)+epsilon(h)*sigma*(Ts^4-Trad^4)`,
`Fs=htc+4*epsilon*sigma*Ts^3` and
`Fh=htc_h*(Ts-Tgas)+epsilon_h*sigma*(Ts^4-Trad^4)`, the same solved balance
`a*(T-Ts)=F` gives
`Qh=A*[Fs*(a*Th+ah*(T-Ts))+a*Fh]/(a+Fs)`.
Factored radiation differences retain near-equilibrium flux. The direct `Fh`
term survives a latent plateau and zero-coefficient onset without dividing by
emissivity. A disconnected bulk retains zero total body tangent while the two
actual baths can still exchange heat.

The original chain/graph, backward Euler/SDIRK stages, covered/exposed face
areas, enthalpy storage and separate bath entropy owners remain unchanged.
Five focused regressions independently solve the exterior balance and caloric
endpoint/stages, verify source query count, latent/onset derivatives, ambient
counterflow and actual free area, fixed-coefficient equivalence, source
identity/domain/ownership refusals and cancellation atomicity. The complete
production thermal library compiled cleanly and all 29 selected native tests
passed; this is numerical verification of the declared closure, not empirical
qualification of a material source.

### Source response at the actual exterior temperature

`ImplicitGrayExteriorLaw` and `ImplicitGrayExteriorState` distinguish original
bulk energy e, solved exterior temperature Ts and actual gas temperature Tgas.
The original caloric object determines whether e is h/Cp or u/Cv. The embedded
`ImplicitGraySurfaceState` energy derivatives are partials at fixed Ts/Tgas;
the additional temperature derivatives hold e/Tgas fixed. Sources return both
partials at the same physical query, retain their own gas/coupled-domain
admission, and must not include the boundary's implicit surface response.

These carriers support the source-integrated material boundary's
`with_exterior_law` closure. They do not alter the existing bulk-only
`GrayAmbientBoundary` or `GrayTwoBathBoundary` contracts. Exterior temperature
bounds need not contain either bath or the bulk; the consuming boundary
intersects them with its physical surface bracket. Nonnegative passive
coefficients may have signed differential response, so consumers must admit
their combined material/surface elimination rather than infer monotonicity.
The original finite caloric storage, heat/entropy accounts and time solver
remain owners; no new surface state, flow state or clock is introduced.

### Two finite caloric coordinates and two actual surface temperatures

`ImplicitSurfaceExchangeLaw` supplies one passive interface flux and separate
partials in the ordered original energies and the two actual material-face
temperatures. Each energy retains its original h/Cp or u/Cv meaning and exact
caloric object. The source owns its complete interface, original coefficient
data, external conditions and coupled support. Source bounds do not authorize
extrapolation or substitute a face temperature for the body's stored energy.

`ImplicitSurfaceExchangeState::from_conductance` forms j=G*(Tsa-Tsb) and
its full product-rule derivatives from source-owned G partials, preserving
zero-conductance onset. Direct state construction remains possible for an
original law that already supplies the flux and its derivatives; the consuming
material layer verifies passive flux/secant consistency at the actual query.
These carriers introduce no surface storage, geometry/view inference or time
integration. The material adapter eliminates its original conductivity
profiles, combines explicit source and caloric derivatives, and enters the
existing `ImplicitLinkLaw` owner. Original finite heat/entropy accounting and
transactional chain/graph stepping remain authoritative.

## Implicit upwind enthalpy with native conduction

`ImplicitConductionChain::with_implicit_advection` and its graph counterpart
bind the original `Advected<()>` mass-flow graph into the existing endpoint
residual. Every Newton trial evaluates actual endpoint enthalpy, including
latent energy at unchanged temperature. Native inlet/outlet enthalpy, common
caloric references, fixed-mass continuity and internal-transfer cancellation
remain authoritative. Source heat and stream enthalpy are each credited once.

The analytic directed tangent is compiled once in joule coordinates, including
the original component-reference derivative under admitted continuity roundoff.
Tridiagonal transport retains the direct chain solver; other topologies use the
existing nonsymmetric sparse FGMRES owner. No dense matrix or numerical caloric
derivative is introduced. The selected backward-Euler or native SDIRK2 term has
no explicit flow CFL restriction; the same endpoints also reach their existing
adaptive h/T/phase controllers. Other observed sources retain their own timestep
restrictions.

This binding is an immutable snapshot of ordinary bulk enthalpy transport or
the original `LiquidEnthalpyCarrier` on a stationary pure-substance skeleton.
Liquid transport retains mixture enthalpy for storage and carries
`h_liquid=max(h,h_liquidus)` on its admitted mobile domain. The same compiled
directed coefficients are multiplied by each CURRENT trial's analytic carrier
slope: zero on the latent plateau and one on the fully liquid branch. The
liquid-side generalized derivative is selected at the liquidus. These dynamic
columns include the component-reference derivative and original cell mass;
they are never added on top of a bulk tangent that would advect solid energy.
Initial fully solid flow refuses. A Newton trial leaving liquid support returns
a material-domain exit for the existing line search; no accepted trajectory
crosses that support boundary. Source identity, fixed density, pure-liquid
inlet and conservative endpoint admission remain owned by the original carrier.

The same implicit binding admits frozen stationary-species graphs, including
inlet-fed recirculation and mixtures with different heat capacities and formation
enthalpies. The existing species solver fixes each cell's composition before
binding. Every residual retains its exact composition-aware caloric law, common
enthalpy reference, inlet-law support and local mass/species admission. The
compiled directed enthalpy tangent is unchanged; each cell's original nonlinear
law still maps its own stored h to temperature. Convection and conduction use
the same BE/SDIRK stages without an explicit flow CFL. This does not evolve
species, invent reaction heat, or change density. Chemical mixing entropy is
available only through original source anchors in the explicit shared reduced
convention described above; caloric data alone never supplies it.

Selected outlets are supported by the same original fixed-h endpoint owner.
Their exact tangent differentiates the component reference but never the
prescribed stream enthalpy. Analytic BE/SDIRK tests exercise this distinction
beyond the explicit flow CFL, including a downstream outlet whose reference
cell differs from its own cell and canonical reversal of connection orientation.
Flow-entropy histories still refuse as an implicit endpoint binding.
Thermal entropy leaves mass-carried
transfers unclassified unless the caller selects
`with_implicit_advection_entropy` after binding an admitted bulk or liquid
upwind flow, a source-anchored stationary mixture including source-bound
selective outlets, or ordinary bulk reconstruction with its thermal entropy
limiter enabled.
No phase mobility, moving inventory, pressure work or momentum feedback is
inferred. The selected method's time error and upwind spatial diffusion still
require independent refinement.

The entropy opt-in reuses the explicit flow owner's original caloric entropy
differences and the actual native endpoint powers. Its
`EntropyAudit::advection` separates net inlet-minus-outlet entropy, upwind
control-volume mixing and the transport arithmetic/continuity residue. The
same immutable material law covers each ordinary flowing component. Stationary
mixtures use each stream's original law and compatible source entropy anchors.
A common relative entropy datum cancels under mass continuity; common linear
species reference constants cancel under the stationary species balance.
Upwind mixing is a spatial discretization/control-volume contribution, not
constitutive fluid dissipation. Ordinary bulk mixing is nonnegative for a
same-law monotone T(h) in exact arithmetic. Selected-stream removal and
unlike-law pressure/composition contributions remain signed.
Isothermal latent transport carries entropy without inventing mixing production.
Native backward-Euler and SDIRK2 quadrature integrates the same endpoint kernel;
analytic mixture tests check both beyond the explicit flow CFL while retaining
their exact ordinary energy trajectories and separate temporal defects.

For a stationary solid skeleton, both explicit and implicit entropy fluxes use
the same original carrier: `h_liquid=max(h,h_liquidus)` and its caloric entropy.
Mixture h and entropy remain in the native storage account. At coexistence,
`s_liquid-s_mixture=(h_liquid-h_mixture)/T_m`; balanced liquid mass cancels this
stationary-phase offset. Superheat can melt the retained solid through the
single original energy equation. No second latent heat or melting-entropy
source is added. Original caloric identity, equal phase density, fixed inventory,
pure-liquid inlet and nonzero-liquid support checks apply at every endpoint.
This extends neither phase mobility nor the admitted stationary-skeleton model.

Backward Euler samples the actual endpoint; native SDIRK2 retains each stage's
entropy with precisely the same positive quadrature weights as its enthalpy
rates. Exact caloric storage, actual-reservoir entropy and conductive production
remain in the original account. Real heat/work sources are sampled separately
before flow powers are combined, so cancelling sources retain their unclassified
status while closed flow loops cease to look like external heat sources.
The signed temporal defect includes time quadrature, algebraic residual and
roundoff; none is converted to physical heat or dissipation. Adaptive comparisons
publish only their accepted fine trajectory. Rejection and late cancellation
publish no entropy or physical prefix.

`Advected::endpoint_transport_with_entropy` exposes the same stateless bulk/liquid
transport query for callers needing instantaneous rates. No observer history
or explicit CFL step is advanced. Nonempty reconstruction without its thermal
entropy limiter refuses the implicit entropy opt-in; the original energy-only
reconstruction remains available. Full entropy bookkeeping runs only on final
audited Newton endpoints. The opted-in face bound itself participates in every
nonlinear residual and analytic tangent, as part of the actual physical flux.

Focused G1/G3/G4 regressions cover exact high-Courant chain/graph balances,
large enthalpy-reference shifts, unequal-mass latent throughflow, nonlinear
directed-cycle Jacobians, temporal/adaptive accuracy and transactional refusals.
The `fs-lbm` flow-enthalpy integration regression binds actual material-driven
SIMPLEC face fluxes to the same implicit owner, checking serial latent inventory
and inlet-minus-outlet joules at Courant 10 in one accepted interval. These are
synthetic numerical checks, not experimental validation or throughput claims.

Focused entropy regressions cover high-Courant closed loops with independent
conductive/spatial/temporal terms, open flow with independently solved BE/SDIRK
stage temperatures, unequal-mass latent transport, fine-trajectory-only adaptive
integration, cancelling genuine sources, and late cancellation/reconstruction
refusal. Liquid-carrier checks retain isothermal stationary-solid inventory,
separate conductive and liquid-mixing production during melting, compare actual
BE/SDIRK phase-crossing quadrature against independent scalar solutions, retain
the pure-liquid limit, and exclude discarded adaptive trajectories. Fully solid
flow, unsupported histories and late cancellation preserve the original state.
These are analytic numerical conformance checks, not full fluid second-law or
experimental accuracy certificates.

## Limited axial convection on the original explicit timeline

`Advected::with_reconstruction` binds a shared `LimitedEnthalpyReconstruction`
to `ReconstructedAdvected<O>`, an observer of the existing explicit thermal
timeline. `with_optional_reconstruction(None)` and empty plans delegate exactly
to the original upwind flux, timestep and audit paths. Nonempty plans use the
same actual centreline distances, signed donor selection, bounded face trace and
optional bulk thermal entropy face limiter as the native implicit residual.
Only the original advective transfer suffix is replaced; inner-observer paired
heat, sources, interval sources, conductance overrides, source transitions,
events and accepted histories retain their original owners. One represented
face power still debits and credits the original cells. The gross transport
audit observes that actual correction; external heat is never counted twice.

For cell lengths `L_previous, L_cell, L_next`, the two limited half-cell
variations use `L_cell/(L_previous+L_cell)` and
`L_cell/(L_cell+L_next)`. Both opposite traces lie between the original neighbor
means, including nonuniform cells and flow reversal. With
`lambda_f=dt*mass_rate_f/cell_mass`, pure bulk transport can be written as
`(1-2*sum_out(lambda_f))*h_cell + sum_out(lambda_f*(h_cell-delta_f))
+ sum_in(lambda_f*h_upstream_face)`. Fixed-mass continuity makes these weights
sum to one; the bound `sum_out(lambda_f)<=1/2` makes them nonnegative. Any
selected thermal entropy face limit only shortens the bounded correction.
The original explicit owner therefore adds one further outgoing-mass decay
budget to its existing combined conduction/advection restriction:
`dt<=0.8/(a+2*b)`, where `a` is the decay represented by the prior native/wrapped
observer bound and `b=max_cell(outgoing_mass_rate/cell_mass)`. The physical
`maximum_courant` diagnostic remains `dt*b`, without doubling its units or value.
Event and transition clipping still happens after this combined restriction.

The default limiter remains `AxialLimiter::Minmod`.
`LimitedEnthalpyReconstruction::with_limiter(AxialLimiter::MonotonizedCentral)`
selects a less dissipative centered trace on the same explicit and implicit
transport owners. It differentiates the quadratic through the actual three
cell centers and bounds the half-cell increment by both neighbor jumps; the
uniform-grid limit is the standard monotonized-central limiter. Bounding the
face increment directly preserves both trace intervals on strongly nonuniform
grids, where simply doubling the minmod physical slope would not suffice.
Both geometric fractions underflowing to zero retains the donor trace.
The active three-cell derivative uses fixed-size storage, with no additional
per-face allocation or new thermal state. Default minmod retains its original
tangent column order and reservation. The optional original-caloric entropy
bound still shortens this same correction, and the combined stability allowance
above is unchanged. Kernel regressions cover nonuniform affine fields, trace
bounds, active derivatives, extreme scales, and conservative signed SSPRK2 pulse
translation; native trajectory and entropy-tangent regressions also exercise
the selection. These numerical checks do not certify arbitrary coupled
nonlinear solves or experimental accuracy.

The immutable plan/flow topology, both-neighbor, junction and external-gradient
admissions run once at construction. Original active-cell, geometry and pinned
mass/material bindings are rechecked at each source epoch, and native caloric
and mass-balance admissions still evaluate the actual state. No per-step
neighbor map or plan clone is needed. The stationary-skeleton
liquid carrier uses its original `h_liquid=max(h,h_liquidus)` with slope in
`[0,1]`; it never transports stationary-solid latent inventory. Stationary
species, selective streams and reconstructed explicit entropy histories remain
unsupported. Bulk thermal entropy face limiting supplies a local flux condition,
not a finite-step entropy account. Explicit time integration remains first-order;
arbitrary nonlinear-source stiffness, temporal convergence and experimental
accuracy are not certified by the geometric face limiter. Native candidate
cloning discards thermal and observer work together on refusal or cancellation.

Focused numerical regressions exercise signed bounded pulse transport against
exact translation at matched time steps, original conductive cooling on the
combined CFL, interval heat and source-transition/event forwarding, unchanged
upwind selections, liquid/bulk limits, caloric refusal and late rollback.

## Limited axial convection in the native implicit graph

`ImplicitConductionGraph::with_implicit_advection_reconstruction` binds the
existing `LimitedEnthalpyReconstruction` to the same nonlinear conduction and
enthalpy residual. Every trial evaluates the actual endpoint face enthalpy and
its active-limiter tangent together. The original upwind/component-reference
derivative remains compiled once; sparse face corrections retain equal opposite
heat accounts and the original cell-mass conversion to joule coordinates.
No finite-difference Jacobian, dense matrix, separate energy store or frozen
starting-state reconstruction is introduced.

The same reconstruction accepts the original stationary-skeleton liquid
carrier. Its minmod samples are `h_liquid=max(h,h_liquidus)`, so an isothermal
coexistence field has zero reconstructed liquid-enthalpy gradient even when
stationary-solid fractions differ. Each active face derivative includes its
own donor or neighbor carrier slope; a neighbor on the latent plateau supplies
no derivative with respect to mixture enthalpy. These face corrections augment
the current liquid upwind tangent, including the component reference, without
adding a bulk-enthalpy tangent. Fully liquid flow recovers the original bulk
reconstruction. A fully solid geometric sample on a zero-flow edge retains
first-order donor upwind; nonzero liquid flow through solid still refuses.
The bounded face trace stays between its original mobile-liquid samples.

The ordinary minmod limiter is continuous and piecewise linear. Equal nonzero
slopes select the left active derivative; a zero or opposite-sign slope selects the flat
branch. These are deterministic generalized derivatives at limiter switches,
not claims of classical differentiability. The existing Newton line search,
interval retries, adaptive controller and local/total joule gates remain active.
The wider stencil uses the graph's nonsymmetric sparse solver; the direct
tridiagonal chain retains its upwind API. An empty reconstruction preserves the
original upwind arithmetic.

Original face-domain, geometry, fixed-mass, junction, external-gradient and
unsupported-carrier refusals apply on each trial. The frozen flow and supplied
observer histories remain unchanged. G1/G3/G4 regressions compare the complete
endpoint Jacobian against finite differences after limiter branches change,
resolve a reversed/datum-shifted latent channel at Courant four against an
independent affine endpoint balance, and refine a circulating cell-average
sine wave against its exact translation and the upwind result. Late cancellation
retains the original thermal history. Backward Euler remains first-order in
time; native SDIRK2 and its adaptive controller use both stages' actual limited
fluxes, with the smooth-solution qualifications above. These bounded checks do
not establish a global TVD, arbitrary-source
maximum-principle, experimental-validation or throughput claim.

Liquid-reconstruction regression source covers stationary-solid latent
preservation, the affine fully liquid limit on unequal cell lengths and
reversed flow, complete endpoint tangents with different donor/neighbor phases,
dry-sample fallback, and an independently solved high-Courant native graph
interval that crosses the liquidus. Native Rust execution of these new cases
was unavailable in the implementation session; separate numerical checks of
their balances and derivatives do not substitute for executing the crate.

### Original-caloric entropy limits for reconstructed bulk flow

`LimitedEnthalpyReconstruction::with_thermal_entropy_limit` optionally bounds
the existing minmod correction using the same material's original entropy
difference and temperature. For a directed internal face A to B, write
`h_face=h_A+c`, `D_A=(h_A-h_B)/T_B-(s_A-s_B)` and
`W=1/T_B-1/T_A`. Its local spatial entropy production per unit mass flow is
`D_A+c*W`. A monotone original caloric law gives the nonnegative donor-cell
gap `D_A`; the limiter retains the original correction when that production
is nonnegative within finite-increment roundoff. Otherwise it reduces only
the paired face correction to `c=-D_A/W`. It neither repairs accepted energy
nor introduces a new heat or entropy source.

The active nonlinear tangent uses `ds/dh=1/T` and each original endpoint's
analytic `dT/dh`, already supplied by the native fused caloric query. It
replaces only the correction columns on the same sparse graph Jacobian;
original upwind and component-reference derivatives remain authoritative.
The two caloric concavity gaps and the temperature jump must be resolved
before dividing. An unresolved active bound falls back to donor upwind,
using a roundoff scale built from finite `delta_h/T`, `delta_s` and correction
terms, never absolute formation enthalpy or entropy references. Equal-temperature
latent plateaus retain their original reconstructed flux and tangent without
an apparent heat capacity or singular quotient.

The actual represented power also obeys the bound: if adding the correction
to a distant physical component-reference flux rounds farther from upwind than
requested, the shared power is rounded toward its original donor value.
An unrepresentable correction retains the exact donor power and tangent.
Both cells always receive the same opposite transfer. The final face remains
inside the original caloric domain. Stationary mixtures, selective streams
and stationary-skeleton liquid carriage require separate entropy closures and
refuse this bulk-only opt-in.

`Advected::endpoint_transport_with_reconstruction_and_entropy` exposes the
same endpoint flux and classified rates without advancing a history. The
native graph can also combine the plan with `with_implicit_advection_entropy`
for backward Euler and SDIRK2. The actual represented correction power times
`W` augments the existing spatial mixing term; boundary entropy remains owned
by the original inlet/outlet calculation. Reconciliation uses the returned
transport powers, and the signed temporal/algebraic defect remains separate.
This is a local semidiscrete flux condition, not a general finite-time entropy
guarantee; it does not widen SDIRK2's passive conduction acceptance gate to
advection or unclassified heat/work sources.

Regression fixtures compare constant-Cp faces with an independent logarithmic
entropy-conservative value in both flow directions and shifted enthalpy datums,
check nonlinear source tangents by perturbing the complete returned powers,
retain isothermal latent transport, and verify smooth spatial refinement.
Distant physical component references exercise rounding and lost-tangent
fallback. Native backward-Euler and SDIRK2 cases check conservative energy,
classified spatial mixing and the separate temporal entropy balance.

### Gray radiative links between finite enthalpy bodies

`implicit::link::GrayRadiativeLink` attaches through the existing chain/graph
`with_link_laws` seam. Two equal-area opaque diffuse-gray faces with full mutual
view exchange radiation and an explicitly supplied parallel gas-gap conductance.
The original ordered caloric objects retain finite storage. Existing contact
coverage excludes these faces from ambient heating. Their nominal and covered
areas must agree.

For each side, `r_i = d_i/k_i(h_i) + Rcoat_i` precedes its radiating surface.
Coating resistances are explicit and ordered; an authored unsplit contact
resistance refuses. One bounded scalar surface solve enforces
`q = (Ta-Tsa)/ra = (Tsb-Tb)/rb = hgap*(Tsa-Tsb) +
epsilon_eff*sigma*(Tsa^4-Tsb^4)`, with
`epsilon_eff = 1/(1/epsilon_a + 1/epsilon_b - 1)` and its insulating zero-emissivity
limit. A dimensionless surface-gap fraction and factored fourth-power difference
retain near-equilibrium heat. The returned secant applies one equal/opposite
power; the original surface residual is retained without reallocating heat.

Both analytic endpoint derivatives include their original `dT/dh` and `dk/dh`,
including phase-dependent conductivity on an exact latent plateau. The
geometry-aware query also includes original density-driven `d'/d` on each
side, with the same fixed-area, fixed-gap reduction described above. The existing
BE/SDIRK residuals, stage quadrature, Newton gates, finite-body thermal entropy
and energy owners remain unchanged. Zero bulk conductivity insulates this link;
undetermined disconnected massless surface temperatures are reported as `None`.
Domain, arithmetic, scalar budget and cancellation failures publish no response.

`GrayRadiativeLink::with_gap_resistance` composes an original `ImplicitContactLaw`
with this same radiative face solve. Its strictly positive area-specific resistance
must represent only NONRADIATIVE gap transport; bulk half cells, explicit coating
resistances and radiation must not already be included. The fixed gap HTC must
be zero, attachment is unique, and the ordered caloric objects must match the
original bodies. The source is queried once per enthalpy trial at its original
ordered bulk enthalpies, with its own temperature/phase chain rules and fixed
external pressure, separation and history coordinates. Solved face temperatures
are not substituted into that source's coordinates.

The same surfaces exchange `Qgap=A*(Tsa-Tsb)/R(h_a,h_b)` in parallel with gray
radiation. Both full enthalpy derivatives retain the additional
`-Qgap*(dR/dh_i)/(R*D)`, where `D=1+ra*dq/dTsa-rb*dq/dTsb` is the positive
surface-response denominator. Normalized factors avoid squaring R or forming an
overflowing reciprocal derivative. Latent-phase resistance feedback survives
when `dT/dh=0`. Source domains and positive finite resistance are checked before
insulation shortcuts; zero R is explicitly outside this finite-gap composition.
The native contact and complete-link attachment rules reject a second owner of
the same face. No extra heat, resistance, reservoir entropy or solver is added.
Focused regressions compare independent surface/tangent and BE/SDIRK endpoint
references, recover the original nonradiative source-contact solution, and check
original source coordinates, query count, identity, domain and cancellation.

`GrayRadiativeLink::with_surface_laws` reuses the original
`ImplicitGraySurfaceLaw` admission and query owner for both ordered bodies.
Each emissivity and its complete derivative are resolved once at the body's
original bulk enthalpy before insulation shortcuts. The source caloric objects
must match exactly; finite source bounds remain authoritative. Solved face
temperatures never substitute for source coordinates. Both fixed emissivities
must be zero before attachment; source convection values and derivatives must
be zero because the common fixed or source-owned gap already owns that heat.

The same surface solve adds
`A*sigma*(Tsa^4-Tsb^4)*(d epsilon_eff/dh_i)/D` to each full native tangent.
Normalized ratios evaluate the radiosity derivative without reciprocals of
tiny emissivities or division by emitted power. A single zero-emissivity face
retains its opening derivative, including on an isothermal latent plateau.
Simultaneous zero emissivities with two nonzero onset derivatives refuse the
nondifferentiable joint corner rather than inventing a Newton linearization.
Constant source values recover the original fixed-emissivity results. Focused
tests compare full tangents and unequal-mass BE/SDIRK chain/graph heat and
entropy against independent fourth-power surface and endpoint balances.

The gray/full-view, cell-state emissivity, cellwise-conductivity and fixed-area
approximations are explicit. This does not replace the existing general-enclosure
radiosity owner or supply surface-temperature-dependent emissivity,
profile-integrated k(T), near-field radiation, participating media or photon
entropy. Focused regressions exercise
independent two-surface balances, reversal, sensible/latent tangents, near-equal
temperatures, insulating limits and original unequal-mass chain/graph BE/SDIRK
energy and thermal-entropy quadrature.

## Implicit endpoint contact coverage on the original heat owner

`BackwardEulerStep::evaluate_with_contact_patches` and its entropy variant
expose a narrow `ContactPatchEndpoint`: the resolver can inspect the current
caloric trial and replace complete expected contact patches, but cannot edit
cell storage, material objects, geometry, boundaries, clocks or schedule state.
The existing batch contact-update checks are shared with the explicit public
replacement operation. A rejected batch poisons the private trial even when
the callback catches the returned error.

Actual endpoint coverage is installed before the original heat-rate and
entropy routines run. Consequently ambient exposed area, reservoir exclusion
and link exchange use one consistent topology. Accepted patch changes persist;
temporary rate/source overrides still restore exactly as before. Events at
the exact target remain pending in the accepted left-limit state. Ordinary
public patch replacement continues to require due events to have been
processed. No extra thermal solver, swept-contact residence law, event-location
algorithm or unconditional nonlinear-convergence claim is introduced.


## Conservative material removal from existing enthalpy cells

`enthalpy::removal` provides `remove_material_candidate` and the explicit
`remove_material_before_pending_events` left-limit ingress. Both stage one private
copy of the original graph. The physical owner supplies the original material
density, removed volume, complete replacement geometry and every incident contact
patch; no density, receding pose, conductance or mass is inferred from hardness.
Stored original and remaining mass/volume must agree with that supplied density
within caller-declared arithmetic tolerances. Positive removal must change both
represented inventories, leave positive material, and retain the exact original
enthalpy-law identity and datum. Unsupported storage, exhaustion, stale cells,
missing patches and incompatible coverage refuse without changing the input.

Specific enthalpy is preserved bit-for-bit. Actual removed mass carries its
current signed `dm*h` out of the retained cells; the same amount is recorded as a
stable negative retained-energy change without subtracting large gross energies.
The audit distinguishes requested from represented volume and reports density and
mass/volume discrepancies. The next heat step therefore uses the actual reduced
mass and caller-supplied geometry. This is material transport, not latent heat or
an additional heat sink. Subsequent debris transport, phase evolution, momentum,
strain and absolute entropy remain with their physical owners.

The ordinary ingress requires due events to be processed. The explicit left-limit
ingress permits only events exactly at the current time, preserving the clock,
schedule and cursor so a coupled owner can apply heat, removal, then the original
event in that order. Five source regressions cover reduced heat capacity, signed
caloric-datum invariance, complete incident coverage, atomic late refusal and
original endpoint-event ordering. These tests were authored and pinned-formatted;
resource constraints prevented native execution in this session.

### Prescribed pump and inlet histories in implicit enthalpy transport

`AdvectionSchedule` supplies bounded, strictly time-ordered, continuous
piecewise-linear inlet enthalpies and a common nonnegative mass-flow multiplier.
Each knot names the same subset of original inlet boundary IDs; unlisted inlet
enthalpies remain unchanged and values remain constant outside the knot range.
The chain/graph `with_scheduled_implicit_advection` attachment binds the original
flow and schedule together. `with_implicit_advection_schedule` adds a history to
an already bound upwind or limited-reconstruction flow.

The existing backward-Euler/SDIRK owner samples one immutable flow snapshot at
each actual implicit stage time, then reuses it throughout that stage's Newton
and line-search work. Native sources, paired transfers, liquid carriage and
entropy evaluation remain authoritative. The full reference-coordinate tangent
scales with the same mass rates; liquid/limited-face derivatives are evaluated
from those actual stage rates without double scaling. BE, SDIRK2 and both shared
adaptive paths clip to the next knot before solving or comparing trajectories.
All original residual/linear/interval budgets, checkpoints and transactional
acceptance rules remain in force. The supplied flow's observer history is never
advanced.

Common scaling retains cell masses, directions and an admitted frozen species
solution, including at zero flow. It does not re-solve an underdetermined stopped
species graph or approximate changing composition/density. Inlet knots are
admitted using each inlet's own original caloric law and datum; actual moving
stage fluxes retain native material, continuity and source-domain checks.
At zero scale there is no boundary mass carriage and no liquid mobility or
reconstruction derivative. A stopped pump can therefore coexist with solid
inventory; every later nonzero stage still requires native liquid support and
fixed-mass balance. A schedule does not supply momentum, arbitrary independent
edge rates, reversed flow, selective streams, discontinuous jumps or a chemical
mixing-entropy closure. Temporal and spatial refinement remain independent.

Focused regressions cover analytic non-autonomous BE/SDIRK inlet/flow states and
entropy quadrature, common-datum shifts, clipping in all four time paths,
first/second temporal refinement, stopped-liquid heating and restart refusal,
the limited-face exact Jacobian, frozen stationary mixtures with distinct inlet
laws, malformed histories, shared work ceilings and late cancellation.

### Prescribed histories on the original explicit convection timeline

`Advected::with_schedule` and `ReconstructedAdvected::with_schedule` attach that
same `AdvectionSchedule` to the original explicit observer. The immutable
schedule and baseline connections are shared across candidate clones. Each
source epoch samples rates and inlet enthalpy at its actual START time, using
the same sampler as native implicit stages. Rates always scale the original
connections, so a zero-flow interval cannot erase later pump restart. Stopped
boundary streams are absent; baseline material, stationary composition,
component references and liquid-carrier bindings remain intact.

The schedule's next knot is composed with the inner observer's next transition;
invalid inner transitions still refuse. Original CFL restrictions and final
event/transition clipping choose the accepted interval. Source owners, interval
heat, inner paired transfers, authored events and accepted histories remain on
that one native timeline. Actual sampled mass rates determine the combined
conduction/advection limit and reported Courant number. Limited reconstruction
uses those same signed rates, without cloning the plan or inner observer at
each substep. Upwind entropy histories describe the actual sampled inlet laws
and fluxes; the existing nonempty-reconstruction exclusions remain unchanged.

Construction checks original inlet identities. Each native advance renews
complete inlet-knot caloric admission against its actual starting history;
native source queries still check every interpolated inlet and actual flowing
cell. Supported knots therefore do not authorize a calorically unsupported
interior sample. A failure after locally accepted work discards both the thermal
candidate and sampled observer, including its source and transport accounts.
Without a schedule, the prior arithmetic and native hooks are unchanged.

This retains first-order explicit time sampling, rather than integrating the
within-step pump or inlet history. `max_dt_s` must resolve their variation,
especially startup from zero flow; knot clipping alone is not a time-accuracy
certificate. Common nonnegative scaling does not infer momentum, flow reversal,
independent branch changes, evolving composition or changing stored mass.
`AdvectionSchedule::flow_at_time` supplies the same admitted frozen snapshot to
other native consumers without advancing state or accepted audits. Those
consumers retain continuity, actual caloric/liquid support and physical-step
acceptance; the snapshot itself is no temporal quadrature rule.

Focused regressions use independent non-autonomous calorimetry and explicit
chain balances, a datum shift, pump stop/restart through source transitions and
events, interval heating, signed reconstructed circulation, actual upwind
entropy, stopped-solid heating into mobile coexistence, and late source-domain
or cancellation refusal. These are authored native numerical checks, not
experimental validation or a coupled temporal-order upgrade.

### Second-order explicit prescribed enthalpy transport

`Advected<()>::advance_ssprk2_to` and the corresponding
`ReconstructedAdvected<()>` entry advance the original cell energy `E = m h`
with SSPRK2: `E1 = E0 + dt Q(t0,E0)` followed by
`Enew = E0 + dt [Q(t0,E0) + Q(t0+dt,E1)]/2`. Native donor or limited face
transfers, fixed-mass continuity, source-caloric and liquid-carrier admission,
conduction and original `T(h)` remain authoritative. No apparent heat capacity
or separate temperature update replaces a latent plateau. The optional
`advance_ssprk2_with_coefficients` callback supplies sources and link/ambient/
reservoir conductances through the existing `BoundaryStepCoefficients<()>`
carrier at each actual stage state; it advances no material or source history.

Both stages use the native combined conduction/advection CFL, including the
extra reconstruction decay budget. A tighter second-stage bound, or a caloric
domain refusal at a speculative Euler endpoint, halves and retries the interval.
Both forward-Euler endpoints and the final mixture must remain on the original
caloric chart. `maximum_stage_evaluations` bounds stage queries and separately
bounds rejected intervals, including predictor failures before another query;
the native maximum-substeps budget bounds accepted intervals. Time must make
representable progress. This is a CFL-controlled second-order method for smooth
solutions, not an adaptive time-error certificate. Nonlinear-source stiffness,
phase corners and limiter switches still require timestep refinement.

Native events and continuous schedule knots cut intervals exactly. The endpoint
stage uses the pre-event topology and authored boundary state; the event then
executes once before the next interval. Coefficient callbacks must respect that
state for discontinuous changes, rather than switching early solely from the
endpoint's numerical time. Every accepted mass-throughput and external-enthalpy
account uses the same two positive stage weights as thermal storage. Internal
transfers remain equal and opposite. Rejected work changes no physical ledger;
the original compensated per-cell/birth-relative final energy gate admits the
complete candidate transactionally.

Only pure `()` flow owners expose these methods. Arbitrary observer histories,
unsteady momentum, interval-integrated heat and kinetics cannot be cloned and
averaged as if they were repeatable coefficient queries. Attached Euler flow
entropy histories explicitly refuse; the existing semidiscrete reconstructed
face-entropy limiter makes no new finite-step entropy claim here. Focused
regressions cover analytic second-order latent-ring transport in both directions,
bounded and less-diffusive reconstructed latent pulses, accelerating pump CFL
retries with exact mass quadrature, forcing/birth event sidedness, and bounded
transactional stage or caloric refusal.

### Adaptive explicit enthalpy time accuracy

The same pure flow owners expose `advance_ssprk2_adaptive_to` and
`advance_ssprk2_adaptive_with_coefficients`. Each proposal compares exactly one
SSPRK2 interval with two fine SSPRK2 intervals using the shared original stage
kernel. None of these three intervals secretly subcycles. Both fine intervals
must satisfy their own stage CFL and original caloric-domain checks. The
accepted state, integrated heat, mass throughput and flow ledger come only from
the fine trajectory; there is no Richardson extrapolation of physical storage.

The mandatory specific-enthalpy norm divides the coarse/fine difference by
`absolute_h + relative * max(|h_coarse-h_start|, |h_fine-h_start|)`. It therefore
resolves latent storage at unchanged temperature without scaling tolerances by
an arbitrary formation-enthalpy datum. An optional temperature norm queries the
original `T(h)` on all three states. The fine local-error estimate uses the
existing second-order method arithmetic: a factor of `1/3` for equal halves,
or `(a³+b³)/(3ab(a+b))` for the actual representable fine spans. The bounded
cube-root controller shrinks rejected comparisons and grows accepted spans.
This is a local smooth-solution estimate; phase and limiter corners may reduce
order, and no certified global accuracy or spatial error bound is inferred.

Native events and schedule knots end a complete comparison. Births and boundary
changes execute once on the accepted timeline, after pre-event endpoint
queries. All discarded coarse/fine stages and retries share one caller stage
budget; an additional trial cap bounds comparisons, including early domain
failures. The native accepted-substeps budget counts both fine intervals.
The caller's minimum dt applies to a comparison; an exact event/target remainder
may be shorter, but a rejected remainder cannot evade the floor. Cancellation,
accuracy-floor refusal and work exhaustion return no partial candidate.
Reports expose accepted energy, actual stage work, rejected trials, normalized
local errors and a suggested next dt without hidden controller history.

Regressions exercise semidiscrete latent-ring convergence in both flow
directions, source relaxation under shifted enthalpy datums, optional original
temperature control, accelerating pump schedules with forcing/birth events,
fine-only mass/heat accounts, and shared transactional work/floor refusals.

## Finite-lifetime luminescent excitation

`luminescence::ExcitationDecay` advances an authored emitting-state energy store
under constant pump power and constant total lifetime with the analytic linear
rate solution. It reports separately integrated pump input, remaining excitation,
radiative-channel release, nonradiative heat and the signed represented-energy
residual. Radiative release precedes photon-spectrum conversion and is not an
assertion of escaped optical energy. The caller owns the emitting-state energy
datum, prompt vibrational/Stokes heat, photon energies, transport and acceptance.
No saturation, stimulated emission, intersystem crossing, measured material law
or general phosphorescence mechanism is inferred.

Small-exposure pump release uses a cancellation-free series; the existing strict
math decay owner preserves initial-store loss and representable late tails.
The prompt limit does not require lifetime-sized thermal steps. The scalar core
is allocation-free and bounded, and returns no candidate for invalid inputs or
unrepresentable required results. Its compensated residual reports arithmetic
loss instead of hiding it in heat. Identical inputs use one deterministic
operation order; no new cross-ISA evidence or interval certificate is claimed.

Adding a nonradiative rate preserves the original radiative rate while changing
both lifetime and yield. Radiative and nonradiative branch fractions are retained
independently, so weak quenching still releases heat when the radiative yield
rounds to one. `nonradiative_yield` exposes that original heat branch directly;
the recycling owner retains it for heat and positive-sink admission. A positive
branch that cannot be represented refuses instead of becoming zero heat.
`ArrheniusQuenching` evaluates an explicitly supplied
attempt rate and per-excitation barrier in SI units; its source owner must enforce
the measured temperature domain and resolve changing-temperature coefficients.
Focused G0/G1/G3 tests cover closed-form pumping and afterglow, tiny exposures,
prompt/large-lifetime limits, scaled late tails, split-step conservation, consistent
quenching, Arrhenius units and admission boundaries.


### Explicit singlet/triplet transfer

`luminescence::transfer::SingletTripletDecay` admits two positive authored
excitation energies (`E_s >= E_t`) and competing radiative, nonradiative and
intersystem-crossing rates. Its reusable fixed-interval plan delegates populations
and integrated event counts to the existing positive `fs-la` reaction/source
owner. Forward crossing deposits the energy gap as heat; reverse crossing
withdraws it from the bath. Both exchanges and both radiative channels remain
separate, including repeated visits, so delayed singlet emission cannot silently
create energy. Radiative channel joules still precede spectral Stokes conversion.
Pump rates are admitted excitation counts per second, with pump-photon relaxation
owned by the caller. The source owner must provide thermally admissible rates and
resolve changing temperatures; no detailed-balance relation is invented. This
linear dilute model excludes ground-state depletion, annihilation and stimulated
emission. Numerical admission, budgets and cancellation retain the `fs-la` owner's
meaning. Candidate energy residuals are reported without correction. Focused
tests cover independent two-state solutions, equal rates, tiny intervals, reverse
crossing, stiff afterglow, stable stores and constant-pump semigroup composition.
Admission refuses a positive ground-state sink rounded away by a unit crossing
fraction. Every candidate additionally checks compensated excitation balance
against the caller's relative tolerance using only initial/final inventories,
injection and ground-state events; reversible crossing counts cannot inflate its
allowance. Public positive channel counts and energies that underflow to zero
refuse, including counts whose later energy rescaling would be representable.

Pinned formatter/parser checks and independent physics/API review passed. Eleven
focused tests are authored; native test execution was unavailable in the shared
editing environment.


### Resonant photon recycling

`luminescence::recycling::ResonantRecyclingPlan` evolves equal-energy excitation
stores with the existing positive `fs-la` reaction/source owner. Each cell has an
authored `ExcitationDecay`; directed absorption fractions, escape and unresolved
transport close each emission partition within an explicit caller-supplied
roundoff allowance, without normalization. Feedback yields are `phi_i * R_ij`.
All receiving material cells must be represented, including unpumped receivers.
The caller enforces the common photon/excitation energy, resolves changing rates
and transport, and owns any primary pump relaxation. Optical flight times are
neglected on this kinetic clock; wavelength conversion and saturation are absent.

The prepared transient includes self-return and arbitrary recycling cycles,
constant pumping and integrated emissions/reabsorptions; it does not require a
steady state or truncate a bounce series. Perfect closed unit-yield recycling
conserves total stored energy and gains only external pump energy. The final
energy admission excludes internal photon traffic from its residual scale.
Positive losses erased by represented feedback and positive weight/channel
underflow refuse. Escape, unresolved outgoing energy and nonradiative heat stay
separate, with an uncorrected residual. `RecyclingLimits` re-exports the existing
numerical work/accuracy limits and 64-cell bound. Independent closed-form tests
cover closed self-return, two-cell exchange and lossy pumping; split-step tests
cover every integrated channel, with explicit fraction, range and cancellation
refusal regressions.

Seven focused core regression groups are authored. Pinned formatter/parser and
independent numerical/API review passed. A native Cargo test attempt stopped
during dependency compilation with `No space left on device`; no luminescence
tests executed in that attempt.

## Spectral singlet/triplet photon recycling and pair annihilation

`luminescence::transfer::recycling::SingletTripletRecyclingPlan` evolves two
real excitation states per active cell with the existing positive reaction
propagator. Spectral photon channels carry photon energies and directed
escape, unresolved, and absorption fractions. An explicitly active receiver
can create at most one singlet or triplet per absorbed photon; every positive
branch must be energetically possible. Passive receivers introduce no state.
Repeated returns, local crossing, reverse crossing, and continuous pump ages
are integrated together without a bounce cutoff or invented leakage.

The output separates source relaxation heat, receiver excess or passive heat,
nonradiative heat, forward-crossing heat, and reverse-crossing bath withdrawal.
Internal emitted and recycled traffic is reported but is excluded from the
external energy scale. All thermal channels enter the external audit as one
compensated net exchange, so balanced internal photon or bath cycles cannot
inflate its error allowance. Authored fractions are validated, never repaired.
The current numerical owner admits at most 64 states, hence 32 active
singlet/triplet cells; passive absorption does not consume those states.

`transfer::annihilation::TripletAnnihilation` supplies the isolated exact
mean-population pair step. Its coefficient is a pair-event coefficient:
events occur at `k * N_T² / V` and consume two triplets. Each event produces
one singlet with the explicitly authored probability; the remaining energy
becomes heat. A productive event requires `E_S <= 2 E_T`, irrespective of
how small its probability is. Active volume is explicit, and population
counts must not be confused with concentrations. This reduced law does not
infer diffusion, spin statistics, density saturation, or a measured material
coefficient. Numerical tests use manufactured parameters.

## Controlled nonlinear phosphorescence evolution

`SingletTripletRecyclingPlan::advance_with_annihilation` couples optional
local pair laws to the original singlet/triplet, pump and spectral photon
network. Positive exponential midpoint prediction and correction resolve
the original stiff linear kinetics together with frozen nonlinear rates.
Only corrector histories enter physical accounts. Original nonradiative and
reverse-crossing channels are reconstructed directly from their original
rates; internal reaction coordinates introduce no reported fictitious bath
exchange. Pair counts, productive singlets and pair heat remain separate.

Inactive laws retain the exact existing linear path. An isolated unpumped
pair law uses its analytic solution. General evolution refines bounded
meshes, compares stored energy and integrated physical channels, and also
checks endpoint radiative power expressed over the interval. A true-time
geometrically refined final segment resolves fast emitting states without
uniformly resolving their lifetimes; its work counts against the authored
budget. Identical meshes provide no convergence evidence. A nonlinear
exposure guard prevents an over-damped predictor from falsely appearing
converged. Budget exhaustion or cancellation returns no accepted candidate.

`AnnihilationLimits` controls componentwise absolute/relative convergence
and maximum segments. The reported discrepancy is a refinement indicator,
not a rigorous error enclosure. Whole-interval conservation independently
uses the original linear owner's relative allowance and compensated net
thermal exchange. Linear coefficients and optical transport stay fixed
within the interval; thermal feedback still requires outer time refinement.


### Current conduction measures on retained material inventory

`ConductionObserver::with_current_geometry` retains original cell, material,
mass, reference geometry and channel bindings while a complete ordered
`CurrentConductionGeometry` supplies current link, reservoir and ambient areas
and positive normal half-distances. The original series resistance, signed
flux, entropy, source refresh, CFL, event and transactional owners consume those
metrics; no caloric or density remap occurs. The independent mechanical or
geometric owner supplies these measures, held for one proposed interval.

Implicit chains and graphs accept the same current measures in backward Euler
and SDIRK2. Their fluxes, conductivity/contact tangents, finite insulating-onset
derivatives, and fixed or scheduled reservoir/ambient derivatives all use the
same current areas and normal distances. Newton trials and both SDIRK stages
hold this mechanical configuration fixed while the original nonlinear material
laws continue responding to enthalpy. The sparse/direct solve, conservative
storage, paired heat, bath entropy and transactional acceptance remain shared.
Original mass and contact coverage are rechecked both at construction and on
continuation, so a stale deformed binding cannot admit a changed inventory.

Complete link, boundary and two-bath laws must explicitly opt in through
`supports_current_geometry`. Its default refusal protects reference-pinned
adapters before any flux query. Compatible link laws receive current flux area
and normal distances with original ordered cells and per-side reference
coverage; boundary laws receive actual area/distance in their existing geometry
argument. Native gray closures opt in while retaining their geometric domains,
including equal-area/full-view restrictions on radiative links. Source-owned
admission may still reject incompatible measures; reference geometry is never
substituted. Unbound ambient Robin channels retain authored reference areas.

Normal two-point flow, unchanged topology and coverage, and first-order
physical-interval/substep refinement remain explicit assumptions. Arbitrary
skewed-tetrahedron Fourier accuracy, opening contacts, remeshing and anisotropy
are not inferred. Phase-density geometry cannot silently compose with these
metrics. Removing the explicit thermal CFL does not resolve mechanical motion
or make the frozen-geometry coupling second order.

Two synthetic regression definitions cover independent mixed-channel heat,
reference inventory, native CFL/entropy, and failed-admission retry. Seven
implicit regressions define unequal-capacity analytic BE/SDIRK relaxation,
branched finite-body balance with asymmetric reference coverage, nonlinear
contact and scheduled-bath Jacobians, insulating onset, integrated Kirchhoff
links/baths, separate bath entropy, compatible gray closures, reference-pinned
source refusal and late-cancellation retry. Focused Rust
parsing/formatting passed; native compilation and integration-test execution
remain unverified in this environment because workspace dependencies are
unavailable.

### Signed heat channels and latent storage

The shared enthalpy rate owner retains Neumaier compensation for each cell
and for external heat through conductive links, separate reservoir/bath
contributions, authored and generated sources, and directed internal transfers.
Only the consuming energy equation rounds the resulting net power. No extra
heat, modified caloric curve, or apparent heat capacity compensates for lost
arithmetic. The correction requires one additional scalar per cell and no
per-channel allocation; the existing sparse traversal and cancellation gates
remain authoritative.

Backward-Euler local residuals and storage likewise use the actual per-cell
`m*(h_new-h_old)` and retained powers before rounding aggregate report fields.
Opposite local defects cannot pass solely because their total vanishes. Explicit
observed advances retain one compensated external-energy ledger across substeps
and bound absolute losses when final cell powers cannot represent a low component.
Equal rounded storage and input reports do not certify a smaller physical
residual. Existing interval-source external-ledger checks remain active.

Source entropy uses these same authored, generated, interval and directed
transfer operands, retaining their correction through division by each cell's
temperature. A latent heat increment contributes `Q/T` without inventing
conductive production or a temporal defect from premature source rounding.
Manufactured regressions cover weak latent heating and entropy between opposing
large channels, reordered link IDs and transfers, generated heat, unequal-mass
inventories, and heat retained across a heating/cooling history.

### Continuous thermal impedance on the native enthalpy graph

`ImplicitConductionChain::prepare_thermal_impedance` and the graph delegate
admit an original per-cell steady heat balance, then retain its complete native
conductive/contact/boundary/source tangent. With peak `exp(-i omega t)` phasors,
the reusable frequency factor solves `(K-i*omega*I) delta E = delta P`, where
`K=-d(net inward watts)/dE` and `E=m*h`. Original caloric derivatives recover
temperature and driving/mutual impedances in K/W. Nonlinear conductivity,
source feedback, integrated surface laws and their existing geometry/contact
derivatives remain in the same Jacobian. No timestep, fitted thermal capacity,
diagonal-only transfer or second material model is introduced.

An interior latent plateau retains its energy perturbation with zero temperature
perturbation; no infinite apparent capacity is formed. DC analysis refuses
unanchored energy-conserving components and singular tangents instead of adding
artificial leakage or selecting an arbitrary temperature reference. Dense work
and numerical payload are explicitly bounded; each frequency shares one native
complex LU across drives, checks original-row residuals and polls cancellation.
Preparation queries original rates without cloning the future event schedule.
No result advances accepted heat, geometry, material state or physical time.

The scope is an equilibrium conduction linearization with external clocks,
baths and prescribed geometry frozen. Advection is not admitted. The caller
must supply locally differentiable constitutive branches and perturbations
that remain within them; a one-sided corner tangent is not a two-sided material
response. Nonlinear stability, positive-realness with active source feedback,
experimental accuracy and continuum error are not certified. The quadratic
temperature/heat-rate pairing is in K W, not a second physical heat account.
Focused G1/G3/G4 regressions compare scalar RC and two-cell mutual impedances,
the passive quadratic identity, native time-domain refinement, nonlinear
source/conductivity feedback, latent storage, DC refusal and cancellation.

## Finite ground-state populations in luminescence

`luminescence::transfer::saturation::SaturableSingletTripletPlan` extends the
existing authored singlet/triplet law with an explicit ground-state inventory.
The original positive `ReactionPlan` evolves conserved mean counts in
`[ground, singlet, triplet]` order under constant per-ground pumping hazards.
Only actual ground-to-excited transitions debit the supplied pump-event energy.
Energy above the destination state is prompt heat. Radiative, nonradiative and
forward/reverse crossing channels retain their original units and meanings.

Excited-population, total-emitter and external-energy balances are checked
independently. The large ground inventory does not inflate the excited-count
allowance, and net bath exchange prevents reversible crossing traffic from
inflating the energy allowance. No clipping or residual-to-heat correction is
applied. Prepared propagators are reusable; native work, accuracy and
cancellation limits remain explicit.

This is a prescribed-field rate model, suitable for an explicitly justified
optically thin or externally maintained field. Optical pumping requires
authored cross sections and local photon flux. Fixed dilute Beer absorption
must not be relabeled as a saturable pump. The model does not solve pump
depletion, carrier transport, stimulated emission, excited-state absorption,
photobleaching, trap distributions or pair annihilation. Radiative energy is
released before the original spectrum, transport and Stokes partition.

Seven focused G0/G1/G3 tests in the module compare independent analytic
saturation, triplet shelving, a finite stable store, pump-off equivalence,
all-channel split-step composition, nanosecond/SI population scaling, and
admission/cancellation/resource failures. These seven tests passed using the
unchanged production math, reaction and thermal source files with the pinned
Rust toolchain. Full-workspace and renderer-integration gates were not run
successfully in this session.

## Photon-conservative depletion of a finite-population pump column

`luminescence::transfer::pump_column::SaturablePumpColumn` couples the original
ground/singlet/triplet kinetics to one ordered, collimated spectral pump. Each
lane can excite both destinations; its photons compete between those channels
and independently declared passive host absorption. Molecular depth is the
original cross section times ray length divided by cell volume, multiplied by
the complete cell's instantaneous ground inventory. No density multiplier or
invented beam area replaces that inventory.

A bounded scalar solve matches the Beer field to exact integrated ground
residence from the shared stiff kinetic propagator. Whole-column temporal
refinement resolves bleaching and the correlation of upstream transmission
with downstream exposure. Only the accepted finer history enters the returned
accounts. Actual pump transitions, independently transmitted photons, passive
heat, stored excitation and emitted state energy are audited without clipping,
renormalization, extrapolation or residual-to-heat correction. Wavelength-order
pump energies remain separate from destination-order kinetic energies, including
weak higher-energy channels whose relaxation heat is smaller than one ULP of
the mean pump energy. Work, accuracy and cancellation bounds are explicit.

This is a spatially homogeneous cell closure with quasistatic photon flight and
fixed constitutive coefficients during each call. Spatial and native thermal
interval refinement remain separate requirements. Pump scattering/reflection,
stimulated emission, excited-state absorption, multiple overlapping beams and
active recycling of emitted photons are outside this column's scope. Thermal
and optical owners supply original geometry and refresh their own material state.

Eight focused native tests cover exact ground residence, passive spectral
transport, the independent saturable-slab fluence solution, two-cell bleaching
and pulse composition, competing spectral/destination channels, weak energetic
branches, steady recovery and bounded work/cancellation. A recovering absorber
retains its independent steady transmission across one million radiative
lifetimes in two temporal substeps, while absorbing over 691,000 photons per
emitter without limiting cumulative excitations to the finite inventory.
These and all 64 preceding thermal
luminescence tests passed against the actual production math and reaction
sources with the pinned Rust toolchain. Renderer/material integration tests
are separate from these 72 executed lower-owner tests.


## Explicit spherical Sherwood conversion

`SphericalDropletGeometry::mass_transfer_coefficient_at_sherwood_m_s` converts
separate molecular diffusivity and diameter-based Sherwood to `Sh D_AB/d`.
It supplies geometry only: correlation selection, source support and physical
applicability remain with the transport caller. Positive finite Sh and
nonnegative finite D are required; positive unrepresentable conductance refuses.
Sh=2 delegates to the original D/r method exactly, and D=0 remains sealed.
A focused regression specifies exact stagnant arithmetic, dimensional
conversion and invalid inputs. This new test has not been executed.

## Simultaneous beams sharing finite luminescent populations

`luminescence::transfer::pump_beams::SaturablePumpBeams` extends the same
photon-conservative cell closure to multiple incoherent or phase-averaged
spectral beams. A path lists original cell indices and its own molecular and
passive optical depths. An original cell has one kinetic inventory and one
accepted history, regardless of how many paths overlap it. All incident
channels contribute to the same ground-state pumping hazards; no beam
independently advances or replenishes a copied emitter population.

The original scalar residence solve and stiff three-state propagator are shared
with `SaturablePumpColumn`. Acyclic path dependencies use a topological sweep.
Counterpropagating and other cyclic paths use bounded physical source iteration
on the common interval-mean ground inventories. Every trial restarts from the
same initial populations. Accepted common-field residence and per-channel
photon residuals must satisfy the original kinetic tolerance. Failed cyclic
grids are discarded and temporally refined. All scalar preparations, including
failed source iterations and coarse grids, count against the aggregate work
budget. Only the converged finer physical history contributes accounts.
Pump-off and fixed passive transport use one exact history per unique cell and
skip source iteration and temporal comparison, even for cyclic geometry.

Each beam retains separate original spectral incident, transmitted, molecular
and passive photon/energy ledgers. Molecular pump events, destination energies,
weak positive relaxation-heat branches, fluorescence/phosphorescence release,
stored excitation and passive heat remain tied to the one shared cell history.
Internal segment passage is not an external energy input. Per-beam photon
closure and whole-system energy closure are audited independently; a residual
is never normalized away or deposited as heat. Empty paths and empty model sets
are supported under explicit positive admission bounds. A path may revisit the
same original cell; every visit keeps its own depth and exact position in the
path while sharing that cell's one evolving inventory. In particular, returns
separated by a passive absorber are not merged. Returning paths participate in
the existing cyclic common-ground solve, and every occurrence counts against
the original segment bound. Invalid cell indices, unsupported photon energies
and exceeded budgets are refused without mutating caller populations.

The physical boundary remains quasistatic directed transport with spatially
homogeneous original cells and fixed constitutive coefficients during each
call. Molecular depth per ground emitter is cross section times ray length
divided by the original cell volume. Spatial refinement and the caller's native
thermal/pulse clock remain separate accuracy requirements. Coherent standing
waves, stimulated emission and excited-state absorption are not modeled. This
owner does not sample scattering or reflection: the geometry caller supplies
complete ordered paths. Active recycling of emitted light is handled by the
separate coupled owner below.

Seven new focused regressions cover identical-path power splitting, original
column equivalence, unique-cell afterglow, passive spectral overlap and empty
paths, weak energetic beams, admission/cancellation/budgets, and two independent
opposed-beam analytic cases. The transient case uses distinct destinations and
photon energies; the steady recovering case spans one million radiative
lifetimes without limiting cumulative pump events to the emitter inventory.
All seven shared-beam tests and all eight original column tests subsequently
passed natively against the unchanged production math and reaction sources.
The opposed transient used 256 accepted substeps, and the stiff steady case
used two while retaining more than 691,000 pump events per emitter. Full
renderer compilation remains a separate gate.

Two further authored return-path regressions specify independent changing-ground
fluence and passive absorption for an active cell, passive neighbor and return
to the original active cell; and stationary locally maintained recycling with
separate first/return depths and native event-energy attribution. The admission
fixture checks invalid indices and counts repeated visits against the segment
limit. The scalar reference equations have been checked independently; these
new Rust regressions are unexecuted because the current workspace has no Rust,
DSR or RCH toolchain. Earlier native pass counts do not cover this extension.

## Finite-inventory triplet annihilation and delayed fluorescence

`saturation::annihilation` couples the original authored `TripletAnnihilation`
law to the conserved ground/singlet/triplet inventory. Pair events occur at
`k N_T^2 / V`, consume two triplets, create `phi` mean singlets and return
`2 - phi` emitters to ground. Productive events must pay their singlet energy
individually; active volume, coefficient and yield are never inferred.

`SaturableSingletTripletPlan::advance_with_annihilation` retains the original
constant per-ground pump hazards while resolving pair density changes. Positive
exponential-midpoint stages freeze `lambda = 2 k N_T / V` inside the same native
three-state reaction propagator as pumping, decay and crossing. Original dark
and reverse-crossing events are reconstructed from native integrated removals;
temporary reaction coordinates introduce no physical bath exchange. Pair
events, produced singlets and direct pair heat have separate fields in
`SaturableSingletTripletStep`. Its energy residual includes this physical heat.

`SaturablePumpColumn::advance_with_annihilation` and
`SaturablePumpBeams::advance_with_annihilation` use the same physical midpoint
schedule with one optional law per original unique cell. Each predictor and
corrector solves its optical field and native kinetic history with the frozen
pair rates together. The returned ground inventory therefore changes actual
pump depletion during the pulse. Pair events are not applied after an otherwise
completed optical solve, and overlapping beams never receive independent copies
of the emitter inventory. A local copy of immutable geometry is made once per
call; only the per-cell kinetic models change between stages.

Both optical and pair convergence allowances constrain the accepted finer
history. Aggregate cell-evaluation and coupling counters include predictors,
correctors, rejected optical iterations and discarded grids. The temporal cap
is the smaller of the optical and pair bounds, including the endpoint tail.
The tail's original-rate bound includes all incident beam hazards at a cell,
as well as excited-state rates, so fast ground recovery under strong pumping
does not require a globally lifetime-sized grid. Pump-off cells still evolve
nonlinear pair encounters and delayed fluorescence; only their optical source
iteration disappears. Passive photon heat and every original per-beam ledger
retain their meanings.

Only corrector histories enter the compensated result. A nonlinear exposure
guard prevents false coarse/fine agreement; bounded refinement compares event
energy, stored energy, every integrated physical energy channel and endpoint
radiative power expressed over the original interval. A graded true-time
endpoint tail resolves fast original lifetimes without forcing a uniform
lifetime-sized grid. Nonlinear encounter dynamics still require refinement.
The reported discrepancy is an indicator, not a rigorous enclosure. No field
is clipped, normalized or converted from a residual into heat.

An absent or zero pair law uses the original linear path. Isolated dark pair
kinetics use the original exact rational solution, including conserved ground
return. Limits and event energetics are checked even at zero duration. The
three independent native regressions passed: isolated branches with yields
zero, one-half and one; a pumped nonlinear equilibrium; and nanosecond-singlet
delayed fluorescence across a one-second dark interval. The latter used 542
accepted nonlinear substeps rather than resolving a billion lifetimes globally.
Three additional optical regression groups passed for independent steady
column/opposed-beam equilibrium, enhanced absorption during an actual pulse
and beam-splitting invariance, and dark shared/unvisited recovery with exact
zero-law behavior and budget/cancellation/retry checks. Together with the
fifteen unchanged column and simultaneous-beam regressions, all twenty-one
focused native tests passed against the production math, reaction and thermal
sources. Renderer execution remains a separate gate.


## Finite-ground spectral photon recycling

`pump_beams::recycling::SaturablePhotonRecycling` couples authored emitted
paths and external beams to the same original ground/singlet/triplet
inventories. Spectral/directional photon fractions sum to one per emitting
state; they are admitted without normalization. Every photon must be
affordable by its source state and each active receiving branch. Source
Stokes heat, destination relaxation, passive heat, escaped photons and
declared unresolved photons remain distinct. Photon flight time is neglected;
paths are ordered, cell-averaged inputs and may return to an original cell.
Every repeated visit uses that same cell's common ground history. The source
and receiver ledgers retain the actual segment order, including passive losses
between visits. Geometry sampling, spectral quadrature and measured material
validity remain caller responsibilities.

Each nonstationary stage solves mean radiative source rates together with
the original common-ground optical closure. Sources are bounded by radiative
rate times the original cell inventory. Safeguarded inverse-Broyden updates
accelerate return modes; rejected trials retain their aggregate kinetic work.
Emitted-source closure retains at most `maximum_secant_history` inverse-Broyden
rank-one corrections (default eight; zero selects safeguarded physical source
iteration). For n radiating states and retained history m, accelerator storage
and application are O(n m). Radiating-state counts can exceed sixty-four while
the original per-cell kinetic node limits remain enforced. A full history
restarts before forming the newest secant correction. Original cell/path/segment
quotas and representable storage sizes are admitted before intermediate
geometry, source-vector and history allocations. Source positivity, physical
closure tolerances and work accounting retain their original meanings.
Optional authored triplet-pair kinetics use the same positive midpoint schedule. Temporal checks include physical event energies,
every attributed heat channel and endpoint radiative power, including when
no pair law is present.

A stationary candidate evaluates physical event flows at the original
populations and cancels same-state self-return before forming net rates.
It checks each state's integrated mismatch against the original
initial/external-energy allowance. Accepted constant histories split the
original radiative event counts consistently through the photon paths.
This avoids subtracting large self-return counts to recover a small escape
rate. It is numerical stationary admission, not an equilibrium certificate.
The final physical energy audit still applies.

Returned external/recycled pump events, photon energies and relaxation heat
are accumulated directly from actual channel events. Internal emitted and
absorbed traffic is excluded from external input and from the global
energy-balance scale. No population is clipped, no coefficient inferred and
no residual deposited as heat. Refusal preserves caller state. Scattering and
reflection sampling, stimulated emission, excited-state absorption and finite
photon flight are not introduced by this owner.

Five focused native regression groups passed against the corrected production
sources: independent exponential-integral self-return decay; exact emitted-fluence
saturation of a finite second cell with separate spectral heat; a near-closed
pumped equilibrium retaining 200,000 radiative events per emitter in one accepted
stationary step; no-active-return parity with cancellation, work refusal and
retry; and an independent stationary recycling/triplet-pair heat balance.
The self-return transient used 525 accepted substeps, with a global residual
of -7.80e-20 J. The near-closed stationary residual was 4.14e-25 J. The mixed
pair/recycling equilibrium used one substep and one cell evaluation, with
a residual of -2.65e-23 J.

Two additional native groups exercise the scalable closure: an independent
65-cell transient recycling ring and a dense-reference comparison of bounded
rank-one inverse updates. The ring uses eight substeps and 54,080 cell
evaluations, with final S/N = 0.18471378925 against the independent E1 reference
0.1847129753070053 and global energy residual -5.81e-21 J. Its changing
populations exclude the stationary shortcut. The accelerator checks multiple
directions, newest-secant satisfaction, bounded restarts, zero history,
numerical reset, cancellation and oversized storage refusal.

Together with the eighteen existing column/shared-beam/pair regressions and
three prescribed pair-kinetics oracles, all twenty-eight focused native tests
passed with unchanged controls and no compiler diagnostics. These use the
actual production math, reaction and thermal sources. The renderer binds this
owner to original geometry, material absorption, native heat and accepted
endpoint opacity. Its separate full-runtime tests remain unexecuted;
lower-owner checks do not establish that renderer gate.

### Prescribed local excitation sharing transported photons

`pump_beams::LocalPhotonPump` declares per-ground `G->S` and `G->T`
hazards in inverse seconds and their original photon energies in joules.
`SaturablePumpBeams::advance_with_local_pumps`, its
`advance_with_local_pumps_and_annihilation` counterpart, and
`SaturablePhotonRecycling::advance_with_local_pumps` combine those fields
with actual transported photons on the same original finite inventories.
Local entries follow original model order. Empty or all-dark input retains
the existing numerical path and an empty local report.

Local hazards enter before the single native population propagation. Actual
local events are the hazard times that same integrated ground residence;
recovery can therefore support more than one absorbed photon per original
emitter without increasing its inventory. `LocalPumpStep` retains these
events, consumed photon energy and excess-energy heat separately from every
beam and internally recycled path. No weak-opacity surrogate ray, beam area,
or incident/transmitted local reservoir ledger is constructed. The caller
must justify a prescribed field as externally maintained or optically thin
and must not declare the same illumination again as an external beam.

The stationary flow admission, positive nonlinear/pair stages, fastest-rate
bound and temporal energy checks all include the same local hazards. The
physical residual includes actual local absorbed energy as external supply,
alongside initial excitation, incident beams and any uphill bath transfer.
Gross photon circulation remains excluded from external supply and the
global error scale. Local contributions are accumulated directly, without
subtracting large transported or recycled totals. Incorrect cell counts,
nonfinite/negative data, or active photons below their destination energy
refuse before state publication.

All 28 focused native pump tests passed with unchanged controls: the 25
existing column/shared-beam/recycling/pair groups and three new independent
local-pump oracles. The exact two-state transient consumed 1.0347737207 photons
per emitter while conserving its inventory, with final S/N=0.60563795311.
A shared local/beam/recycling/pair equilibrium used one substep and one cell
evaluation; the independent direct local/pair branch used two substeps and
six evaluations. Their energy residual magnitudes were below 7e-23 J.
The transient also checks inactive-field parity and invalid shape/energy
refusal. The prior trapping and 65-cell ring results and work counts stayed
unchanged. These tests executed actual production math and thermal sources
in the existing harness with no compiler diagnostics or ignored tests.

## Stationary and harmonic prescribed throughflow

The native implicit operating-point and thermal-impedance owners now admit
prescribed first-order whole-enthalpy upwind flow. Every residual uses the
original endpoint transport, continuity/material/source-datum checks, paired
transfers and external enthalpy powers. The already compiled native directed
energy tangent is reused without symmetrization or temperature-flux fitting.
Schedules sample rates and inlet data at the frozen operating clock; no event,
inventory, flow history, heat/work credit or entropy history advances.

Stationary storage must reach an active original bath or a positive original
enthalpy outlet along directed transport and bidirectional conduction. A linear
reverse traversal checks this topology; actual native factorization and every
physical heat row still decide numerical acceptance. Closed circulation and a
stopped ungrounded pump do not create an energy reference. Whole-enthalpy
throughflow can determine stored energy on a latent plateau while dT/dh stays
exactly zero; harmonic energy response is retained even when temperature
response vanishes. Inlet data and flow rates are held fixed for heat-drive
phasors; flow modulation is not inferred.

Reconstructed/limited advection, separately mobile liquid carriers and selected
prescribed-enthalpy outlets require their own stationary constraints and refuse
this path. Frozen stationary species retain their original per-endpoint source
admission; no composition evolution or momentum closure is added. Six native
G0/G3/G4 support-graph tests passed. Six full-owner G1/G3/G4 regressions cover
analytical two-way throughflow, directed DC/frequency response, latent storage,
frozen schedules and mass, closed/zero-flow refusal, budget/cancellation and
unsupported reconstruction. Those full-owner tests are added but unexecuted
here because the real runtime dependency compilation is resource-blocked; no
full-crate, continuum, stability or forward-error certification is claimed.


### Source-declared harmonic phase corners

`ImplicitConductionLaw::admit_harmonic_operating_point` lets the original
material refuse known derivative corners before ordinary or closed-component
harmonic preparation. Finite time stepping and stationary solves retain their
original one-sided constitutive tangents. The default preserves the caller's
existing differentiable-branch obligation and does not certify smoothness.
The pure-fluid PR adapter rejects its exact saturated-liquid/vapor endpoints;
interior latent points retain stored-energy response with zero temperature
response. The additional native throughflow regression distinguishes those
three cases. It is source-reviewed and parsed, not a completed full-runtime
execution claim.

## Physical boundary drives on the original equilibrium thermal response

`ImplicitConductionChain::prepare_driven_thermal_impedance` and the same
method on `ThermalOperatingPoint` compile bath-temperature, inlet-enthalpy and
common fractional mass-flow inputs around the original accepted operating
field. The original constitutive endpoint, nonsymmetric continuous energy
tangent, row-scaled complex LU and complete original-row residual gates are
shared. One frequency factor accepts many combinations of physical drives.

An inlet contributes its actual scheduled mass rate times the supplied
enthalpy phasor. A common fractional change of every prescribed connection
contributes the native net advective watts, using the original datum-reduced
fluxes and continuity-roundoff convention. It preserves fixed mass and any
admitted frozen stationary composition. Independent valve changes require
their own continuity/species solve. A stopped pump remains stopped; peak
fractional amplitudes reaching a zero-flow or reversal corner refuse.

Ordinary ambient/reservoir temperature inputs use the original conductance
at fixed cell energy. A complete nonlinear boundary without an explicit
bath-input derivative refuses that selected drive; its original state tangent
still participates in heat and inlet/pump responses. Bath temperatures must
remain positive, and retained original inlet laws admit peak enthalpy extrema.
These are branch/domain checks, not finite-amplitude linear-error bounds.

All inputs use peak exp(-i omega t) phasors and original IDs. Missing physical
inputs mean zero; heat phasors are additional prescribed watts, with nominal
source feedback already in the tangent. Compensated deterministic input
assembly checks cancellation, finite representability and explicit resources.
The added preparation requires 64*n*n+512*n+128*k numerical bytes and
64*(n+k) input work allowance beyond the original dense-work bound. No clock,
event, material inventory, flow history or heat/work credit advances.

Six G1/G3/G4 regressions specify directed combined forcing, original mass and
bath geometry, frozen schedules, latent energy with zero temperature response,
formation-datum invariance and input/domain/budget/cancellation refusals.
They are source-reviewed but not compiled or executed in this session because
the development executor is disconnected.
