# CONTRACT: fs-system

> Status: version 0. The Component contract in fs-ports is v0 until five kernel
> slices pass on it (epic `frankensim-rc-1109-root-3gc12`, beads `.6.1` to
> `.6.7`), and this driver changes with it. Physics owners that implement the
> contract so far: the gas cavity in fs-gas, the linear modal structure in fs-nlmodal, and, in
> fs-mbd, the rigid body and the damped rigid-body
> network. The crate's own `elements` are ideal lumped elements, and
> `frankensim system <scene.json>` runs devices made of them from a data file.

## Algebraic owners in monolithic composition

A registered `StageResidual` owner may declare zero internal unknowns.
Its residual and stage vectors are empty, but its output equations still
participate in the coupled interface solve, and its acceptance receipt is
audited before atomic publication. Automatic tier selection and an explicit
monolithic policy both admit these owners, including wholly algebraic loops.
Stationary analysis uses its configured monolithic solver after tier admission
even when a loop mixes stateful and algebraic owners.

A monolithic execution plan reports both the internal stage unknowns and the
retained interface inputs in its solve dimension and variable names. The
shared residual conformance check also covers zero-unknown owners: output
derivatives, cancellation, acceptance and conservation remain checked.
Regression coverage is in `tests/driver_modes_and_store.rs`.

## StageResidual on default lumped elements (bead `frankensim-rc-1109-root-3gc12.11.7`, Stage 2b)

Eleven lumped elements implement `fs_ports::component::StageResidual`:

- Stateful (`unknowns() == 1`, exact midpoint/discrete-gradient residual matching `propose` to roundoff, with analytic `jvp_z` and exact transpose `vjp_z`): `Mass`, `Spring`, `HeatCapacity`, `Storage`.
- Zero-unknown (`unknowns() == 0`, residual outputs and `accept` matching `propose`): `Damper`, `Force`, `Conductance`, `Heater`, `Admittance`, `FlowSource`, `Gyrator`.

The dual laws support common-effort networks with the same source clocks,
work, heat and entropy accounts as their ordinary proposal path. A publishing
flow source retains both endpoint values in its signal output. The gyrator
keeps its lossless cross-domain power relation.

`LoweringGraph::add_with_tangent_and_residual` registers components carrying both
`PortTangent` and `StageResidual`; the default lumped and acoustic cavity
constructors in `lowering.rs` use it. A cavity's integrated volume is its storage
unknown, so the existing cavity/duct graph can participate in monolithic solves
and stationary-analysis admission. Other elements such as `GapField`,
`ControlledSource`, `Slider`, `Restriction`, `GasVolume` and `ConductanceMatrix`
still lack this tier. Direct `SystemBuilder::add` and `add_with_tangent` callers
do not attach `StageResidual`, even when the concrete type implements it.

`tests/elements_stage_residual.rs` checks residual conformance, derivative
identities, current-driven capacitor transients and the independent equilibrium
`V = I/G`, and refusal when the tier was not registered. The acoustic lowering
test compares cavity/duct evolution with an independent midpoint recurrence,
the interface solver, and the energy account.

## Assembly source admission

The lowering registry checks source-history dimensions against the receiving
port before constructing a source or a scheduled boundary. Power ports require
their domain's canonical effort or flow units according to the input role;
heat ports require temperature in kelvin or heat rate in watts. Scalar histories
cannot broadcast into vector ports or mixed-unit material streams. Signal ports
currently declare no physical dimensions, so only dimensionless scalar source
histories are supported there. A source refusal names its declaration and target
port. An assembly declaring a controller refuses by name until controller
lowering exists; its feedback law is never silently dropped. These admission
checks leave the existing kernel signal/controller components available to
callers constructing their graphs explicitly. Regression coverage uses the
lowering source tests and the real mechanical, thermal and electrical graphs in
`tests/lowering_registry.rs`.

A declared query-point key takes the axis spec that every claim bounding that
axis gives it: a dimension-only spec when the dims match, or absolute
temperature for a kelvin value. A key the claims type differently, or whose dims
differ from their spec, is refused by name. Under any other semantic kind the
key stays untyped, and the query refuses that axis.

## Owner groups and owner constructors (bead `frankensim-rc-1109-root-3gc12.11.2`)

`LoweringRegistry::register_group` adds a named group. Before any per-item
constructor runs, each group's selector picks parts, interfaces and sources of
the draft, and its constructor lowers each selection whole. Selected items skip
their per-item constructors, and their receipts name the group. Registering a
name again replaces that group in its place. An item two groups select is
refused, naming both groups.

With `slices`, `lowering::owners::register_owner_defaults` adds three groups.
They lower into the physics owners' own Components through those owners' public
constructors:

- `owners::mbd`: rigid bodies joined by `viscous-bushing` interfaces or
  card-less `ShearFilm`s lower to one fs-mbd `ViscousNetwork`. A film's normal
  is the bushing's z axis. Its damping is `mu A / h` on x and y, `mu J / h` on
  rz with `J = A^2 / (2 pi)` (a full disk), and 0 on z, rx and ry, each
  receipted with its reason. `mu` is read from the third body's card at the
  declared query point and held there (one-way coupling). A film's receipt
  has a null `card_hash`, family `shear-film`, and its fluid as
  `third_body_card`. Bodies joined by `relaxation-bushing` interfaces lower
  to one `RelaxationNetwork` instead, and a network that mixes the families
  is refused, naming both sets. Such a card holds the damping and, for each
  branch `k`, a diagonal stiffness `stiffness_<k>_x` .. `_rz` (factor
  `L = sqrt(K)`) and `relaxation_time_<k>`. The memory has no default: the
  interface's state declares each branch's elastic displacement
  `displacement_<k>_x` .. `_rz`, held as `a = L x`.
- `owners::thermal`: thermal solids whose material card names a registered
  caloric law lower to the cells of one fs-thermal `ThermalTimeline`. The
  card picks the path: a model-bearing card on a thermal-only `Solid` moves
  that part from the kernel's lumped lowering to the timeline. A
  `thermal-boundary` interface from a cell's face to a face of a stationary
  `DrivenBody` frame is the cell's exchange with the ambient. A cell outside
  the product law's `T` band is refused at lowering or stops the run; nothing
  clamps. The time plan's `scheme` is the timeline's; the mechanical owners
  run beside any scheme. Under `backward-euler` or `sdirk2` the cells step on
  fs-thermal's implicit conduction graph, with the step, tolerance and work
  limits above and the graph's other controls at their defaults. Each law
  gives its implicit form; the product law's needs its `T` band and the
  card's `thermal_conductivity`, read at the declared query point and held
  for the run (`thermal_conductivity.held`). A law without one, and a
  `thermal-boundary`, are refused by name.
- Loss routes: a card-less `HeatRoute` interface (bookkeeping, not a material
  interface) from a part that exactly one bushing joins to a cell carries its
  declared fraction of that bushing's heat. The rest goes to a sink named
  `unallocated`. Its receipt's `card_hash` is null, its family `heat-route`.

Interface sides are surfaces or regions, so each of these joins faces; a
route's bushing is found by part, since a `.fsim` gives each interface its own
faces. Owner groups lower before any part, so the frame is recognised by its
declared kind, not by the kernel's ground mark.

The groups override no kernel default, so a draft naming none of these
families or laws lowers as before. A body's pose and motion come from the
part's declared initial state (`DraftPartDecl.initial`), and a body without one
is refused. A cell's temperature comes from the initial state, or the ambient
when none is declared. fs-cli lowers every assembly through the kernel
defaults plus these groups.

fs-couple's `rigid_damping_heat` target keeps the hand-built graph as the
oracle and lowers the same pairing declared as a draft. The lowered network,
cells and routes match that graph bit for bit in every step size, body state,
cell temperature and ledger account. The lowered ledger closes to 1e-12 of the
dissipated energy. The bodies' mass properties are held to the boxes' own.
The beds carry the fixture's reservoir, which exchanges nothing, so a nonzero
exchange is not compared. Replacing any one group with an empty one is refused
by name. Only the constant source is declared: the fixture's
contact and its 4 W event at 0.35 s are thermal events a draft cannot
declare, so both sides of the compare leave them out. Fixture values; not a
seed-dataset authority.

## Live passive constitutive laws for scalar dissipation

A Damper can follow an `Arc<dyn DissipationLaw>` through `Damper::following`.
The law receives the LIVE heat-port temperature and signed relative flow
f_a - f_b at every proposal and tangent evaluation. It returns a finite
nonnegative dimensionless factor and its two finite partials, ordered
[temperature, relative flow], multiplying the original reference resistance.
The resulting effort is R(T,u) u and the heat output is -R(T,u) u^2.
Both JVP and VJP carry the dependence through force and emitted heat.

The law may explicitly return zero dissipation. Negative or nonfinite
factors, unrepresentable responses, and declared source-support failures
refuse the trial; negative differential resistance alone is not treated
as negative dissipation. This remains a Smooth component and confers no
monotonicity or global uniqueness guarantee. A factor cusp at zero flow
is admissible only when its effort and heat products have the supplied
continuous derivatives there.

The complete authored law contributes identity words in addition to the
reference coefficient. Stateless checkpoint encoding is unchanged.
Dampers without a law retain their original arithmetic and identity.
The hook adds no material dependency to fs-system: a physics/material
owner supplies its source law, while the shared component owns the heat,
entropy, tangent and composition behavior.

## Diffuse-gray radiative exchange

`ConductanceMatrix::gray_radiation(area_m2, [emissivity_a, emissivity_b])`
constructs a two-port heat interface for equal-area, diffuse-gray surfaces
with full mutual view. It uses the existing conductance component's proposal,
heat delivery, entropy receipt, tangent and identity paths; it adds no
integration driver or separate surface state.

The heat rate from `a` to `b` is
`sigma * area * epsilon_eff * (T_a^4 - T_b^4)`, where
`sigma = 5.670374419e-8 W/(m^2 K^4)` and
`epsilon_eff = 1 / (1/epsilon_a + 1/epsilon_b - 1)`.
The effective emissivity is evaluated without taking reciprocals of tiny
emissivities. The temperature factor uses the fourth-power secant, preserving
small temperature differences. A 1 K reference normalizes the law's factor
to a dimensionless value and its base conductance to W/K; its two exact
partials feed the ordinary JVP/VJP. Every proposal uses the current port
temperatures. The optional `ConductanceLaw::rate_partials` supplies direct
normalized-rate derivatives; other laws retain the product-rule fallback.
Radiation uses direct cubic partials so a large temperature contrast does not
cancel the colder surface's derivative.

Admission requires finite positive area and emissivities in `(0, 1]`.
Unrepresentable coefficients or law evaluations refuse; the interface does
not clamp temperatures. Identity retains the original area and ordered
emissivities. Fixed authored properties confer no material-card provenance.
There is no participating medium, general enclosure/view-factor solve,
coating/half-cell resistance or separate surface storage.
`tests/radiation.rs` contains independent Stefan-Boltzmann, reversal,
equilibrium, tangent and finite-body energy-balance checks.

## Recovery and residual decrease in smooth Newton blocks

For a feedback block containing a `Smooth` component and otherwise only
`Smooth` or `ExactLinear` components, the driver retracts a generated
Newton update when an owner refuses it or its successful proposals fail
a sufficient residual-decrease test. It retains the preceding valid inputs
and halves the same correction. Each retry consumes the existing
`max_iterations` budget; rounding that returns every input to its starting
value stops the search.

The search merit is the maximum residual above a frozen absolute/roundoff
allowance, normalized by scales from the preceding valid iterate and the
full Newton correction. These remain fixed on that correction's search,
so a growing trial cannot disguise its error with its own large outputs.
The accepted trial satisfies merit <= (1 - 1e-4 * fraction) * prior_merit.
A coordinate with zero prior scale and zero correction is omitted only
from this merit: curvature can open it, and the ordinary residual checks
still judge it immediately.

If that search exhausts after any typed owner refusal, its first such
refusal is returned. If all proposals succeeded, `NotConverged` reports
the actual last trial residual. This does not publish candidate state.
An exact fixed point still requires no iteration budget.

Recovery starts only after valid initial proposals. Initial-input
refusals, cancellation, tangent failures, malformed shortening requests and
`Proposal::Shorten` retain their meanings. Load-line blocks and scalar
bracket plans retain their policies. Residual propagation, conservation,
admission and atomic publication still govern physical-step acceptance.
Sufficient decrease does not establish global convergence or uniqueness.

Numeric workspace is allocated at sealing; the merit reuses existing
scale/floor scratch. A recovered typed error still uses the existing boxed
error representation, so that failure path is not claimed allocation-free.
`tests/smooth_support.rs` compares logarithmic and exponential feedback
with independent analytic/bracketed references, with declared/differenced
tangents, finite overshoot, initially zero nonlinear coordinates, budget
exhaustion, rollback/replay, cancellation and shortening.

## State-triggered thresholds

`System::step_to_event(dt_max, config, guard, cx)` locates a directional zero
of a caller-supplied continuous guard on accepted physical state or time.
The initial ordinary step respects scheduled events and owner step limits;
only its actual accepted interval is searched. Bisection replays that same
kernel from the entry checkpoint. The driver retains one ordinary audited
step at the crossing-side endpoint once its duration bracket is within
`EventConfig::time_tolerance_s`, or at a sampled exact zero. Its report
includes the bracket and guard values, retained `StepReport`, and trial count.
Discarded trials contribute no accepted steps or ledger rows. The state,
compensated clock, warm starts, port values, configuration, ledger, last-step
receipts and interface allowances are restored together through the public
checkpoint. Cancellation, nonfinite guards, owner refusals,
spent bisections and unresolvable brackets restore the complete entry state.

An exact zero on entry returns `AtStart` without advancing, independently of
direction; the caller must handle or disarm it. No bracket returns `Advanced`
with the original ordinary step. A shortened refinement is sampled at its
actual duration and refuses if it cannot advance the lower bracket endpoint.
The guard, boundaries and proposals must be externally pure, except for
cancellation requests. A guard may not depend on last-step means, counters or
solver diagnostics. This assumes at most one zero per trial interval; hidden
multiple crossings, grazing roots and discontinuous guards are not detected
reliably. It locates the numerical one-step trajectory, supplies no physical
reset/impact law, and makes no trajectory-accuracy or allocation-free claim.
The initial full ordinary trial must be admissible. An owner that refuses an
invalid endpoint causes this call to refuse; the driver does not first recover
and locate an event preceding that endpoint.
`tests/events.rs` covers ballistic and heated-body time oracles, scheduled
caps, direction, replay, retained accounting and diagnostic allowances, and
atomic cancellation, refusal, trial-budget and clock-resolution boundaries.

`System::step_to_events(dt_max, config, directions, guards, cx)` extends the
same driver to a nonempty set of competing thresholds. The callback fills one
finite output for each direction; unwritten outputs are reset to NaN and
refused with their index. Every output is checked before returning all exact
entry zeros as `EventSetStatus::AtStart`, independently of direction. Each
guard must satisfy the same continuity and at-most-one-zero assumptions.

All guards share one duration bisection and `EventSetConfig` trial budget.
Its lower endpoint precedes every requested crossing; its upper endpoint
reaches at least one. Once the common bracket is within tolerance,
`EventSetStatus::Located` contains every requested crossing reached in that
bracket, in declaration order. These events are numerically unresolved within
the bracket, without a claim of exact simultaneity. A later root just beyond
the retained endpoint is not included, even if its distance from the earliest
root is less than the tolerance. An exact sampled upper zero collapses the
bracket only if it is the sole crossed guard: an exact later zero must not
hide another guard's earlier crossing. Final guard values and all located
brackets correspond to the single retained ordinary step. Increasing the
guard count adds guard work, not a separate physics search per guard.

The scalar API uses this same search and preserves its original report and
`EventError` variants; the batch API returns `EventSetError` with indexed
invalid outputs, an explicit empty-set refusal, and every unresolved crossing
on budget exhaustion. Both have identical complete rollback and cancellation
semantics. No reset, priority-dependent reset ordering, or impact law is
introduced. `tests/event_sets.rs` covers competing physical thresholds,
declaration order, exact later roots, unresolved simultaneous brackets,
initial-zero validation, scalar parity, shared physics work and full rollback.

## Scalar monotone feedback

A `Class::Monotone1d` component without `PortTangent` can close a feedback loop
when it has one unknown input and removing that input's incoming dependencies
makes the rest of the block acyclic. Sealing retains that evaluation order.
Each trial evaluates the original components and connections in that order;
an expanding, safeguarded scalar bracket locates the feedback solution.
Refused exploratory inputs narrow the search domain without inventing a
residual sign or publishing a candidate. The initial owner's refusal or
shortening request remains authoritative, and cancellation always propagates.

The resulting candidate still passes the ordinary per-wire and propagated
residual checks, independent conservation audits, and atomic admission.
Bracket probes and subsequent Newton corrections share `max_iterations`.
The path allocates its plan at sealing and reuses the original step buffers.
The native nonlinear-storage regression compares with an independently solved
monolithic cubic and checks energy, finite support, declaration order,
shortening, cancellation, exhausted work and exact checkpoint replay.
Multiple remaining feedback coordinates and general set-valued relations do
not acquire a scalar solve; they retain their existing tier requirements.
This supplies a numerical solve for an admitted continuous relation, without
a global uniqueness, global error or physical model-accuracy certificate.

## Purpose and layer

Layer: L3. The composition kernel: a system of `fs_ports::component::Component`s
joined by connections and stepped as one. It owns what each coupled pairing
used to write for itself: the state store with atomic publication, the
equations of the connections, the order of evaluation, the interface solve,
step negotiation, and the energy, measured mass and entropy ledger.

Direct dependencies: fs-ports (the contract) and fs-exec (`Cx`). It sits
outside fs-couple so that the kernel compiles and tests in seconds and so that
physics owners below fs-couple can be composed without it.

## Public types and semantics

- `SystemBuilder`: `add` (mandatory tier), `add_with_tangent` (component offers
  `PortTangent`), `add_on_line` (component offers `LoadLine`: a relation that
  is not a function of the flow, see "The load line" below), `port`, the
  connections below, `boundary`, `seal`.
- Connections, all generated by the kernel from data:
  - `common_effort(setter, members)`: the setter returns the effort every
    member receives and receives minus the sum of the members' flows.
  - `common_flow(setter, members)`: the setter returns the flow, every member
    receives minus that flow, and the setter receives the sum of the members'
    efforts.
  - `bond(a, b)`: the two-port case of either.
  - `transformer(a, b, ratio)`: `a` receives effort and returns flow; `b`
    receives flow and returns effort. The equations are `f_b = -ratio f_a`
    and `e_a = ratio e_b`. A signed finite nonzero ratio converts a lever,
    gear, piston area or electrical turns ratio, with equal port widths.
    Unlike `bond`, explicit transformations may cross power domains; the
    ratio supplies the conversion in the ports' SI units. Two inertias that
    both receive effort still require a constraint solve or reflected-inertia
    lowering, not this connection.
  - `common_flow_ratios(setter, [(member, ratio)])`: a common-flow junction
    whose members each sit behind their own ideal transformer. The setter
    returns the flow `f`; member `i` receives `-ratio_i f` and the setter
    receives `sum(ratio_i e_i)`. One member is `transformer`; all ratios
    equal to one is `common_flow` between power ports. Members may be in
    different power domains. This is how one inertia carries several loads
    at different ratios: a wall of area `A` between two gases under a force
    is the wall's mass with one gas at `A`, the other at `-A`, and the force
    at one. `common_effort_ratios` is the dual: member `i` receives
    `ratio_i e` and the setter receives `-sum(ratio_i f_i)`. Both need power
    ports of equal width and finite nonzero ratios.
  - `power_map(a, b, P)`: the same roles, with `P` a finite row-major matrix
    of `b.width()` rows and `a.width()` columns. The kernel generates
    `f_b = -P f_a` and `e_a = M_a^-1 P^T M_b e_b`; rectangular and
    rank-deficient maps are allowed. A `PortSpec::field_power` declares a
    fixed dense symmetric positive-definite metric for `e^T M f`; ordinary
    ports use the identity. The builder computes the effort adjoint with
    the metrics' admitted Cholesky factors; stepping performs no metric
    solve. Non-finite adjoint arithmetic refuses as `MapAdjoint` before any
    port is claimed. With two ordinary ports the original shared-coefficient
    transpose path is retained. Heat, stream and signal ports are refused.
    Bonds, ordinary junctions, ratio junctions and scalar transformers need
    equal power metrics; an explicit identity matches an ordinary port.
    Unequal metrics require `power_map`. Interface power magnitudes and
    input-error energy allowances use the declared metric, and a non-finite
    interface energy allowance refuses the step.
  - `power_map_members([(a_k, P_k)], b)`: several effort-receiving ports onto
    one flow-receiving port as one connection. The kernel generates
    `f_b = -sum(P_k f_a_k)` and, member by member, `e_a_k = M_a_k^-1 P_k^T
    M_b e_b`: the map `[P_1 ... P_n]` of the members' coefficients taken
    together, such as the lever from two bodies' force and torque ports onto
    one contact. `power_map` is its one-member call, with the same drivers
    term for term; `common_effort_ratios` is its scalar case, `P_k = r_k I`.
    Each member is admitted as `power_map` admits its port, before any port
    is claimed; no member is refused as `Empty`. Limit: the map is constant
    from admission. A lever whose arms turn with its bodies keeps the ratios
    of the pose it was admitted at, so a graph that moves its bodies through
    a step takes their arms at the admitted pose.
  - `heat_route(source, sinks, remainder, temperature)`: heat leaving the
    source is delivered to the sinks by fraction; a remainder goes to a named
    `HeatSink` outside the system. `RouteTemperature::SinkEquivalent` gives the
    source `sum(f) / sum(f / T)` of its receivers, so the route generates no
    entropy; `RouteTemperature::Fixed` is for a source whose law ignores
    temperature and makes the route one-way. A receiver may be named by
    several routes: it receives the sum of their shares. Routes that share a
    receiver are audited together, in the ledger row of the first of them;
    the others' rows stay zero. A source belongs to one route, and a
    receiver of a route takes no other connection.
  - `signal(from, to)`: one-way copy.
  - `signal_map(from, to, map)`: a one-way dense linear signal map, with
    `to.width()` rows and `from.width()` columns in row-major order. This
    selects individual outputs of a vector controller for scalar actuators,
    or mixes signal coordinates. Negative and zero coefficients, zero rows,
    rectangular and rank-deficient maps are admitted; nonfinite coefficients,
    wrong shapes/roles and an already-driven input refuse before changing the
    builder. One output may feed several mapped or ordinary connections.
    Exact-zero coefficients add no dependency. The fixed map participates in
    the ordinary execution plan and loop tangent, and in durable problem
    identity. It contributes no energy, storage, delay or physical transformer
    law; coordinate units and gains are authored by the caller.
  - Streams (`Carrier::Stream`, crates/fs-ports/CONTRACT.md) need no
    connection of their own. A storage's stream port is the setter of a
    `common_effort` (or one end of a `bond`): every transport element joined
    to it receives what the storage holds is like, and the storage receives
    minus the sum of what crosses into them. Mass, enthalpy and each carried
    quantity therefore cancel across the junction by the same linear
    equations as a power junction's flows. A `common_flow` of stream ports
    with more than one member is refused (`StreamCommonFlow`): one flow
    cannot be handed whole to several storages. Two storages, or two
    transport elements, do not join directly. Ratios, transformers, power
    maps and heat routes refuse stream ports. Each position of a stream port
    is a tolerance group of its own in the interface solve, because its
    coefficients are different quantities. Stream connections are audited
    separately for mass cancellation as well as enthalpy balance.
  - Field heat ports (`PortSpec::heat_field`, crates/fs-ports/CONTRACT.md)
    join through `bond` or `common_effort` to ports of the same width, point
    by point, each in its owner's order. A heat route refuses one
    (`HeatRouteWidth`, the port and its width named): a route moves one heat
    rate at one temperature. In a connection's energy audit a field heat
    port counts `dt sum |Q_i|` toward the magnitude the tolerance is
    relative to, not `|dt sum Q_i|`: a trace's loads may cancel in their
    sum. The solve's slack and its entropy share are taken per point, as for
    every input. Two traces meet through a `ConductanceMatrix`.
    Identified stream ports additionally require the same medium/source laws,
    enthalpy reference, and ordered quantity identities and dimensions. A
    mismatch (including an unidentified peer) refuses as `StreamBasis` before
    claiming any port. An explicit adapter is needed to change basis or datum.
    `GasVolume::identified` binds its gas-constant and heat-capacity-ratio bits,
    its absolute-zero energy datum, and the caller's named carried masses;
    `Restriction::for_stream` retains a source port's exact basis on both ends.
    `Restriction::for_gas` takes its ports from a `GasVolume`'s stream port the
    same way and keeps the gas, whose state function gives the entropy the
    restriction reports.
    Anonymous constructors remain available with their original caller
    obligations. Basis declarations authenticate no physical source law, and
    the ledger's ordinal `carried_residual` still requires corresponding
    meanings across all devices being summed, including disconnected devices.
    The identified/unidentified flow-reversal runs in `carried_accounting.rs`
    agree bitwise in trajectory and have equal ledgers; admission regressions
    cover species, order, units, material and energy-reference mismatch.
- `Component::stored_mass_kg(state)` optionally measures total stored mass
  [kg] from the typed state. A measured inventory is finite and nonnegative;
  `Some(0.0)` declares zero storage, while the default `None` means unknown.
  For measured owners the driver checks the change between accepted and
  candidate inventories against signed mass entering through stream ports.
  This check is independent of the energy receipt and uses the mass
  contract's rounding allowance, not a component-supplied mass tolerance.
  `elements::GasVolume` validates its contents and reports `mass_kg`;
  `elements::Restriction` reports zero. Other owners must opt in before the
  kernel can claim their inventory balance.
- `Component::stored_entropy_j_per_k(state)` (fs-ports) optionally measures
  stored entropy [J/K] from the typed state; the default `None` means not
  declared. The driver asks every owner, at the accepted and the candidate
  state, and when both are declared checks with fs-ports'
  `audit_entropy_step` that the entropy entering through heat ports
  (`dt q / T`) plus the receipt's `entropy_generated_j_per_k` equals the
  measured change, within the fs-ports rounding allowance. The allowance is
  relative to the reported inventories, so an owner measures its entropy
  with rounding relative to the value it reports: `C ln(T / T_ref)` taken as
  `ln` of the quotient is refused near `T_ref` (`tests/entropy_accounting.rs`).
  `elements::Conductance` declares zero: it stores nothing, so each step's
  generation must cancel the entropy its heat carries through its ports.
  `HeatCapacity` does not declare: it books no generation and publishes its
  step-mean temperature, so its measured change would differ from its port
  entropy at third order in the step.
- `boundary(port, f)`: the input side of a port prescribed as `f(t0, dt)`.
- `scheduled_boundary(port, f, event_times_s)`: the same input with immutable,
  finite, strictly increasing absolute event times. The step negotiation
  respects the next future event using the kernel's compensated clock.
- `System::step(dt_max, cx) -> StepReport`, `state`, `last_input`,
  `last_output`, `last_allowance`, `guess`, `plan_shape`, `plan_loops`,
  `time_s`, `steps`, `ledger`, `set_config`. `plan_loops` lists, per loop
  of the plan, how many inputs it solves for and how many of those are
  solved together; the others are copies and follow by substitution.
  `last_allowance(port)` is what the solve of the last accepted step told
  the interface audit it had left in each input of the port, in the
  input's unit: rounding and the configured absolute tolerance, never the
  relative one (for a block with a load line see "The load line"). A
  checkpoint retains it together with `last_receipts`, so both describe the
  restored last accepted step rather than a discarded trial.
- `System::step_toward(end_s, dt_max_s, cx) -> StepReport` caps one ordinary
  step at an absolute endpoint using the compensated clock. It retains the
  clock's low bits when finding the remaining interval; when a full cap would
  leave only clock roundoff before the endpoint, it takes at most half the
  remaining interval. Component limits, scheduled events, trial shortening,
  audits and atomic publication still use `step`. The endpoint must be finite
  and strictly later than `time_s()`, and the cap finite and positive. Repeated
  calls support nonmultiple horizons without a rounded-time subtraction tail.
- `System::step_adaptive(dt_max, config, metric, cx) -> AdaptiveReport`:
  error-controlled stepping of the same coupled graph. From one checkpoint,
  compare one coarse step to two equal half steps at the same end time. The
  caller supplies a norm over typed states or account increments and declares
  the least time order of the coupled method. `ErrorNorm` combines explicit
  `ErrorTolerance { absolute, relative }` per coordinate; invalid values or
  an empty norm refuse instead of disappearing in a maximum. Divide the
  weighted coarse/refined difference by `2^order - 1`, accept at most one
  tolerance unit, and retain the two audited half steps without extrapolation.
  A shortened coarse step sets the comparison interval; a shortened half
  discards the whole attempt and restarts from the original checkpoint.
  Component refusals, cancellation, invalid estimates, exhausted attempts,
  and the minimum-step/clock-resolution floor restore every accepted state,
  warm start, clock and ledger. Discarded trials add no accepted ledger rows.
  The report counts all completed trial solves and suggests the next step
  using a bounded proportional controller. `System::steps` counts the two
  accepted substeps. This is an estimated local error target, not a global
  bound: order, covered accounts, source smoothness and event resolution
  remain caller obligations. Source-work norms should compare interval
  increments, so a large past account cannot hide a new error. External
  side effects of boundary, metric or component calls cannot be rolled back.
- `System::checkpoint() -> Result<Checkpoint, CheckpointError>` captures
  accepted typed states, interface warm starts and outputs, last-step receipts
  and interface allowances, the exact compensated clock, step count, complete
  ledger and solve configuration.
  Mass and entropy accounts, coverage counts and both words of every account
  are part of that ledger and rewind with the other accepted accounts.
  `restore(&checkpoint)` atomically rewinds the originating sealed system;
  the same checkpoint may be reused. A checkpoint from another system is
  refused before anything changes, even when the graphs have the same names
  and types. `Checkpoint` exposes read-only `state`, `time_s`, `steps`,
  `ledger` and `config`; it keeps no trial candidates. Capture after a failed
  or cancelled step still describes the last accepted state. Boundary
  closures and their captured external data are not copied: deterministic
  replay requires identical boundary responses and pure component proposals.
- Durable snapshots (bead `.5.4`, item 4):
  `System::snapshot(&context, limits, cancellation) -> Result<SealedSnapshotV2, CheckpointError>`
  seals the accepted generation in fs-exec's v2 envelope.
  `System::resume(&opened)` makes an `OpenedSnapshotV2<SystemImage>` the
  accepted generation of a system, in the process that sealed it or in another.
  - `SystemImage` is the envelope's state type. It is one state family of this
    crate; its charter's identities are pinned and checked by every seal.
  - Per member, the image holds the identities of the member's codec charter
    and the plain words of its state (`Component::encode_state`). It also
    holds the accepted wire values, the compensated clock, the step count, the
    solve configuration and every number of the ledger. Names are not in it;
    they come from the resuming system.
  - The image's charter is at codec version 3. Since version 2 (item 5) each
    ledger account is written as the hi and lo words of its normalized
    double-double total, and the absolute accounts sit beside the signed ones.
    A decode refuses a pair that is not normalized. The carried accounts
    (06a4a2d4ba) and the per-node entropy flag and account (17c5438df5) were
    added to the words under version 2's number, so version 3 has version 2's
    schema and another codec identity. fs-exec refuses an image of version 2
    as `WrongStateSchema` before reading its words. An image of codec
    version 1 (7a35133254, a sum and its compensation per account) names the
    same state type and another schema and codec, so fs-exec refuses to open
    it as a `SystemImage`.
  - The context names the problem. `System::problem_id()` hashes:
    - the `Debug` text of the description the system was sealed from: names,
      ports, classes and tiers, load lines, connections, driven ports and heat
      sinks;
    - each member's charter identities;
    - the words of each member's `Component::identity_words`, its parameters.

    Every float in a description is admitted finite and printed exactly by
    the derived `Debug`, so the text tells descriptions apart. A toolchain
    that prints one differently gives another identity: it refuses older
    snapshots and never makes two descriptions one. Boundary functions are
    not in the identity, so a resumed run follows the same trajectory only
    under the same boundary responses.
  - `snapshot` refuses a context of another problem (`ForeignProblem`) and a
    member without a codec (`NoCodec`).
  - `resume` decodes every member and the ledger before it changes anything.
    It refuses:
    - a snapshot of another problem (`ForeignProblem`): another graph, another
      member charter, or a parameter one ulp away;
    - a member that refuses its words (`Decode`);
    - an image of another shape (`Image`).
  - One opened snapshot resumes into several systems, which then go on alone
    (a fork). It also rewinds the system that wrote it.
  - As after `restore`, the allowance of the last step is zero until a step
    is taken.
  - Opening belongs to fs-exec. In the process that sealed, use the root it
    kept (`unseal_v2_expected`). In another process, use an admitted
    authority over the envelope's subject (`unseal_v2_authorized`), with the
    caller's verifier and policy.
  - `fs-blake3` is a normal `[dependencies]` entry of this crate, not a
    dev-dependency: `problem_id` hashes with it.
  - Codecs of the kernel's elements:
    - `GasVolume`: `3 + carried` words. A decode runs the volume's contents
      check. Identity: gas constant, ratio, carried count.
    - `Restriction`: no words. Identity: conductance, carried count, and its
      gas's identity if it has one.
    - `Mass`: two words, position and velocity; identity: inertance.
    - `GapField`: two words, the field coordinate and the gap; a decode
      refuses a gap outside the open range. Identity: the kind (capacitor
      or inductor), the coupling scale and both limits.
    - `Storage`, `Spring`, `HeatCapacity`: one word, respectively integrated
      flow, extension and temperature. Identities retain capacity, both spring
      coefficients, or thermal capacity.
    - `Force`, `FlowSource`, `Heater`: three words, both source-clock words
      and signed energy received. Identity retains the waveform kind and every
      exact constant, sine parameter or ordered table knot.
    - `Damper`, `Admittance`, `Conductance`, `Gyrator`: no words. Identity
      retains the original damping, conductance or conversion coefficient.
    - `HeatNode`, `Resistance`, `DiracJunction`, `ModulatedTransformer`: no
      words. Identities retain all conductances, the complete resistance or
      incidence matrix, and the Dirac schema's coordinate and clock bindings.
    - `ConductanceMatrix`: no words. Identity retains the number of interface
      points and, of each, both sides' points and weights and its conductance,
      and then, if it follows a law, the law's words, the number of guard
      rows and each row's points and weights.
      Signal-modulated ratios belong to the restored inputs and outputs.
    - `ControlledSource`: the same three-word source account, with power
      domain and effort/flow role bound by the graph's port description.
    - `Slider`: eight words retain slip, stick state, both prescribed-end
      histories and their mirrored-span counters. Invalid option/boolean tags,
      nonfinite values and impossible counters refuse. The accepted counter
      that causes a subsequent span refusal survives restart unchanged.
      Identity retains both friction force limits and prescribed-end layout.
    - With `constitutive`, `NonlinearCapacitor` and `NonlinearInductor` each
      store their one exact charge/flux coordinate. Decode reuses the original
      source law's support and energy checks. Identity includes every knot of
      the original curve, without resampling or fitting.

    Native codecs preserve floating-point bits, including signed zero and
    source-clock compensation; they reject extra or missing words, nonfinite
    state and unrepresentable stored energy. A heat capacity additionally
    requires positive absolute temperature. A failed encoding appends nothing.
    Power domains and optional signal ports remain in the graph identity, so
    changing an element's domain or publishing layout also refuses resume.
    Other components still require their own codecs before a containing system
    can be sealed. `ledger::Sum::to_words` and `from_words` give a sum's two
    words exactly, for the codec of a state that holds one.
- `Ledger`: per component (port energy in, stored change, dissipation,
  numerical loss, entropy generated, accepted balance defect and its absolute
  sum), per connection (energy defect, entropy generated, and the absolute
  sums of its energy and mass defects), per boundary (energy and entropy in),
  per heat sink (heat); `energy_residual_j()` is the first-law residual of the
  window.
  - Every account is a `ledger::Account`: a double-double total
    (`fs_math::Dd`; `fs-math` is a normal `[dependencies]` entry of this
    crate). `value()` is the total rounded to f64, `dd()` both words.
  - The readers that combine accounts (`energy_residual_j`,
    `boundary_energy_in_j`, `stored_change_j`, `sink_heat_j`,
    `numerical_loss_j`, `connection_defect_j`, `entropy_generated_j_per_k`,
    `mass_residual_kg`, `carried_residual`) add the accounts' double-doubles
    and round once at the end. A residual far below the rounding of its
    accounts survives: in `tests/ledger_accounts.rs` two accounts of 1e6 J
    leave 2^-40 J, which the reader gives exactly and rounding each account
    first gives as 0.
  - The clock, the source clocks and `Ledger::duration_s` stay compensated
    `ledger::Sum`s, and the step audit is unchanged: it combines one step's
    port values, not accounts, and decides acceptance. No time, state or port
    value moved with the change to double-double (item 5's line compare of
    the kernel users).
  - The absolute accounts (`NodeAccount::absolute_defect_j`,
    `MassAccount::absolute_defect_kg`, `CarriedAccount::absolute_defect`,
    `EntropyAccount::absolute_defect_j_per_k`,
    `ConnectionAccount::absolute_energy_defect_j` and
    `absolute_mass_defect_kg`) add the size of each accepted step's defect
    where the signed account adds the defect, at publish only. A refused step
    adds to neither. Defects of opposite sign cancel in the signed account and
    not in the absolute one.
  - `largest_absolute_energy_defect_j()` is the largest absolute energy
    account of a component or connection. It is a report: `System` refuses
    nothing on it or on any absolute account. The hand-coupled
    `rigid_damping_heat` of fs-couple refuses on its own absolute sums against
    `energy_tolerance_j` at the end of each `advance_with_source_driver` call
    (`crates/fs-couple/src/legacy/rigid_damping_heat/mod.rs`, lines 554 and
    715-728). Whether the kernel should refuse the same way is a question for
    the owner.
  `NodeAccount::has_streams` records whether a node participates in stream
  transport, and `NodeAccount::mass: Option<MassAccount>` carries measured
  inventory accounting. `MassAccount` has `port_mass_in_kg`,
  `stored_change_kg`, `balance_defect_kg` and `absolute_defect_kg` accounts,
  plus `audited_steps`. `NodeAccount::entropy: Option<EntropyAccount>` is
  the same for a declaring owner's entropy: `port_entropy_in_j_per_k`,
  `entropy_generated_j_per_k` (the audited steps' receipts),
  `stored_change_j_per_k`, `balance_defect_j_per_k` and
  `absolute_defect_j_per_k`, plus `audited_steps`; `None` for an owner that
  declares nothing.
  `BoundaryAccount::mass_in_kg` records signed mass entering through a stream
  boundary, and `ConnectionAccount::mass_defect_kg` records the accepted
  uncancelled mass at a connection; non-stream rows contribute zero.
  `mass_residual_kg() -> Option<f64>` returns stored mass change minus boundary
  mass in minus connection mass defect when every stream owner has measured
  coverage for every accepted step in the window. Missing or incomplete
  coverage returns `None`, rather than a zero conservation residual. The
  known boundary and connection accounts remain available in that case.
- `SystemError`, `AdmissionError`, `OwnerError` (the component's own error,
  recoverable with `downcast_ref`).
- With the `constitutive` feature, `constitutive::NonlinearCapacitor` and
  `NonlinearInductor` wrap the existing `fs-phs` charge/voltage and
  flux/current curves. Their scalar state is charge [C] or flux linkage [Wb].
  The capacitor's port receives current and returns voltage; the inductor's
  receives voltage and returns current. Both integrate the conserved
  coordinate from the held input, use the source curve's exact energy change
  and discrete conjugate mean, and expose its analytic endpoint derivative
  through port JVP/VJP. Authored knots are retained without resampling;
  extrapolation is refused with the original `CircuitError`. The component
  class is `Monotone1d`: at a coincident knot, its derivative is the owner's
  one-sided Newton representative. No unique derivative at that knot,
  hysteresis, thermal evolution or material provenance is inferred. Read the
  retained curve for endpoint voltage/current; the port output is a step mean.
- `elements`: ideal lumped elements, each a `Component` with `PortTangent`
  (the `Slider` excepted: dry friction, a `LoadLine`, see "The load line"),
  constant SI parameters and a closed energy account. `Mass` (port `p`
  receives the net force, returns the step-mean velocity; state `Motion`),
  `Spring` (ports `a`, `b`; linear plus optional cubic term; force is the
  discrete gradient of the potential), `Damper` (ports `a`, `b`, `heat`; all
  absorbed power leaves as heat at the temperature `heat` receives, with
  entropy `P dt / T`), `Force` (port `p`, a prescribed `Waveform`), and
  `HeatCapacity` (port `heat` receives heat rate, returns the step-mean
  temperature; `publishing()` adds the signal output `temperature`, the
  temperature at the start of the step, which carries no power and leaves
  state and receipt unchanged, so a law that reads it is held over the step:
  first-order feedback), `Conductance` (heat ports `a`, `b`, entropy
  `G (Ta - Tb)^2 / (Ta Tb) dt`, stored entropy declared zero), `HeatNode` (a thermal node of no capacity
  joined to `n` bodies by conductances, with a heat source at the node:
  port `source` receives the heat rate put into the node and returns its
  temperature `T_n = (P + sum G_k T_k) / sum G_k`; `body<k>` receives that
  body's temperature and returns `G_k (T_k - T_n)`; each branch generates
  `G_k (T_k - T_n)^2 / (T_k T_n) dt`; exactly linear; it refuses a source
  that leaves the node no positive temperature. Without a source it is the
  conductances in series; with one it is the plane a heat is generated on
  between bodies, and the heat divides by where the bodies put the node, not
  by a prescribed fraction; `frankensim system` offers it as the scene kind
  `heat_node`), `Heater` (Joule heater, power never
  negative). `Force`, `FlowSource` and `Heater` keep a `SourceAccount` whose clock is a
  compensated sum like the kernel's, so they read the time a boundary
  receives bit for bit; their stored energy is the energy that entered them
  (negative when they deliver), so a driven device closes its first law.
  `Waveform` is a constant, a sine, a piecewise-linear table or a
  piecewise-constant held table with
  validating constructors (parameters finite and at most 1e150 in size);
  `mean(t0, dt)` is the step average in closed form, over the interval the
  two floating-point ends span, and NaN at a time that is not finite, which
  `Force`, `FlowSource` and `Heater` refuse.
  These sources advertise every future table knot as their next event, so
  the coupled driver ends a step at each change of slope, including the
  transitions to and from held endpoint values. The interval retains the
  source clock's compensation. If a requested cap would leave only clock
  rounding before the next knot, the source halves the remaining interval
  instead of creating a vanishing follow-up step; the cap remains strict
  and genuinely distinct adjacent knots remain distinct events. A motion
  source publishing its span can therefore carry a stuck slider across a change of
  acceleration without an artificial catch-up force or friction heat.
  Constants, sine waves and one-knot tables have no scheduled event.
  `Waveform::held_table` takes each new value at its knot and holds it
  until the next knot; its first value also applies before the first knot.
  Its events start at the second knot, and its mean weights the constant
  intervals without interpolating across jumps. The interpolation law has
  its own identity tag, so the same knots as a ramp describe another problem.
  Raw boundary closures still require the caller to cap steps at their
  own input events; a span alone does not locate an interior change.
  `Mass`, `Spring`, `Damper` and `Force` have `in_domain(kind, ..)`
  constructors that give the same four laws in another power domain of
  `PortKind`: inertia, torsion spring, rotary damper and torque about a
  shaft; inductance, 1 / capacitance, resistance and source voltage in a
  series current loop. `Gyrator::new(kind_a, kind_b, ratio)` joins two
  domains without loss (`e_a = ratio f_b`, `e_b = -ratio f_a`): the ideal DC
  machine between a current loop and a rotor. Port kinds keep the domains
  apart: a junction across two of them is refused (`Mismatch`).
  Two-port signs: joined to masses by common flow, a port receives minus the
  mass's velocity, so port `a` on mass A and `b` on mass B extend at
  `v_B - v_A` and return the force on each end. A `HeatCapacity` asks for a
  half step (`Proposal::Shorten`) instead of crossing absolute zero.
  The dual elements support common-effort networks:
  - `Storage::in_domain(kind, capacity)` has a flow-input port `p` and an
    `f64` state `q`, with `q_new = q_old + dt * flow`,
    `effort = (q_old + q_new) / (2 * capacity)`, and stored energy
    `q^2 / (2 * capacity)`. Its default domain is mechanical compliance;
    in the electrical domain it is the capacitor setting a parallel node's
    voltage. Capacity is positive and finite.
  - `Admittance::in_domain(kind, conductance)` has effort-input ports
    `a`, `b` and a temperature-input `heat` port. It returns opposite
    flows `conductance * (effort_a - effort_b)`, releases all absorbed
    power as heat, and accounts for entropy `dissipation / temperature`.
    Conductance is positive and finite. Its default is mechanical mobility;
    electrical conductance is the reciprocal of a branch's resistance.
  - `FlowSource::in_domain(kind, waveform)` receives effort at `p` and
    returns the prescribed step-mean flow. Flow is positive into the source:
    a negative electrical waveform supplies positive current to the node.
    It shares the validated waveform and compensated clock of effort sources.
    `publishing()` adds the signal output `span` of width two: the
    waveform's value at the two ends of the interval the mean is taken
    over, `t0` and `t0 + dt`, in the waveform's own sign. A mean does not
    say where a flow ends its step; a `Slider` on a side the source moves
    is told through it.
  Each supplies analytic port JVPs and VJPs. `Storage` and `FlowSource`
  are exact-linear; `Admittance` is smooth because its heat output is
  quadratic. The `new` constructors use the mechanical translation domain.
  - `Resistance::new(kind, matrix)` is the memoryless resistor N-port: one
    flow-input port `u` of width `n` and a temperature-input `heat` port. It
    returns the efforts `R u` and releases the absorbed power as heat, with
    entropy `dissipation / temperature`. `R` is row-major `n x n`, every entry
    finite, symmetric (a signed zero equals its mirror), and positive definite:
    the constructor
    refuses `ElementError::Asymmetric` or `ElementError::NotPositiveDefinite`
    (the first Cholesky pivot that is not positive and finite). The power is
    computed as `|L^T u|^2` from the factor `R = L L^T`, so it is never
    negative; the
    work in, `dt u . R u`, equals it within `(n + 2)` roundoff units of the
    balance scale `sum |R_ij u_i u_j| + (sum sqrt(R_ii) |u_i|)^2`. A diagonal
    `R` returns each row the `Damper` effort `r u` exactly. Smooth, with
    analytic JVP and VJP. `tests/resistance.rs` checks conformance and
    tangents over a sweep, the efforts against a sum by hand, the heat against
    `u^T R u`, the diagonal case against `Damper`, the refusals, and a mass
    on a common-flow junction with it, whose velocity falls by the trapezoidal
    factor `(1 - r dt / 2m) / (1 + r dt / 2m)` per step while the heat it books
    equals the kinetic energy lost and warms the body on its heat port. The
    kernel test has one row.
  - `ConductanceMatrix::new(width_a, width_b, points)` is the linear
    conductance between two heat traces, the many-point `Conductance`: field
    heat ports `a` and `b` receive the temperatures at their points and
    return the heat loads [W] into the element. The interface is a list of
    `InterfacePoint { a, b, conductance_w_k }`, stored as given: at each,
    either trace's temperature is interpolated from up to three of its
    points (`(point, weight)` pairs), and the heat `G (T_a - T_b)` crosses
    there from `a` to `b` and goes back to each side's points by the same
    weights. The loads are `K T` with `K = D^T W D`, one row of `D` per
    interface point (`a`'s weights and minus `b`'s) and `W` the diagonal of
    the conductances; nothing dense is stored, and `matrix()` assembles `K`
    for tests. The constructor refuses, as `ElementError::Parameter`: a
    width of zero, no interface point, a weight's point outside its trace, a
    weight that is negative or not finite, a side's weights that do not sum
    to one within 8 roundoff units, and a conductance that is not positive
    and finite. The weights make each interpolated temperature a mean of
    temperatures above zero, so the entropy the element reports,
    `dt sum_q G_q (T_a - T_b)_q^2 / (T_a T_b)_q`, is never negative: no
    field is refused for it. The loads sum to zero within the
    `(2 d + n + 18)` roundoff units of `dt sum_q |G_q (T_a - T_b)_q|` the
    receipt claims, `d` the most interface points that load one trace point
    and `n` the points of both traces. It stores nothing and books no
    dissipation. Exact-linear; JVP and VJP are `K` itself. With one point on
    each side and one interface point it is `Conductance`: the same loads
    and the same entropy to the bit. A contact between two P1 traces is the
    quadrature of `(T_a - T_b) (v_a - v_b) / R''` over their common
    refinement, one interface point per quadrature point
    (`fs_conduction::component::ConductionBody::matching_contact` and
    `mapped_contact`). `tests/conductance_matrix.rs` checks conformance and
    tangents over 200 drawn cases of four interface points over three and
    two points at any temperatures, each load against `K T` by hand, `K`
    symmetric with rows that sum to zero; the one-point case against
    `Conductance` to the bit over 100 draws; 2,000 drawn fields on a
    matching patch with every point at its own temperature between 250 K
    and 1500 K, all generating entropy; the refusals; a heat route from a
    trace refused by name; and two traces of two lumps each, stepped 50
    times through the kernel, against the monolithic backward-Euler solve
    of all four within 1e-11 relative, with the ledger closed and entropy
    generated.
  - `ConductanceMatrix::following(law)` makes the same element's
    conductances follow the temperatures: every interface point's
    conductance is its constant `G` times `law.factor(T_a, T_b)` at the two
    temperatures interpolated there, asked in every proposal. It is the
    one element, not a second one: without a law nothing changes, to the
    bit, and a law that answers one with no slope gives the constant
    element's loads and receipt to the bit. `ConductanceLaw` is the hook:
    `factor(kelvin_a, kelvin_b)` returns a positive finite factor with its
    two partial derivatives [1/K], or a reason it has no answer;
    `identity_words` joins the element's identity after its interface
    points. At a knot of a piecewise law the partials are the inward
    segment's, the convention of fs-conduction's
    `linearize_temperature_for_newton`. Under a law the element is
    `Class::Smooth`. Its JVP applies, at each interface point, the block
    `G [f + (T_a - T_b) f_a, -f + (T_a - T_b) f_b]` on the two interpolated
    differentials, and its VJP that block's transpose along the same
    weights: assembled point by point, nothing iterative, not symmetric,
    and exact transposes of each other. A law with no answer, or a factor
    or partial that is not finite or a factor that is not positive, refuses
    the proposal with `ElementError::Law { reason }`; the element never
    clamps a temperature or extrapolates a law, and the kernel leaves the
    state where it was. The law is asked at the interface points: a trace
    point outside the law that reaches every interface point through a
    mean inside it is answered, unless a guard row reads it.
    `guarded_at(rows)` adds guard rows to an element that follows a law. A
    `GuardRow` is the two interpolation rows of an interface point with no
    conductance: the law is asked there too, in every proposal and before
    the interface points, and a refusal there is the same
    `ElementError::Law`. No heat crosses at a guard row, it is no part of
    the tangent, and inside the law a guarded element steps as the
    unguarded one to the bit. Rows on an element with no law, and a row
    whose weights are not an interface point's, are refused at
    construction. Which rows a contact needs is its builder's to say (for
    matching traces, each pair of trace points that face each other). The
    entropy rule and the receipt's tolerance are the constant element's,
    with the scaled conductance, which is positive.
    `tests/conductance_matrix.rs` checks, with a law
    authored for the tests (`1 + a (T_a / 300)^2 + b T_b / 300` inside a
    temperature window, no material's): conformance and the tangent over
    200 drawn cases, each load against the law by hand, and
    `<J d, bar> = <d, J^T bar>` within 64 roundoff units; the factor one
    against the constant element to the bit over 100 draws, and the law's
    and the guard rows' words in the identity; the refusals, through
    `propose` and through the kernel with the lumps' temperatures kept to
    the bit; a lump above the law that the one interface point averages
    away, stepped without guard rows and refused through the kernel with
    them, the lumps kept to the bit, and ten steps inside the law equal to
    the bit with and without rows; and two traces of two lumps each,
    stepped 50 times through the kernel, against a backward-Euler Newton
    solve of all four written out by hand, within 1e-9 relative. No claim:
    no law is supplied here (a card-backed one is fs-conduction's), no
    guard rows are chosen here, and a law that follows anything but the
    two temperatures at a point is not this hook.
  `tests/parallel_elements.rs` checks component conformance and tangents,
  time reversal of the storage element at each of its 400 points in the four
  domains (`R` the identity on the state, bounds derived from the rounding of
  `Storage::propose` and its stored energy; a reversal that negates the state
  is refused at the outputs),
  floating capacitor nodes with parallel resistor branches, lossless LC
  oscillation, and a driven parallel RLC circuit including source work,
  shared heat and the total energy ledger.

**One step.** (1) The step is the largest that every `limit_dt_s` and
`next_event_s` allow. (2) Blocks are evaluated in dependency order; a block
without a loop is proposed once; a looped block is solved by Newton iteration
on its interface inputs. (3) A component may answer `Proposal::Shorten`; the
whole trial restarts with the shorter step. Inside a looped block that holds
for a proposal at an iterate of the solve, not for the kernel's own perturbed
proposals, and with one exception for a component on a load line ("Requests
to shorten inside a loop" below). (4) Every candidate is audited with
`fs_ports::component::audit_step`, with stored energy re-measured from the
states, and any measured mass inventory is audited against the stream
transfer. (5) Each connection must cancel interface energy and stream mass
within their respective tolerances, and the step must not destroy entropy.
(6) Every `admit` is called. (7) All states are published by swapping the two
generations, then the clock and the ledger advance. Mass rows advance only
with this accepted publication.

**Plan.** `seal` gives every scalar input one driver, builds the dependency
graph between components, and orders its strongly connected blocks (iterative
Tarjan, index order).

**The load line.** Dry friction holds any force up to its limit at zero slip,
and a closed contact pushes with whatever keeps its gap closed. The effort
such a component returns is not a function of the flow it receives, so the
block solve cannot iterate on that flow: there is nothing to differentiate.
What the component's neighbours do with an effort is a function, and frozen at
an iterate it is affine, `flow = released - C effort`. A component added with
`add_on_line` (`fs_ports::component::LoadLine`) is stepped on that line: it
is given the released flows and the matrix `C` and returns the efforts its own
relation allows there.

- Each owner's line is found once per step, before anything is judged, from
  the block's starting values. All load-line owners' outputs are held at what
  the last step accepted while the other components are evaluated. The
  block's interface Jacobian excludes every load-line owner's columns and
  is factored once; a back-substitution per effort gives that owner's
  self-response with the remaining unknowns eliminated. `C` is minus that
  response. The other load-line owners' efforts stay independent unknowns,
  not part of `C`: the full Newton solve resolves their cross-coupling.
  Each line's drivers add only `C * own effort`, so its owner still receives
  the line defined by the unchanged `LoadLine` contract. The matrices are
  frozen for the step; a nonlinear neighbour contributes its slope at the
  step's start. The one-owner case keeps its original arithmetic order.
- What finding the line evaluated is not evaluated again in the solve's first
  pass. A component other than the load-line owners that declares its tangent,
  and whose inputs are bit for bit those finding the line proposed it at, is
  not proposed again in the first pass: its candidate, outputs and receipt
  are already those of that proposal. The first assembly takes its columns
  from finding the lines. Between the two only the load-line owners' inputs
  are written (the released flows, then `release`), the step length is the
  same, and an accepted state changes only when a step is published; a
  proposal is pure (`fs_ports`: equal arguments give bitwise-equal results)
  and a tangent depends on the inputs, the accepted state and the step only.
  A differenced component is proposed again: differencing writes its
  candidate. After the first update every component is proposed and every
  column taken as before. A mass held by a slider is proposed once in a step
  taken at its first check (before: twice) and differentiated once where
  that check assembles (before: twice); the slider is still asked in every
  pass (`load_line.rs`, bead `.12.5`).
- During the solve the line's flow inputs hold the released flows, and the
  drivers of those inputs carry the extra terms `C * own effort`, set anew in
  every step. The residual of a line input is therefore still the residual of
  the flow that crosses, and Newton runs as for any block, on a map that is a
  function: released flow in, effort out. The line component is differenced
  through `propose_on_line`. Both convergence tests and the refusals are
  the ones above. What a wire is measured by is not: a block with a load
  line has a rule of its own, found piece by piece as refusals of steps
  that were in order (four reviews, each finding the part the fix before
  had left). The audit of a junction knows that junction only and is told
  rounding, never the solve's relative tolerance. The line component is
  differenced, so an update leaves every wire of its block with what the
  differencing left of the step's change. Where the solve then takes a
  wire by another measure than its junction's, the audit refuses what the
  solve took.
  - Every wire of the block is relative to its own size, in both tests,
    and not to the largest of its kind in the block. For a flow on the
    line that is the flow that crosses and what the other components give
    (not the released flow its input holds, nor the terms `C * own effort`
    of its driver: a slider that holds a mass lets nothing cross while the
    released flow is `c` times the push). Refused before: a mass held
    under 4.9 N rising to 4.9196 N in 0.4 s, at its second step,
    -3.2258e-19 J against 6.3977e-20 J (review of ef1631e0d2); the same
    beside a second mass at 0.2 m/s on a spring, 1.7830e-19 J against
    6.4675e-20 J (review of 74529fda7f); a block riding a cart, when the
    friction between them passed through zero, after 248 steps,
    3.7002e-15 J against 1.1119e-15 J (review of 0e3df343f2).
  - A wire that is not a flow on the line is allowed, in both tests, the
    rounding of the largest wire of its kind in the block, `16 eps` of it,
    and the audit is told that for the wire. A wire's own size is not its
    floor: the block's unknowns are solved together, and the rounding of a
    large wire is in the small ones. (The last digit of a held mass's
    force balance of 11.7 N moves a wire of 1e-25 N beside it by a
    thousandth of itself from one iterate to the next; relative to its own
    size alone that block did not converge in 25 iterations.)
  - A flow on the line is allowed the rounding of the difference it is
    (below), and the solve takes it only at half of what the audit is then
    told for it. Taken at all of it, a step at the edge of the solve's
    allowance was over the audit's by what forming the residual rounds off
    (3 N rising to 3.3 N in 20 s: refused after 1728 steps, -2.4416e-20 J
    against 2.4397e-20 J; review of 74529fda7f).

  What this can hide. The audit multiplies what it is told for a wire by
  the wire's paired output. For a wire of a block with a line that is
  `16 eps` times the largest of its kind in the block times that wire's
  paired output. A junction whose efforts are below about 3.5e-6 of the
  block's largest is therefore audited to that rounding and not to 1e-9 of
  itself, and a junction at rest gets nothing from it. It is the rounding
  of what the wire is solved with, not a tolerance. A block without a line
  is measured and audited as it was. (0e3df343f2 allowed every junction of
  a block with a line the rounding of the block's largest power instead,
  flat: a junction below 3.55e-15 of that power was not audited at all.
  Removed, with its accessor `last_block_rounding_w`.)
- When the block is accepted, the flows on the line, `released - C effort`,
  are written to the line's inputs. The audits, the ledger, `last_input` and
  the next step's start see what crossed, not what the solve iterated on.
- Such a flow is a difference of two terms and carries their rounding, not
  that of its own size: a held slider's flow is `r - c (r / c)`, zero for
  some numbers and a unit in the last place of `r` for others. That
  rounding, `16 eps (|released| + sum |C effort|)`, is in the allowance of
  the flow's wire in the solve, and so in what the interface audit is told
  for a line input, on top of what the solve allows every wire. (Told
  nothing, on an earlier kernel, the audit refused a mass held at rest for
  282 of 2000 pushes:
  6.8e-24 J against 6.8e-33 J.)
- Within a row of the component's relation the step a block accepts does
  not depend on `C`: the component settles where the flow on the line is
  the flow its neighbours give, whatever passive matrix it was told, and
  only the path of the iteration differs. So no closed form of a motion
  within a row can show a wrong `C`; a test reads the matrix directly
  (below), and another publishes flows through its off-diagonal entries.
  Where the relation jumps from one row to another, `C` decides. The
  component is handed the released flow `flow + C effort`, so between two
  rows whose efforts differ by `d` the released flow of one and the same
  state differs by `(C - c) d`, with `c` the neighbours' true answer. With
  `C = c` the rows meet. A `C` larger than that leaves a band of states
  with no point to settle on, where the solve cannot converge; a smaller
  one leaves a band with two, and the iteration's path picks. A slider has
  two such jumps: at break-away when its static limit is above its sliding
  force (`d` is their difference), and, stuck with a start slip, between
  the force that ends the step at zero slip and the hold on the mean,
  whatever its coefficients. (Read off `Slider::on_line`. ef1631e0d2 said
  the accepted step does not depend on `C` at all, and 74529fda7f named
  the first jump only; each found wrong by the review of that commit.)
- What `C` is: the neighbours' answer at the step's start as the kernel
  has it. It is exact for neighbours that are linear within the step and
  declare their tangent, and the slope at the step's start for others that
  declare one. A neighbour that is differenced gives it to the accuracy of
  the difference, each entry resolved against the rounding of its output or
  counted (below, "Differenced columns"). (74529fda7f said `C` is
  exact for every neighbour that is linear within the step. Found wrong by
  its review. Until 2026-11-12 an entry whose change the output's rounding
  took was handed as it came, a count of rounding units or nothing, a stated
  limit: beside a spring of 400 N/m added without its tangent, a mass at
  rest handed the slider `dt / (2 m)` = 0.00125 where the two answer with
  0.0012484, a `C` too large, and the step in which the mass should break
  away did not converge; a mass of 2 kg moving at 1 m/s under 1e-6 N handed
  a line no mobility at all.)
- Where the component's relation jumps, the differenced column is not a
  slope. The kernel differences the line component at a released flow
  moved up by `sqrt(eps)` of it. In the window that wide below a jump, on
  the side the move goes to, the moved proposal is on the other row, and
  the column is the difference across the jump, and an update made with
  it leaves the force wire with far more than rounding. It is the
  remainder of differencing, largest here; no step is refused for it, and
  the kernel is not changed for it. A tangent the component declares for
  the row it takes would remove it; the tier has none yet. (Found on
  2026-11-12 by a probe written to show a stall there, which there is
  not: a slider that holds 9.8 N and slides against 5.88 N, under a push
  that creeps up to the limit at 1e-9 N a step of 1 ms on a mass of 2 kg.
  While the component was first asked at the released flows of the step
  before, inside the window (1.46e-7 N, `sqrt(eps)` of the limit), pushed
  the way the move goes, the slider's force was off the push by up to
  8.4e-8 N, alternating in sign, the mass crept at up to 2.1e-11 m/s, and
  it let go with the push still 6.65e-8 N inside the limit, 7e-9 of the
  limit. Since every input of the component is released before its first
  ask, it is first asked at the flow the mass gives it, and that scene no
  longer reaches the limit: the mass is held both ways, inside the window
  and out, to the push's end, 5e-8 N inside the limit, the force within
  1.8e-15 N of the push. Asserted in `tests/load_line.rs`.)
- `C` need not be symmetric. The step is refused (`SystemError::LoadLine`,
  naming the component and the reason, nothing published) when the other
  components do not determine a line (their Jacobian, with all load-line
  owners' columns left out, is singular or not finite), or when the symmetric part
  of `C` has an eigenvalue below minus 1e-6 of its largest entry: the
  neighbours then answer an effort with more flow, not less, and a relation
  such as friction has no point to settle on. The margin is what `C` is known
  to when a neighbour is differenced.
- A flow of the line that the block does not solve for (a prescribed end)
  has a zero row: it does not answer the component's efforts. Its column is
  zero too with the connections there are, because nothing in the block
  reads the effort returned at a prescribed end. A component all of
  whose line flows are given has no line: it is stepped by its mandatory tier
  at those flows, and if it is in a loop through another port (a slider
  between prescribed velocities whose heat warms a body) it is differenced
  there as a function of its other inputs.
- Several components may be on load lines in one block. Their relations are
  solved together; declaring two rubbing surfaces on one body is no longer
  an admission refusal. This does not select a force distribution for
  redundant active constraints: a singular interface Jacobian returns
  `SystemError::Singular` with the block's component names and publishes
  nothing. A fixed point already satisfying the equations can be accepted
  without factorization; no uniqueness claim is made. `AdmissionError::TwoLines`
  remains an API variant but is no longer emitted by sealing.
  `add_on_line` still refuses a line that is not over power ports receiving a
  flow, in increasing order (`AdmissionError::LinePorts`), and claims no name.
- The end of a step on the line. A port's flow is its mean over the step,
  and the elements of this crate end a step at twice that mean less its
  start (the midpoint rule); a component that is exact for a held effort
  does not. What a component promises is `Component::end_of_step`
  (crates/fs-ports/CONTRACT.md), and at sealing the kernel holds every port
  of a load line to it:
  - A line port that names no end input (`LoadLine::line_ends`) is joined
    only to neighbours whose flow ends the step by the midpoint relation.
    `Slider` names none.
  - A line port that names one is joined only to neighbours that tell the
    end of their flow. The kernel makes the wire: the end input is driven
    by the sum that drives the port's flow, taken over the neighbours' end
    outputs, with the same signs and ratios. Nothing else is connected to
    it.
  - The end input is an unknown of the block like the flow, so the factors
    and back-substitutions that give the line also say how it answers each
    effort. The line keeps that as one more row per end coefficient, after
    its square matrix; the input's driver carries the row as the flows'
    drivers carry theirs, the component is handed the released end flow,
    and the end that crossed is written back when the block is accepted. No
    proposal is made that was not made before, and a block without end
    inputs runs the code it ran.
  - Refused (`AdmissionError::LineEnd`, the port named, and what does not
    go together): a neighbour that answers `Untold`; a neighbour that tells
    its end beside a port that takes none; an end input beside neighbours
    that all keep the midpoint promise; neighbours of both kinds in one
    port's flow; an end input that is already connected, or is not a signal
    input of the port's width; an end input on a port whose flow is
    prescribed (a prescribed flow is held over its step, or the component
    is told its span).
- A block with one or more components on lines keeps no factors between steps and
  costs, per step, one assembly of the neighbours' columns (a differenced
  column one proposal of its component, and up to five more where it is
  resolved, below), one factorization and one back-substitution per effort,
  on top of its solve.
- Nothing of the line crosses a step: the matrix, the drivers' extra terms and
  the released flows are all formed in the step that uses them, so a run and
  its replay from a checkpoint agree (test below).
- The kernel locates no event. A component on a line may: it answers a
  proposal with `Proposal::Shorten` and the step that ends at the event, as
  a friction component does whose slip would pass zero inside the step
  (`fs_tribo::tangential::component::FrictionInterface`). `elements::Slider`
  does not: a step in which its relation switches (an arrest, a break-away,
  a reversal) is first-order accurate, and what happens in it is its own
  law of a step.

**Differenced columns.** A component that declares no tangent is
differenced: one input moved by `sqrt(eps)` of its scale (its own size or
the largest of its group in the block, one unit where both are nothing), the
column the change of its outputs over the move. Where the input is small and
an output is not (a mass that moves with little force on it, a spring at
rest that holds a load), the move changes the output by less than its
rounding, and the entry is a count of rounding units, or nothing. Every
differenced column is resolved against that rounding: those that form a load
line's matrix, the neighbours' at the step's start, and those a block
solve's updates are made with. The shares below are numerical policy of the
differenced branch:

- If the initial forward probe is refused, the kernel tries the same-sized
  backward probe: a material law can admit an inward step from the edge of
  its support without admitting the kernel's outward perturbation. Cancellation
  is returned immediately. If neither side steps, the first typed owner error
  is retained when present; when both request shortening, the shorter request
  is retained. A refusal at the actual iterate still rejects the trial.
- An entry is resolved when its change is at least `2^-32` of its output,
  `2^20` units in its last place: a change of `2^20` units puts the column's
  rounding under about 1e-6 of the entry. An entry whose output is nothing
  is resolved.
- While one is not, the column is differenced again on the same side of the
  iterate, at the step that aims the least resolved change at `2^-26` of its
  output, `2^26` units, `2^6` times what resolves it (the change of a linear
  answer is in proportion to the step; for a change that rounded to nothing,
  `2^26` times the step), at most four times, and never at a step longer
  than `2^-10` of the input's scale, or of one unit where that is less. Each
  longer step is one more proposal of the component.
- Each entry keeps the difference that first resolved it: a column that the
  first move resolves is the column it was, to the bit, and an entry of a
  curved answer is not taken at a step longer than its own needs.
- An entry resolved at a longer step is differenced at the longest of them
  on the other side of the iterate as well, one more proposal, and kept if
  the two agree within `2^-10` of it. Where they do not, the output does not
  answer in proportion over the two steps (a kink inside one of them), and
  the entry is left as the first step had it: never the secant across the
  kink. So is an entry that the component does not take on the other side
  (it asks to shorten there, or refuses). A kink whose secant lies within
  `2^-10` of the slope on the other side is taken that close to it. An
  output that rounds by large terms differs on the two sides by its rounding
  as well, and can be left by this (below, "No-claim boundaries").
- A longer step that changes an output by more than `2^-13` of it has gone
  past where the output answers in proportion (a jump), and the entry keeps
  what it had: a jump the other way as far on the other side would agree.
- An entry still unresolved at the longest step, or where the component does
  not take a longer step, keeps its latest difference.
- Each entry left unresolved is counted in the step's report
  (`StepReport::unresolved`, the trials of the step summed). An output that
  does not answer the input at all is among them, at every step: the kernel
  cannot tell it from one whose change rounds away, and its entry, nothing,
  is right. The step is taken and audited all the same; a caller that needs
  the line component's answer to more than the step's balance declares the
  neighbour's tangent, which the line is handed as declared, to the bit.
- The block solve's own columns are resolved the same way, since
  2026-11-14. An entry off by about `1/N` of itself, after a change of `N`
  units, leaves about `1/N` of the residual after an update; `2^20` units
  close the first test's `2^-40` of a scale in the two updates an exact
  block takes, where `2^10` would take four or five. Differenced once,
  as they were before, an update made with an entry the rounding took is a
  fixed-point iteration in that entry's place, which need not contract:
  two capacities near rest behind a conductance (`tests/elements.rs`), and
  three cells of a backward-Euler timeline near rest whose entries were 0 to
  2 units (fs-couple `friction_heat`, two contacts sharing a cell, refused
  at its step 19 at the iteration limit). Most columns differenced again are
  ones with an output whose change is exactly nothing, which a longer step
  cannot tell from a change that rounded away, at one proposal each.
  Counted on ts2 with a probe build, 2026-11-14: the fs-system tests
  differenced 1059209 block columns again in 1213654 block solves, at
  1309118 more proposals and about the same run time; fs-couple's
  `friction_heat` 5208 columns in 1682 solves (its tests 37 s to 61 s); the
  Hertz graph 165 columns, its patch asked 31 and 51 times where it had been
  asked 16 and 26. A component that declared its
  structural zeros would spare those proposals; none does. (On a scratch
  kernel of 2026-11-12, resolved block columns took the drawn graphs below
  from 695778 updates to 674215 on a Mac; before longer steps were kept off
  jumps they gave a sliding slider's force, which does not answer its
  released flow, columns of 120.5 and 25.2 once a longer step carried it to
  sticking, and the drawn graphs took 790554.)

**Requests to shorten inside a loop.** A component of a looped block is
proposed at every iterate of the block's solve, and, if it declares no
tangent, at inputs the kernel has moved to difference it. Not all of these
proposals are the step's. Both rules below were found by the first component
that locates an event (a friction interface that cuts the step at a stop,
2026-11-12); on the kernel before them it could not take a step.

- A perturbed proposal that asks to shorten does not shorten the step. A
  column of a component's tangent is a forward difference: one input moved
  by `sqrt(eps)` of its scale. Where the component takes the step at the
  iterate and asks to shorten at the moved input (the iterate ends a step
  at rest, and the moved input lies beyond the stop), the request is about
  the kernel's probe, not about the step. The column is then differenced on
  the other side of the iterate, `base - step`, where the component takes
  the step. Only if it asks on both sides is the step shortened, to the
  shorter of the two. Where the forward proposal takes the step the
  arithmetic is what it was. (Taken as a request, the probe cut a step
  that ended at a stop a second time, 5e-10 s short of it; the step ended
  with a slip of 1.3e-8 m/s where it should have ended at rest, and the
  step after it was refused, `StepUnderflow`.)
- A flow of the line, or the end of one, that no effort of the line
  component moves (its row of the line's matrix is zero: a side that a
  velocity source in the block moves) is not differenced. Its column is
  held at zero and no proposal is made at a perturbation of it; the
  component reads it as a given end. A perturbation of it is a jump of that
  end that nothing in the block makes, and a component that locates the
  stop of the slip a jump starts asks on both sides for the step that ends
  there: `fs-tribo`'s friction interface stuck to a belt moved by a
  velocity source had every step cut to `2^-27` s, the stop of the
  kernel's own `2^-26` m/s, and crawled (bead `.5.8.1`). Where a line has
  such a flow the arithmetic is not what it was before (2026-11-14): the
  factors hold a zero where the difference was, and the rounding estimated
  for the rows that read the component lacks that column's part.
  `tests/load_line.rs`: a slider beside such a side is asked as for one
  differenced column, beside a free mass as for two.
- The component on a load line is first asked at the released flows of the
  step. A block's solve starts from the flows the last step ended with, and
  before the first step from none; for the line's component those are not
  the released flows of this step. A component that enters a step with a
  slip those flows show nothing of (a body thrown onto it with no starting
  flows given) would see a stop inside the step at them, and the shorter
  trial it asks for would start from the same values and ask again: the
  first step of a mass of 2 kg thrown at 1 m/s onto 3 N of friction was
  refused after 64 shortenings. (7af293b67e named a prescribed side that
  has jumped as a second such case, here and in its message. It is not one:
  the new flow of a prescribed end is among the step's inputs from the
  start, and a block on a side that jumps is dragged and caught on the
  kernel before that commit as after it. Not run before it was written.) So
  before any load-line owner is first asked in a trial, every input of each (the
  line's flows, the ends of its flows, and any other, such as the span of a
  side it is told) takes one update with every owner's efforts held at what the last
  step accepted, by the factors the search for the line has just made (the
  block's Jacobian without any load-line owner's columns). That update gives
  the released flows of the step, and what the others give its other
  inputs, exactly where the other components are affine, and their first
  Newton iterate where they are not. Every other unknown of the block keeps
  its starting value, a `guess` included. An `elements::Slider`, which
  never asks, takes the first row of its law at those flows as well.
  (Until this rule the update was made only after the component asked to
  shorten at the starting values. A component that answered there was
  answered at the last step's crossing flows plus the line times its last
  efforts, and the block's first update corrected it: in a bowed string's
  run every trial's first ask was there to the bit, and a held step took
  1.998 iterations on average. With the rule it takes 1.000, chain fsv.)
- What a component asks at any proposal of a trial, the first included,
  shortens the step at once. An event located at an iterate that has not
  converged is located from flows that are off by what the solve has left
  there; the shorter step is taken or cut again, and the component locates
  from there. Not covered: a bound on how often.
- A block in which no component asks to shorten from a perturbed proposal,
  and whose line component does not ask at the starting values, runs as it
  did, to the bit (the example scenes and every test that was there).

## Invariants

- Nothing is published unless every component proposed, passed the audit and
  admitted. A failed step leaves every accepted state, the clock, the ledger
  and the warm-start values bitwise unchanged.
- A step repeated after a failure is bitwise equal to a first attempt.
- Mass measurement and both owner and connection mass audits precede every
  publication. Unaccounted measured mass beyond the allowed roundoff is
  refused by name as `Contract` with `ContractViolation::MassBalance`, even
  when energy balances. Uncancelled connection mass beyond tolerance is
  `InterfaceMass`. Either refusal leaves all states, clock, warm starts and
  ledger rows unchanged, including nodes audited before the failing one.
- Entropy measurement and the owner entropy audit also precede every
  publication. A declared change that port entropy and booked generation do
  not explain is refused by name as `Contract` with
  `ContractViolation::EntropyBalance`, with the same guarantee: nothing is
  published.
- Junction equations use coefficients of exactly plus or minus one, so the
  connection itself adds no rounding. A junction always closes a loop between
  its setter and its members, so its interface energy cancels to the solver
  tolerance, and the remainder is reported per connection. A one-way heat
  route multiplies by fractions and cancels to roundoff.
- A refused builder call changes nothing: every port of a connection is
  checked, including for use twice in the same call, before any is claimed.
- A transformer or power map uses an adjoint pair of linear drivers in the
  same block plan, Newton Jacobian and energy audit as a junction. Its powers
  cancel to roundoff and interface-solve tolerance; no integration is added
  to a device. A rejected ratio, map shape or coefficient claims no ports.
- A convergence test compares each residual with the size of the terms its
  driver sums, not with their total, so a balance of large efforts that
  nearly cancel in the driver can still converge. Each scalar signal input is
  judged against its own size only.
- A residual is also allowed 16 units of its first-order rounding,
  `eps * sum |d residual / d input| |input|` over the block's own inputs,
  taken from the tangents of the previous iterate. Without it a block could
  not converge once a component's output became a small difference of large
  inputs: the heat rate between two bodies near one temperature, the force of
  a spring whose ends move together. The interface audit allows the energy
  that this rounding estimate and the configured absolute tolerance amount
  to, on top of its own relative tolerance; it is not told the solver's
  relative tolerance, so a solve accepted too loosely is still refused there.
  Inputs the block does not solve for (boundaries, earlier blocks) are not in
  the estimate.
  - A step taken at its first check tells the audit the estimate of the
    assembly that judged it, made at its own iterate (see the second test
    below). It used to tell the audit the allowance formed before any
    assembly, which has no rounding in it, and then a junction whose flow
    had gone to zero was refused for a rounding-sized difference of large
    terms: the audit's relative tolerance goes to zero with the energy
    that crosses. `tests/elements.rs`: a pad of 50 J/K at 400 K on a
    `HeatNode` of 2 and 6 W/K whose other body port is held at 300 K, with
    30 W on the node, 4000 steps of 0.5 s. The pad follows the midpoint
    rule's closed form to 1e-11 of 305 K at every step, settles at
    `300 K + P / G`, and more than a thousand steps are then taken at
    their first check. Before the change it was refused at step 1660
    (2.5e-14 J uncancelled against 2.2e-18 J); the kernels of 42f2803e8d
    and f5db075f3b refuse it too, at steps 1485 and 1660.
  - What such a step tells the audit is that estimate and nothing looser.
    A second test brings the same graph to rest, plants a heat rate on the
    pad's input and loosens the solve to 1e-6 so that the step is taken at
    its first check and the verdict is the audit's. The audit's tolerance,
    read from the refusal of a planted microwatt, is between 1 and 1.5
    times 16 units of `eps` times the terms the pad's heat depends on
    through the plane (1.5 W/K of 305 K and a quarter of 30 W, 1.65e-12 W
    in all); a hundred times that is refused by that junction's audit by
    a factor between 80 and 120; a quarter of it is taken with no update.
    Against a scratch kernel that tells the audit a thousand times the
    estimate the test fails.
  - Not covered: a junction at rest whose large terms all come from inputs
    the block does not solve for has nothing in the estimate and can still
    be refused there.
- An energy-rate unknown (a heat rate into a receiver, a stream's enthalpy
  rate into a storage) is also allowed 16 units of `eps` times the largest
  power crossing a power port of its block. The first-order estimate above
  misses a heat rate that is a product of two rounding-sized factors: a mass
  held at rest by forces that cancel moves at the rounding of those forces,
  and its damper reports the square of that, which has no size of its own to
  converge relative to. The interface audit and the second-law check are
  told the same allowance (as energy, and as that energy over the receiver's
  temperature). The allowance is an absolute floor of about 4e-15 of the
  block's largest power; a heat rate smaller than that is not resolved.
  `tests/step_cost.rs`: a mass held by a spring against a constant force,
  its damper sharing a body's temperature, steps 40,000 times to rest (it
  used to be refused after about 20 s at a residual of 4.6e-37 W), ends at
  `F / k`, and the body holds `F² / (2 k)` to 1e-10.
- A residual is judged a second time, by what it does. A wire that passes
  the test above can still be too open for what reads it. An effort is a
  gauge quantity: a temperature closed to 1e-12 of 300 K is open by 3e-10 K,
  and a conductance that reads it returns a heat rate wrong by that much
  times its conductance, while the heat wire agrees with what the
  conductance returned and passes too. So the residuals are carried once
  through the components that read them, by the slopes
  `t_ij = d(target of wire i) / d(input j)` of an assembly of the block's
  Jacobian made in the same step, and for every wire `i`

      sum_j |t_ij| max(0, |r_j| - f_j)  <=  rel S_i + f_i

  - The first test is `|r_i| <= rel S_i + a_i`. `S_i` is the largest, over
    the block's inputs of the wire's kind (a carrier and side; each
    coefficient of a stream port; each scalar signal alone), of an input's
    own value and the terms its driver sums. `a_i` is the absolute
    tolerance, plus 16 units of `eps` times `sum_k |t_ik u_k|` from the
    step's most recent assembly (nothing before the step's first), plus for
    an energy rate 16 units of `eps` times the block's largest power.
  - Only the part of a residual above its own rounding floor `f_j` is
    carried. The floor is the absolute tolerance plus 16 units of `eps`
    times the sum of two sizes: the larger of the wire's own value and the
    terms its driver sums, and `sum_k |t_jk u_k|` (for an energy rate also
    the block's largest power). Without the floor a force balance near
    rest, which sits at the rounding of forces that cancel, would keep a
    body's velocity wire from ever closing.
  - Both sums over slopes run over the stored entries of the assembly, one
    per term of a driver that reads an output of the component an input
    belongs to. Where two terms of one driver read the same component,
    their magnitudes add, in what is carried and in the floor alike.
  - What it bounds, for every wire: its residual, and the first-order image
    of all residuals under one more pass through the block, by
    `rel S_i + f_i` plus each in-floor input's floor carried once,
    `sum_j |t_ij| f_j`. For a heat-rate wire that bound times the step is
    the energy a step can leave between a component and what it is joined
    to; for a power port, times the conjugate effort too.
  - What it does not bound: the whole first-order error `(I - T)^-1 r` of a
    loop whose gain exceeds one. Measured on the kernel of f5db075f3b at
    the accepted iterates of every graph that existed then, in an
    instrumented scratch copy that is not in the tree (73 tests and 13
    scenes with a looped block): the full correction `J^-1 r` exceeds a
    wire's tolerance in six tests and one scene. Once it is on a wire with
    a physical size, 1.8 times (2.9e-11 W on a heat wire of the friction
    graph's sixth case). The others are systems at or near rest whose
    wires are below 1e-15 in absolute terms, where the ratio reaches 4e10
    with a correction of 2.8e-21.
  - The slopes are the step's own. At a step's first check they are formed
    at the iterate being judged: the block is assembled once, and if the
    test passes the step is taken with no update and counts no iteration.
    At a later check of the same step they are those of the assembly that
    produced the update, one iterate earlier. A block whose residuals are
    all exactly zero is taken without assembling. An affine block that
    holds kept factors for the step length (below) adopts them and their
    slopes in place of an assembly. A check is two passes over the stored
    entries, a third when it refuses, and calls no component.
  - Differencing a component writes a perturbed proposal to its candidate.
    A block taken at its first check therefore proposes each component
    that has no declared tangent once more at the accepted inputs. A
    declared tangent writes to scratch only.
  - What it costs: a block that is not affine and is taken at its first
    check pays one assembly in that step, and one more proposal of each
    differenced component, where the first test alone paid neither. An
    affine block pays neither while its step length stays the same; it
    keeps the factors of one length, so a run that alternates lengths
    assembles at every change. `tests/elements.rs` counts the proposals
    asked of the warming plane's differenced source in each of 80 steps:
    three in a step taken at its first check (the check, the assembly that
    judges it, and one more at the accepted inputs), `2 k + 1` in a step
    of `k` updates; and the source publishes, to the bit, its proposal at
    the temperature the step accepted. An exact fixed point would ask one.
    Run against a scratch kernel that judges a first check by slopes left
    in the workspace without assembling, that test fails at step 27 with
    one proposal; the replay tests below pass on that kernel, so they do
    not tell which slopes made a verdict, and this count does.
  - What a block keeps from an earlier step does not decide a step. A
    checkpoint does not carry a block's workspace and a discarded trial
    does not roll it back, so a verdict or an update that read it would
    differ between a run and its replay. Two things did. f5db075f3b judged
    a step's first check by the last step's slopes; and since 20c8e96867
    an affine block with a differenced part kept factors made from
    difference quotients, which are rounded at the point they are taken at
    and are not the numbers a fresh assembly forms (see the bullet on
    affine blocks below). `tests/elements.rs` and `tests/step_cost.rs`
    compare a run with a copy that takes a step of another length and is
    restored before every step, in iterations and in bits: the warming
    plane below over 80 steps, of which more than ten are taken at their
    first check (against f5db075f3b the run took 0 iterations at step 27
    where the copy took 1); and a mass pushed against a spring that is
    added without its tangent, over 2,000 steps (against 42f2803e8d and
    against f5db075f3b the copy took 2 iterations at step 1 where the run
    took 1). A third, a loaded mass to rest over 12,000 steps, is a block
    that is never taken at its first check and is not affine: it passes on
    those older kernels too and guards neither defect.
  - `tests/elements.rs`, at the default settings: two bodies on a
    `HeatNode` whose source falls with temperature, each losing heat to a
    room through a film (the warming plane), 80 steps of 0.25 s. After each
    step the test solves the step's own equations for the states the
    kernel published and measures each body's balance in joules: 6.6e-12 J
    at worst (allowed 1.6e-11 J). The same test against the kernel of
    42f2803e8d: 5.2e-10 J at step 73. The same measure over 240 steps
    whose length changes at every step (0.25, 0.125, 0.0925 and 0.5 s in
    turn): 1.3e-11 J at worst (allowed 2.7e-11 J), with 150 steps taken at
    their first check and 101 iterations; against 42f2803e8d 1.1e-9 J at
    step 139; f5db075f3b, which had no slopes at a new step length, took
    none of those steps at its first check and spent 251 iterations.
    `crates/fs-couple/tests/friction_heat/graph.rs` runs at the default
    settings: the hand-written module, which accepts 1e-10 J, refused the
    third case's ninth step with 2.4e-10 J before the second test and sees
    4.1e-12 J at worst with it.
  - A block refused through this test alone names the open input that
    moves the furthest-moved wire most, with that input's residual
    (`NotConverged::worst`). `tests/elements.rs`: two bodies 1e-10 K apart
    through a wall of 1e6 W/K, both of the wall's temperatures guessed at
    300 K, no iteration allowed: the refusal names `wall.b` with a
    residual of 1e-10 K, inside the first test's 3e-10 K. 42f2803e8d takes
    that step with no iteration; f5db075f3b refuses it and names `wall.a`,
    whose residual is zero. And at a check after updates: a body whose heat
    capacity rises with temperature behind a wall of 1e5 W/K, the far side
    3 K warmer, needs four updates; allowed three, the refusal names
    `wall.a`, open by about 4e-11 K. A scratch kernel that does not carry
    the open input to the refusal at a later check names `body.heat`
    there, the block's first unknown, which is closed. An input is ranked
    by all its entries on the moved wire together.
  - What changed for existing graphs when the test came in (A/B of
    f5db075f3b on 42f2803e8d, every target that uses the kernel): the 14
    example scenes gave byte-identical reports with the same iteration
    counts; no test that pins an iteration count changed.
- A looped block keeps its workspaces between steps, and an accepted step
  allocates nothing (`tests/step_cost.rs`, counted per thread, damped and
  undamped lines). A refused step may: it drops the step's lists and builds
  its error.
- The linear system of a Newton update is the interface Jacobian in every
  scalar input of the loop, but it is not eliminated densely. At sealing the
  inputs are split by the pattern of that Jacobian (row `r` has an entry in
  column `j` when the driver of `r` reads an output of the component `j` is
  an input of). An input *follows* if nothing it depends on depends on it;
  where inputs depend on each other in a ring, the ones that depend on the
  most members of the ring are *kept*, no two of them neighbours in one
  round, until no ring is left, and of a ring that is one plain cycle one is
  kept. The followers are substituted in dependency order (their part of
  the matrix is unit lower triangular), and the kept inputs are eliminated
  together with partial pivoting. This is the same linear system, so the
  Newton iterates are those of a dense solve up to the rounding of the
  elimination, and every input keeps its residual, its tolerance and its
  place in the audit. For a line of masses, springs and dampers the kept
  inputs are the force sum on each mass and the body's heat sum; for modes
  that meet at one contact, the contact's rate alone, which is the scalar a
  hand-written impact loop solves for. `tests/step_cost.rs` asserts those
  counts through `plan_loops` (counts, not which inputs: with one damper
  the tie between its temperature and the body's heat sum goes to the
  lower index) and compares each graph with a hand-written loop to 1e-12.
  It also pins the iterations: one a step for the undamped line and the
  modes, at most two for the damped line, and for a damped line whose heat
  is routed to two bodies at their entropy-equivalent temperature (a driver
  that is not linear, its slopes remade for every assembly) exactly the
  iterations of the line whose dampers share one body's temperature.
  `src/solve.rs` tests the elimination against dense elimination on named
  patterns (lines, a star, rings, a chain, a full pattern, a follower
  reached both directly and through others) and on 600 drawn at random,
  half of them shaped as the kernel makes them; the counts the rule gives
  where keeping every candidate at once, or counting what an input depends
  on outside its ring, would give others; its refusal of a singular matrix;
  and that a plan refused or factored without a solve is as good as new
  once cleared. The rule that picks the kept inputs is not a search for the
  fewest: keeping more than needed costs time only.
- An assembly that is cut short leaves no factors behind: the kept step
  length is voided before the first entry is written. `tests/step_cost.rs`:
  in a line of six masses one spring refuses to give its tangent on demand;
  thirty times in ninety steps the line is asked for a step of another
  length while it refuses (the refusal names the spring, and clock and step
  count stay), and every step of the usual length that follows equals a
  line never interrupted to the bit, in one iteration. Cut short by a
  refusal only: an assembly ended by a singular matrix takes the same exit
  in the code and is not run by a test. (A request to shorten can end an
  assembly only through a differenced component, and a block with one keeps
  no factors.)
- A block whose components are all `Class::ExactLinear` and declare their
  tangents, and whose drivers are all linear, is affine: its interface
  Jacobian is taken to depend on the step length only. Its factors are kept
  and used for the first update of every later step of the same length, and
  their slopes for that step's second convergence test; the rounding
  estimate is remade from the kept entries at the current inputs, in the
  order it was first formed. A block with a component that has no declared
  tangent keeps nothing, whatever its class: the kernel differences that
  component, and a difference quotient is rounded at the point it is taken
  at, so factors made from it in one step are not those of another. A
  component whose slopes move with its state must declare `Class::Smooth`.
  If one declares `ExactLinear` all the same, what the block keeps is
  stale: an update that falls short is caught, and the next iteration
  derives the Jacobian as before, but the second convergence test judges
  by the stale slopes, so a step can be taken outside that test's bound,
  and a run and its replay need not agree. No component in the tree is
  like that: of the ten `ExactLinear` elements of `src/elements.rs` and
  `src/elements/node.rs`, nine have tangents that take neither the state
  nor the inputs, and the spring, which is `ExactLinear` only with a zero
  cubic coefficient, reads them in a term multiplied by that coefficient.
  That was read from the source, and it says nothing of the three elements
  of the class added since: `DiracJunction`, `ControlledSource`, and
  `ConductanceMatrix` without a law. The conformance kit now holds a
  component to the claim: fs-ports' `check_port_tangent` asks
  `check_exact_linear` of every component of the class, which takes the JVP
  along each unit input direction at the case and at the case's stale state
  with other inputs, and refuses a slope that differs in any bit. This
  crate's tests ask it of all thirteen elements of the class, at 5,486
  pairs of points in each feature set when this was written (counted with a
  kit made to print the type it is asked of). A mass and a heat capacity
  whose slopes were scaled by `1 + eps x` of a state `x`, and a conductance
  whose slope was scaled so by an input, each fail `tests/elements.rs` at
  that law. The law compares two points at a time and is no proof; it asks
  nothing of a component with no declared tangent, or of one whose tests
  do not run the tangent tier. The
  kept factors solve to the bit of the elimination they replace (unit test
  in `src/solve.rs`, against an elimination that carries the right-hand side
  along), and a line of masses and springs stepped with them equals, bit
  for bit over 400 steps, a copy forced to derive them afresh for every
  step by a detour through another step length and a checkpoint
  (`tests/step_cost.rs`).
- Once cancellation has been requested, any refusal during a trial is
  reported as `SystemError::Cancelled`.
- A graph without loops is evaluated once per component per step, in
  dependency order: it reproduces a hand-staggered run of the same components
  bitwise.
- Every port that exchanges power or needs an input is in exactly one
  connection or boundary, except a heat receiver, which may be in several
  heat routes; `seal` refuses anything else.
- A sink named by several heat routes is one ledger row; naming it again at
  another temperature is refused (`SinkTemperature`).
- Numerical loss is accounted separately from physical dissipation.

## Error model

`AdmissionError` for descriptions that cannot run, including
`TransformerRatio`, `MapShape` and `MapCoefficient` for invalid power maps,
`Empty` for a junction or a power map without members,
`HeatRouteWidth` for a heat route that names a field heat port,
and `LinePorts` and `LineEnd` for load lines. `TwoLines` remains an API variant
but multiple load-line owners are admitted.
`SystemError` for a failed
step: the component's own error unchanged (`Component`, with the phase), a
contract breach (`Contract`), a bad negotiation answer, `NotConverged` naming
the block and the worst input, `Singular`, `Unsupported` naming a component
that is in a loop without a tangent, is not smooth and was not added on a load
line, `LoadLine` naming a component whose neighbours give it no usable line,
`GuessRequired` for a
temperature in a loop without a positive starting value, `Interface` for
uncancelled energy, `InterfaceMass` for uncancelled stream mass, `SecondLaw`,
`StepUnderflow`, `Step`, `Cancelled`. An invalid description or a
failed step never panics. A `PortRef` or `NodeId` taken from a different system
is a caller error and is not checked.

## Determinism class

Deterministic. Sequential; components, connections and blocks are visited in
index order; no hash-ordered container, clock or random number.

## Cancellation behavior

`step` polls `Cx::checkpoint` on entry, before every block, in every Newton
iteration and immediately before publication, and returns
`SystemError::Cancelled` with nothing published. Components poll the same `Cx`.

## Unsafe boundary

None. The crate denies `unsafe_code`.

## Feature flags

- `constitutive`: source-owned nonlinear electrical storage; enables the
  optional `fs-phs` dependency. Off by default. `tests/constitutive.rs`
  checks conformance and tangent action away from knots, named support
  refusals, nonlinear RC relaxation with heat, and nonlinear LC exchange
  against independent piecewise or continuous-equation references. At its
  five conformance points both adapters pass `check_time_reversal`: `R` is
  the identity on charge and negates flux linkage. The bounds are derived in
  the test from the rounding of fs-phs's scalar curve, for a curve whose
  knots, slopes and knot energies are exact in binary. A capacitor reversal
  that negates the charge is refused at the outputs, and so is an inductor
  whose curve is not odd: reversed, its motion follows `-I(-flux)`, which is
  another law.

`slices` (off by default) adds fs-mbd and fs-thermal as dependencies. It
enables the owner slice tests, which compose real physics owners on the
kernel, and `lowering::owners`, the only library code that uses either crate.
Without the feature the crate and its tests build only fs-ports and fs-exec.

`modal-slices` (off by default) adds fs-nlmodal as a dependency for the owner
slices of its linear modal structure and its joint loaded basis, apart from
`slices` so that neither build pulls in the other's owners. The library
itself does not use fs-nlmodal.

## Conformance tests

`tests/adaptive.rs` compares error-controlled oscillator trajectories against
closed-form motion, checks tighter-tolerance refinement, and compares retained
states, source clocks, warm starts and every ledger row with two explicit
half steps after an error rejection. It also exercises a late half-step
shortening, invalid norms, exhausted attempts, cancellation after both fine
steps, rollback/retry, and per-coordinate tolerance scaling.

`tests/checkpoints.rs` checks reusable restore/replay of a driven dissipative
graph, source and kernel clocks, interface warm starts, every ledger field
and solver configuration; restore after failed or cancelled trials; and
foreign-system refusal without mutation. Checkpoints are in memory and bind
to one sealed system's allocation identity; they are not durable snapshots.

`tests/snapshots.rs` (bead `.5.4`, item 4) seals two native graphs and drops the
sealing system: the device of `tests/checkpoints.rs` at step 73 (a sine force
whose account holds a clock, a moving support, damper heat routed to a body,
a recorder whose words grow with the run, and a window), and the two vessels
of `tests/elements.rs` with their restriction at step 250. A child process of
the same test binary receives only the envelope bytes (3,560 and 2,912 bytes
with the ledger of codec version 3, measured 2026-11-24) and the
subject id the sealer kept, on its stdin. The child first refuses a wrong
subject and a flipped byte with nothing changed. It then shows, to the bit,
what the uninterrupted run shows after the split and after each of the 184
and 350 later steps: time, step count, every state, the last values on each
port, and the ledger. In one process, one snapshot resumes into two systems
that then go on alone, and it rewinds the system that wrote it. These
refusals each leave the system showing what it showed before: another graph,
a spring stiffer by one ulp, a member under another charter, a member without
a codec (which also cannot be sealed: `NoCodec`), a last member whose decode
refuses another width (`Decode`, after the members before it decoded), and a
context of another problem at sealing. Truncated bytes, and a byte flipped at
the start, middle or end, never open. A gas constant one ulp lower or a
conductance one ulp higher is another problem. The image charter derives its
pinned ids: codec version 3, the state type of version 1, version 2's schema,
a codec other than version 2's, and a schema and codec other than version 1's.
An image sealed under version 2's identity holding a current image's words
opens as version 2, a version 3 reader refuses it as `WrongStateSchema`, and
the system that sealed it is unchanged. Mechanical, thermal and source elements use their
production codecs, without test-only serialization wrappers. Native codec
regressions also change every core physical coefficient, every waveform
parameter and knot, the power domain and optional signal layout; a resumed
system refuses these changes without mutation. Invalid state and width tests
exercise the same native codecs directly, including source-clock low bits.

The snapshot suite also resumes a belt-driven friction/heat graph after a
velocity reversal, preserving all native states, source clocks, port values and
ledger rows through later sticking and sliding. The slider's refusal-history
counter is separately retained across its codec. General junction tests bind
the original coordinate frame as well as incidence coefficients. With
`constitutive`, the nonlinear LC tests resume the original charge and flux
laws through varied step lengths with bitwise state and snapshot equality;
changing any source knot refuses atomically.

`tests/ledger_accounts.rs` (bead `.5.4`, item 5) drives members without ports
whose receipts state planned stored changes, so each step's defect is known
exactly. Defects of +0.25 J and -0.25 J leave the signed account at zero and
the absolute one at 0.5 J. A step refused in `admit` leaves the ledger's Debug
text as it was. Over 200 changes spanning 97 bits, the signed and absolute
accounts equal exact integer sums in units of 2^-70 J, and the residual is
that sum rounded once. Two accounts of 1e6 J leave a residual of 2^-40 J,
given exactly, where rounding each account first gives 0. On the turning flow
of the two vessels, every absolute account (energy, mass and each carried
coordinate per component; energy and mass per connection) equals a judge's
sum of the size of each step's change of the signed account, within a
derived double-double allowance, over 600 steps.

G0, `tests/system_g0.rs`, on one-port elements defined in the test:

- a bonded mass and spring against the monolithic implicit midpoint rule over
  5,000 steps, with constant stored energy and a closed ledger;
- a driven mass, spring and damper on a common-flow junction, with the damper's
  heat routed 0.6 / 0.3 / 0.1 to two thermal masses and an ambient sink:
  trajectory against the monolithic rule over 10,000 steps, every joule of
  loss found by the declared split, first-law residual of the window, entropy
  generated by the damper and none by the route;
- a one-way graph (oscillator, boundary force, fixed-temperature route) has no
  looped block, spends no iteration and equals a hand-staggered run bitwise;
- a cubic spring with its tangent hidden converges by differencing and
  conserves energy;
- a conductance at its material law's temperature limit steps inward with
  its tangent hidden, against the declared tangent and an independent cooling
  recurrence; two refused probes preserve the owner's first error, and
  cancellation during a probe publishes nothing;
- a refusal in `admit` by the last component publishes nothing, returns the
  owner's error unchanged, and a retry equals a first attempt bitwise;
- cancellation requested during a step publishes nothing;
- a component bounds the step in advance (`next_event_s`) or during the trial
  (`Proposal::Shorten`), and the whole system takes the short step;
- a non-smooth component added without a load line and in a loop is refused
  by name; a temperature in a loop without a starting value is refused by
  port;
- nine admission refusals;
- a loaded damped oscillator settles to rest over 10,000 steps under the
  default tolerances, and an undamped one started at its loaded equilibrium
  stays there (both stopped with `NotConverged` before the convergence scale
  was corrected);
- a refused junction or heat route leaves the builder usable and leaves no
  sink behind; a heat rate cannot be shared by several members of a common-flow
  junction;
- a signal loop `x = cos(x)` reaches its root beside a signal of 1e12 in the
  same block; two masses on one spring (a three-port common-effort junction)
  follow the monolithic rule;
- `NotConverged`, `Interface`, `Singular`, `Negotiation` (including a NaN
  limit beside a valid event), a bad `Shorten`, and `SecondLaw` are each
  produced and name what failed; cancellation noticed by a component is
  reported as cancellation.

Elements, `tests/elements.rs`:
- every element passes `check_component` and `check_port_tangent` at 200
  generated operating points (three spring laws), with outputs compared to
  independent formulas (mean velocity, linear spring force, damper heat equal
  to absorbed power, conductance antisymmetry);
- time reversal: at each of those 200 points the mass, the three springs and
  the heat capacity pass `check_time_reversal`, 1,000 checks in all, a count
  the test asserts.
  - For the mass, `R` negates the velocity. For the springs, `R` is the
    identity on the extension, and for the heat capacity on the temperature.
  - The bounds are derived in the test from the rounding of `Mass::propose`,
    `Spring::propose`, `HeatCapacity::propose` and their stored energies.
  - A reversal that also negates the mass's position is refused, and only
    through the test's view of the state: the position reaches neither the
    output nor the stored energy.
  - A heat capacity that publishes its temperature is refused before it
    steps: a signal port has no parity under time reversal;
- unit rescaling: at each of those 200 points the mass and the damper pass
  `check_unit_rescaling` in units changed in each base dimension alone and in
  all six together, 2,800 checks to the bit, a count the test asserts.
  - So do a mass at rest under no force and a damper with both ends moving
    together, steps whose every energy is zero.
  - The test builds the mass (kg) and the damping (kg/s) in the new units, and
    scales the mass's position (m) and velocity (m/s).
  - A damper built with its damping left in the old units is refused at the
    outputs;
- subnormal energies: a mass and a damper at speeds of 1e-175 to 1e-155 m/s,
  where the energies of a step lie below the normal range of f64, pass
  `check_component` at 2,000 steps each, a count the test asserts. An
  element's receipt claims only roundoff relative to the size of its terms
  (`64 eps` times their sum); fs-ports' `audit_step` adds the f64 underflow
  allowance;
- refusals: non-finite or out-of-range parameters, empty or non-increasing
  tables, a heater whose power could turn negative, a temperature input that
  is not positive, a step that would cool a body below 0 K (`Shorten`);
- waveform means against composite Simpson quadrature, additivity over
  adjacent intervals, and hand values;
- an anchored mass and spring rotate by `2 atan(w dt / 2)` per step (the
  implicit midpoint rule's exact discrete solution) and keep their energy;
- two masses on a hardening spring keep their momentum, their centre of mass
  moves uniformly;
- two bodies through a conductance relax by `(1 - h) / (1 + h)` per step,
  keep their energy, and the ledger's entropy matches `sum C ln(T / T0)`;
- a heat node: over 200 drawn sources and body temperatures it passes the
  conformance and tangent checks, sits where its balance puts it, returns
  each branch's heat, and generates exactly the entropy that leaves through
  its ports less what enters; it also conforms with its bodies within a
  nanokelvin of one temperature and exactly at it (its first energy
  tolerance was sized by the rates, which vanish there, where the rounding
  is that of the terms `G T`; a graph refused such a step, and the tolerance
  is now sized by those terms); it refuses no conductance, a conductance
  that is not positive, a source that leaves it no positive temperature and
  values that are not numbers, writing nothing. Two bodies through a node of
  2 and 6 W/K with nothing on its source keep step with the same bodies
  through one conductance of 1.5 W/K to 1.1e-13 K over 400 steps (allowed
  2.4e-9 K from the kernel's solve tolerance) and generate the same entropy
  to 1e-12. With 30 W on its source every joule is in the two bodies; the
  cold body's share of the source is 5.70 in the first step (it also takes
  what the hot body gives up) and 0.707 at the end, `C2 / (C1 + C2)`, with
  the node 4.35 K and 3.49 K above the bodies;
- a driven damped oscillator built from the library matches the monolithic
  midpoint rule, the body's account equals the damper's dissipation, the
  force's clock equals the kernel's bit for bit;
- a heater in a body behind a conductance to a fixed ambient approaches
  `T_amb + P / G` by the discrete factor, and heat kept plus heat out through
  the boundary equals heat given.
Three planted defects were each caught: the spring's far end pulling the same
way (momentum and conformance), the damper's heat tangent without its factor 2
(tangent check), a force sampled at the start of the step instead of its mean
(oscillator reference).

Unit tests cover the block ordering (including a 200,000-node chain without
recursion), the equivalent-temperature gradient, the dense solve, the
condensed solve (which inputs are kept for a line, a star, a chain, a plain
cycle, an input that depends on itself, full patterns and a cycle fed from
another ring; agreement with dense elimination to 1e-13 on nine named
patterns and to 1e-12 on 600 drawn at random; the same bits from kept
factors; refusal of a singular kept matrix and of a solution that
overflows; reuse after a refusal) and the compensated sum.

Power conversions, `tests/power_maps.rs`:
- positive and negative gear ratios give the reflected stiffness `ratio² k`
  and the exact implicit-midpoint phase, conserving the oscillator's energy;
- a pressure source drives a piston through its declared area, reproducing
  constant acceleration and pressure-volume work in a different power domain;
- full-rank, rank-one and zero rectangular maps match an independent
  monolithic two-coordinate solve, including conservation and null modes,
  with one Newton update per step;
- invalid ratios, shapes, coefficients, roles, carriers and occupied ports
  are refused before mutation; a repaired builder seals and steps;
- a junction with one ratio gives the transformer's states bit for bit after
  400 steps, for a positive and a negative ratio;
- a mass between two springs behind different levers (one reversed) and a
  levered constant force follows the oscillator with the reflected stiffness
  `r1² k1 + r2² k2` about the shifted rest position, at the exact midpoint
  phase, with one Newton update per step; the springs' extensions are
  `-r_i x` and the source's account is minus the levered force times the
  travel;
- the dual: a capacitor feeding two coils through turns ratios (one
  reversed, with initial currents) rings at `w² = (n1²/L1 + n2²/L2)/C`, each
  coil following the junction current by its own share, energy conserved;
- ratio junctions refuse no member, a zero or non-finite ratio, a wrong
  role, a heat or signal port, unequal widths and a member named twice,
  claiming nothing; the same ports then join and the system seals and steps.

Power maps from several ports, `tests/power_map_members.rs`:
- members `r_k I` on width-two banks (one ratio reversed) give the ratio
  junction's states bit for bit after 800 steps;
- a width-two and a width-one member onto three springs follow the midpoint
  equations of their stacked map, solved by elimination in the test, for
  1,000 steps with one Newton update per step: a full-rank stack, a rank-two
  stack with an unconstrained inertial mode, and a member the map does not
  reach, which keeps its motion exactly; energy conserved;
- a field member, an ordinary member and a field common port: each member
  takes its own metric adjoint, the trajectory follows the weighted midpoint
  equations for 600 steps and the three ports' received power sums to zero;
- no member, a wrong role (the common port named as a member included), a
  map of the wrong shape, a non-finite coefficient (its index counted
  through the members' maps) and a member named twice or already connected
  are refused, claiming nothing; the ports then join and the system seals
  and steps.

Field conversions, `tests/field_maps.rs`:
- nonmatching 2-to-3 field oscillators with distinct off-diagonal metrics,
  full-rank and rank-one transfers, compared for 600 steps with independently
  condensed midpoint equations and physical energy;
- metric mismatch refusals for all direct junction and transformer forms,
  mixed ordinary/field maps in both directions, and explicit identity bonds;
- non-finite metric-adjoint refusal without claiming ports, then repair;
- finite extreme-scale work, and atomic refusal when an infinite
  absolute-term power magnitude cannot supply a finite audit allowance.
These are discrete connection oracles, not a field-discretization error bound.

Owner slice, `tests/rigid_body_slice.rs` (`--features slices`): the native
fs-mbd rigid body (`fs_mbd::component::RigidBody`), the first physics owner on
the kernel.
- At one state under force, torque and gravity it passes the conformance
  check, and the work the kernel books from the wire values equals the
  integrator's own `applied_work_j` to 1e-13 of the paired magnitudes.
- At 400 generated states (heights to 500 m, fall speeds to 100 m/s, a
  principal frame rotated against the body frame, free fall, a torque across
  the spin, general wrenches, two step sizes) it passes the conformance check.
- A body released at the origin under gravity follows the parabola for 5,000
  steps to 1e-10, with the stored energy unchanged to 1e-10 of the kinetic
  energy.
- A tumbling body bonded to a three-axis spring follows the monolithic
  midpoint rule for translation over 4,000 steps with a closed ledger.

Owner slice, `tests/viscous_network_slice.rs` (`--features slices`): the
fs-mbd damped rigid-body network (`fs_mbd::component::network::ViscousNetwork`,
the owner's `step_viscous_network_with_loads` behind the contract), three
bodies and two connectors, one with off-centre arms and rotational damping.
- Under forces and gravity at three step sizes it passes the conformance
  check; the energy the kernel books equals the owner's applied force work
  minus its connector heat to 1e-12 of the paired magnitudes; heat outputs,
  stored change and numerical loss are the owner's numbers bit for bit; the
  entropy is `sum(D_e / T_e)`.
- A propose returns the owner's own states, bit for bit, and the mean
  velocities of its displacements.
- Refusals: no bodies, a connector to a missing body or to one body twice, a
  temperature that is not positive, a state of the wrong length, cancellation.
- The damping-heat pairing as a declared graph: connector 0's heat goes 60 %
  to one `HeatCapacity`, 30 % to another and the rest to a named sink,
  connector 1's heat all to the second (a shared receiver). Over 400 steps the
  network's states equal the owner marched alone, bit for bit; each cell's
  account, the sink and the network's dissipation match the fractions of the
  owner's summed connector heat to 1e-11; the ledger closes.
- With gravity and a held force on every body, the boundary work the kernel
  books equals the owner's applied force work to 1e-11, and work in equals
  mechanical change plus heat.
- Built with a damping law (`ViscousNetwork::temperature_dependent`), under
  held forces and gravity: the step conforms and equals the owner's step with
  each connector's form taken from the law at its signal's temperature, bit
  for bit; entropy is generated at the heat ports' temperatures, not the
  damping ones; the network without the law moves differently; a temperature
  outside the law and signals of 0, -5, infinity and NaN are refused.
Five planted defects were each caught: half of a connector's heat leaving,
velocities written to the wrong body, entropy at the other connector's
temperature, a later route replacing an earlier delivery, and routes that
share a receiver audited alone.

Owner slice, `tests/modal_structure_slice.rs` (`--features modal-slices`):
the linear modal structure of fs-nlmodal (`fs_nlmodal::component::ModalStructure`,
fs-modal's exact held-force step behind the contract), three modes with one
loss factor, a one-row strike attachment of mode shapes and a three-row
identity (modal) attachment.
- At a full sample, a half sample and an uncached duration it passes the
  conformance check, and the tangent check at relative 1e-6 with the absolute
  allowance of `tests/elements.rs` (the roundoff of an output at the kit's
  difference step); the energy the kernel books equals the owner's work minus
  its viscous heat to 1e-12 of the exchanged magnitudes; the heat output,
  dissipation and entropy `D / T` are the owner's numbers; the numerical loss
  is below 1e-12 of the exchanged magnitudes plus the stored change.
- A propose returns the owner's own states bit for bit and the mean flows of
  their displacements, and the stored energy is the owner's step total bit
  for bit.
- A free decay of 4,410 samples with the heat routed to a sink: the states
  equal the owner marched alone, bit for bit; the sink, the structure's
  dissipation and its energy loss match the owner's summed heat; the ledger
  closes; each mode's energy follows `exp(-2 zeta omega t)` within the
  within-cycle ripple bound `2 zeta / (1 - zeta)` of a damped oscillator.
- A half-sine strike through a boundary: states bit for bit, boundary work
  equal to the owner's summed work to 1e-11, work equal to heat plus stored
  energy, the ledger closed.
- Refusals: weights of the wrong count or not finite, a port named `heat`, a
  temperature that is not positive, a state of the wrong length, cancellation.

Owner slice, `tests/loaded_modal_slice.rs` (`--features modal-slices`): the
joint loaded modal basis of fs-nlmodal (`fs_nlmodal::loaded::LoadedModalBasis`)
run as one undamped `ModalStructure`, with the loss of the set, `Phi^T C Phi`,
an `elements::Resistance` joined to its modal port by a common-flow junction
and both heats routed to sinks.
- One mode under added mass 0.012 and loss 6.9 1/s: both columns of the
  kernel's step matrix equal the held scheme worked by hand (each loaded mode
  stepped exactly under a held force, the force the resistance's answer to the
  step's mean velocity) to the solve tolerance. The step's determinant is
  `(1 - e) / (1 + e)` with `e = (c' dt / 2) sinc^2(w' dt / 2)`: the energy
  envelope decays by that factor per step (within one step the ratio depends
  on the phase, between its square and 1: both ends are checked as the
  eigenvalues of the kernel step matrix's Gram matrix). The decay rate it gives
  exceeds `held_loss_factors` times `c'` by the factor `atanh(e) / e`, so by
  at most `e^2 / (3 (1 - e^2))`.
- Two coupled modes near a bar's first mode and its tube's, with a mutual
  added mass and a mutual loss: a free decay of 4,410 samples equals the held
  scheme summed by hand to the solve tolerance per step; the resistance's heat
  equals the modes' energy loss to 1e-10 of the start, the sink takes it, the
  undamped modes book none, and the ledger closes.

Owner slice, `tests/thermal_timeline_slice.rs` (`--features slices`): the
fs-thermal explicit enthalpy timeline
(`fs_thermal::component::ThermalTimeline`, the owner's `advance_observed`
behind the contract), two linked cells with a law each.
- At three steps and heat rates it passes the conformance check, with the
  returned temperatures equal to explicit Euler by hand; again with the two
  cells' enthalpies at +1e9 and -1e9 J/kg, where the stored total is small
  and its terms are not.
- `limit_dt_s` is the owner's documented bound `0.8 C / G` and what the owner
  proposes when driven directly; a step of one bound, and one of two and a
  half bounds in three substeps, agree with the owner driven directly with
  held rates to 1e-11 K and 1e-12 of the energy.
- At a clock time of 2^30 s, where an ulp is 0.24 us: a 1 us step (four ulps
  on the clock, one substep the target clips), a 1.1 us step under a 1.1 us
  bound (five ulps, one substep the bound limits) and a 10 us step under a
  0.3 us bound (42 substeps of one ulp, each integrating 0.3 us) all store
  the kernel's `rate * dt` to the rounding of the enthalpy. A step under
  half an ulp is refused.
- A step cut to an event whose time does not round-trip in floating point
  ends on the event; a shorter step does not.
- Refusals: no ports, a cell twice, a port on a missing cell, a rate that is
  not finite, the owner's own refusal (from the step limit and from the
  step), cancellation on entry and cancellation raised inside the owner's
  step, a birth at the end of a step and inside one; without the environment
  port, a film, a reservoir, authored sources that cancel, a film at ambient
  temperature, an event that adds a film, and a state admitted by a timeline
  that has the port.
- Published temperatures are those at the start of the step whatever heat
  the step delivers, also over three owner substeps and after an environment
  port; after a step they are the temperatures that step ended on; a cell
  published twice (in one call or two) and a published cell that is not
  active are refused.
- With the environment port the returned rate is the film's and the authored
  source's by hand, also on a timeline with that port alone, and exactly
  zero while nothing exchanges. In a kernel run with fan changes at 0.4 s
  and 0.7 s, a step ends on each change, the cells match the owner alone to
  1e-10 K and the heat in the named sink is the owner's boundary heat to
  1e-9.
- The damping-heat pairing with both real owners: the fs-mbd network's
  connector heat routed 60 % and 30 % to the two cells and the rest to a
  named sink, the second connector's all to the second cell. The thermal
  bound (10 ms) is the system's clock. Over 300 steps the bodies equal the
  owner network marched alone bit for bit, the cells equal the owner
  timeline fed the same heat by hand to 1e-11 K, and every joule is in the
  cells or the sink to 1e-11.
- Durable snapshots (`ThermalTimeline::restoring`), five tests on a test law
  with a caloric identity. The owner's identity words are equal for a
  timeline built again and differ for another law name, storage basis, slope
  bound or budget limit, another field of a born cell or of a link, and
  another time or kind of a later event; an unnamed law gives none. Its state
  words are read back on a simulation still before its start, uneven covered
  areas included; a state past a separation and a boundary change is read
  back and goes on to the bits of the run never saved, also on a kept
  simulation that still holds the link; fourteen kinds of words are refused.
  The component writes and reads a state on the timeline it keeps, with an
  identity that follows its ports, published cells, step limit, environment
  port and timeline; it does not write a state of another timeline, and does
  not read one with an event due, without a port cell, or with an exchanging
  boundary and no environment port; with no kept timeline, on an unnamed law
  and with an implicit driver it has no codec. Through `System::snapshot` and
  `resume`, a system sealed past a separation is resumed by a second one
  that held the start's link: it then holds the two cells and no link, and
  both go on to the same bits; a timeline that keeps none is refused with
  `NoCodec { node: "cells" }`.
- The implicit driver (`ThermalTimeline::implicit`), backward Euler and
  SDIRK2: a quarter-second step conforms and equals the owner's own solve
  with the rates held (one interval, cells to 1e-11 K), backward Euler also
  by hand to 1e-9 K, and the two schemes differ; on a 400 W/K link, where
  the explicit bound is 10 ms, a half-second step leaves the cells
  `76 / 405` K apart, as one backward-Euler interval does; a constant 3 W
  source owned by the conduction graph, which no cell's boundary shows,
  comes out on the environment port and is stored with the port heat.
- A rejected interval: with a conductivity that rises with temperature and
  two Newton iterations allowed, a half-second step is answered
  `Shorten { 0.25 }`; a half of exactly `minimum_dt_s` is still asked for and
  one under it is the owner's `NotConverged`; with one halving allowed the
  owner alone refuses the half second and so does the Component after the
  quarter second. In the kernel the run takes as many steps as the owner
  alone accepts intervals and as many shortenings as it rejects, and reaches
  the same cells.
- A loose solve (two iterations, 1e-4 J a cell): the receipt's residual is
  the owner's (above 1e-6 J, inside the tolerance), the total is kept to
  rounding and the environment port is zero. Over 20 us on the 0.2 W/K link
  the owner accepts without iterating: the stored change is exactly zero and
  all 8e-5 J delivered is the receipt's numerical loss.
- An implicit run with a fan change at 0.4 s ends a step on it (seven steps
  for the second, as the owner alone), reaches the owner's cells to 1e-10 K
  and puts the film's heat in the room's sink; a step that reaches a birth
  is the owner's refusal. Refusals: an event pending at the start, a cell
  twice, a longest interval of 0, -1, infinity or NaN, a rate that is not
  finite, controls the owner refuses, the law's refusal, cancellation
  raised inside the owner's call.
- The damping-heat pairing on an implicit timeline, both schemes, at steps
  of an eighth of a second (twelve and a half explicit bounds): over twelve
  steps the bodies equal the owner network marched alone bit for bit, the
  cells equal the owner's solve fed the same heat by hand to 1e-10 K, and
  every joule is in the cells or the sink to 1e-11.
- End-of-step temperatures (`returning_end_temperatures`, backward Euler
  only): the step, the receipt and the environment port are those of the
  default timeline and each cell port returns the end temperature to the
  bit; the explicit rule and SDIRK2 refuse the option. Two cells with no
  link inside, the kernel's `Conductance` of 0.2 W/K between their ports,
  against the same cells joined by the owner's link and stepped by the
  owner, 40 steps of 0.5 s: 1.1e-13 K apart (5.7e-14 K before a block's
  differenced columns were resolved; allowed 4.1e-10 K from the two solves'
  tolerances); the same graph with the default mean-temperature
  ports is 5.5e-2 K away, the trapezoidal rule across the wall.
- A contact carried through the ports
  (`ConductionObserver::with_links_carried_outside`, fs-thermal): two cells
  in contact, each also losing heat to the room through a film on its
  exposed face; the conduction graph names the contact as carried outside
  and the kernel's `Conductance` of 0.2 W/K joins the two end-temperature
  ports. Against the owner's backward-Euler solve with the contact's
  conduction inside, 40 steps of 0.5 s: 1.7e-13 K apart (5.7e-14 K before a
  block's differenced columns were resolved; allowed 4.1e-10 K);
  16.758 J to the room in the owner's audit and in the kernel's sink, to
  1e-9 of it; the same cells not in contact, whose films act on the whole
  faces, end 0.19 K away, which is what losing the covered area would be.
- Entropy (`an_implicit_timeline_books_its_owners_conduction_entropy`): the
  quarter-second step with a 280 K film on cell 10 books, under both schemes,
  the owner's own `conduction_generation_j_k` to 1e-12 of it. Under backward
  Euler that is the link's `G (T_a - T_b)^2 / (T_a T_b)` plus the film's
  `q (T_0 - T_a) / (T_a T_0)` at the end temperatures. The explicit rule books
  exactly zero.
- Declared exchange
  (`an_implicit_timeline_books_a_declared_exchanges_production_and_none_below_zero`):
  a coupled source moving 0.05 W from cell 10 to cell 20, declared closed with
  a zero tolerance, books the owner's conduction plus exchange production, to
  1e-12 of it, under both schemes. 1e-3 W moved the other way inside a 1 W
  tolerance leaves a negative production in the owner's audit and books the
  conduction alone.
Thirty-seven planted defects in the explicit component were each caught, and
twenty-three in the implicit driver; none was planted for the end-of-step
option or for the contact carried outside. Five more were planted in the
entropy wiring, and each failed the entropy test: the implicit step booking
zero, its storage change, or its generation plus its reservoir input; the
explicit rule booking heat; the receipt ignoring the solve. Two more, in the
declared exchange's booking (the conduction alone; a negative production
booked), each failed the declared-exchange test.

Against the hand-written coupling,
`crates/fs-couple/tests/legacy/rigid_damping_heat/graph.rs` (run with
`cargo test -p fs-couple --test rigid_damping_heat`): the same graph, built
from that module's own test fixtures (four bodies, two connectors, both
routed to both cells, nonlinear calorics, cells on a bed), against
`fs_couple::rigid_damping_heat::advance_candidate`. In four cases (with and
without an authored source, a contact and a source change at 0.35 s; targets
on and off the thermal step) the kernel takes as many steps as the module
takes thermal substeps, the bodies are equal bit for bit, the two thermal
clocks are equal, the cells agree to 1e-11 K, and the module's report
(dissipation, unallocated loss, boundary heat, heat stored in the cells,
mechanical energy change) equals the kernel's ledger to 1e-11 of the
dissipated energy. Five planted defects in the graph side were each caught.

The same file compares the module's temperature feedback
(`rigid_damping_heat::temperature::advance_candidate`) with a graph in which
the timeline publishes cell 0's start-of-step temperature as a signal and the
network built by `ViscousNetwork::temperature_dependent` reads it through the
module's own synthetic law (`c = 0.7 (1 + 2 (T - 300 K))`, the cell warming by
0.7 K and 1.7 K). For targets of 0.5 s and 0.2 s at a 1/32 s thermal step: as
many kernel steps as thermal substeps, the cell to 1e-11 K, the bodies'
velocities to 1e-9 m/s, and dissipation, stored heat and mechanical energy to
1e-11 of the dissipated energy. The kernel solves the signal to its relative
tolerance of 1e-12, so the law sees the module's temperature to that; in this
fixture the bodies come out equal bit for bit, which the kernel does not
promise. Not compared: the module refuses a whole advance when its law has no
form at the final temperature, where the graph takes that step and refuses
the next; a run with thermal events; two connectors reading one cell.

The same file compares the module's implicit entry points with graphs on an
implicit timeline, on the module's own stiff fixture (two cells on a
1000 W/K link, which the explicit rule cannot step within its substep
allowance), for backward Euler and SDIRK2:
- `advance_implicit_candidate`, both connectors' heat into cell 0, two
  quarter-second intervals: the same intervals, the bodies bit for bit, the
  two thermal clocks equal, the cells to 1e-9 K, and dissipation, stored
  heat and mechanical energy against the ledger;
- `temperature::advance_implicit_candidate`, cell 0's temperature setting
  the connector's damping, eight sixteenth-second intervals: the cells to
  1e-9 K, the bodies' velocities to 1e-8 m/s, the accounts to the solve's
  closure;
- through rejected intervals (calorics that curve, two Newton iterations
  allowed): the module halves an interval and re-runs its mechanics over the
  half; the graph's timeline answers `Shorten` and the kernel re-steps the
  network with it. With the heat released at a fixed temperature, as in the
  module, the graph takes the module's accepted intervals one for one, one
  shortening for every rejection (20 intervals and 64 rejections for
  backward Euler, 15 and 40 for SDIRK2), and the bodies are equal bit for
  bit. With the heat released at the receivers' temperature the kernel
  solves the network and the cells together and the timeline can answer on
  an iterate: the run is finer there (23 and 76, 16 and 44) and agrees with
  the module's to 6e-6 K.
Bit-for-bit bodies rest on the kernel's negotiated step being the thermal
clock's `end - start`, which holds on every step of these fixtures and is
not a property of every grid.

Owner slice, `tests/relaxation_network_slice.rs` (`--features slices`): the
fs-mbd network with Maxwell memory
(`fs_mbd::component::relaxation::RelaxationNetwork`, the owner's
`step_relaxation_network_with_loads` behind the contract; bodies and the
complete elastic memory are one state, and the memory is stored energy). The
list of what it shows is in `crates/fs-mbd/CONTRACT.md`.

Against the hand-written Maxwell-memory coupling,
`crates/fs-couple/tests/legacy/rigid_damping_heat/relaxation/graph.rs`: the module's
own pair (Kelvin viscosity 0.4, one Maxwell branch preloaded to 0.6 sqrt(J),
closing at 2 m/s, no actuator), its connector's heat routed 70 % and 20 % to
two timeline cells and the rest to a sink.
- `relaxation::advance_candidate`, explicit thermal steps, targets 0.25 s and
  0.3 s: as many kernel steps as thermal substeps, the bodies and the memory
  bit for bit, the cells to 1e-11 K, and dissipation, unallocated loss, heat
  stored and the change of mechanical and memory energy against the ledger to
  1e-11 of the dissipated energy.
- `relaxation::advance_implicit_candidate`, backward Euler and SDIRK2, with
  two Newton iterations allowed so that the solve rejects intervals: the
  module's accepted intervals one for one, one kernel shortening for every
  rejection, the bodies and the memory bit for bit, the cells to 1e-10 K.
Each of the two is run with the held effort zero and no actuator in the
graph, and with 0.7 N and -0.7 N held on the connector as a boundary on the
Component's `drive0` port; the module's signed actuator work is then the
energy that entered through that boundary, to 1e-11 of the dissipated energy.

Against the hand-written driven couplings,
`crates/fs-couple/tests/legacy/rigid_damping_heat/driven_graph.rs`: the same
Component with no Maxwell branch and an actuator on every connector.
- `driven::network::advance_candidate` and `advance_implicit_candidate`
  (backward Euler and SDIRK2): the module's cycle of three bodies and three
  connectors plus a free fourth body, held forces `(F, -F / 2, 2F)` for
  `F = -0.35` and `0.4` N.
- `driven::advance_candidate` and `advance_implicit_candidate`: the module's
  pair with 0.7 N s/m, held forces of -1.4 and 1.4 N.
The module calls fs-mbd's driven-pair and driven-network owners and the
Component the relaxation-network owner, which share connector geometry,
impulse application and rigid drift and differ in how they assemble and
solve the source. So the bodies agree to rounding and not by bits: measured
1.7e-15 in position and momentum (allowed 1e-13), the cells equal to the bit
(allowed 1e-12 K), and dissipation, unallocated loss, heat stored, mechanical
energy and each actuator's signed work against the ledger to 7.3e-15 of the
dissipation plus the work (allowed 1e-12). The thermal intervals are the
module's one for one.
Five planted defects in the graph side of these comparisons were each caught.

With temperature feedback, same file: `driven::advance_temperature_candidate`
and `driven::network::advance_temperature_candidate`, each with its implicit
twin. The timeline publishes cell 0's start-of-step temperature and the
Component's damping law (the module's own synthetic law) reads it on one
signal per connector. The module's pair (0.6 N held, c = 0.7 (1 + 2 (T -
300))) and its cycle (forces 0.2, -0.1, 0.4 N, c = 0.7 (1 + (T - 300) / 2)),
sixteen explicit thermal steps and two implicit steps of a quarter of a
second in both schemes; cell 0 warms by 0.39 to 1.37 K. Measured agreement:
5.9e-15 for the bodies, one ulp of 300 K for the cells, 5.2e-15 for the
energy accounts and each actuator's work (allowed 1e-12, 1e-11 K, 1e-12).
The cycle is run to 0.25 s and not to the 0.3 s of the module's own test:
sixteen steps of 0.3 / 16 s do not sum to 0.3 s in floating point, the
module's thermal clock lands on its target and the kernel's ends one ulp
short.

With calibrated kinetics, `relaxation/graph.rs`:
`relaxation::temperature::advance_candidate` and its implicit twin on the
module's own fixture (the pair preloaded to 2 sqrt(J), 0.7 N held, an
Arrhenius shift about 300 K read at cell 0). The Component follows a shift
law that calls the module's calibration; both sides call the owner's
held-kinetics entry. Bodies and memory equal to the bit in the three runs
(allowed 1e-12), with the shift factor ending at 0.77 (explicit) and 0.89
(implicit).

The allowances of these two comparisons hold the kernel to what it does
today, delivering a start-of-step signal exactly. Its contract is the solve
tolerance; if it ever uses that here, they become 1e-9. Four planted defects
in the graph side of these two comparisons were each caught.

Every entry point of `fs_couple::rigid_damping_heat` that advances a state
now agrees with a declared graph: sixteen of them (explicit and implicit,
each for frozen damping, temperature-dependent damping, Maxwell memory,
Maxwell memory with calibrated kinetics, the driven pair, the driven network,
and the last two with temperature-dependent damping). The comparisons hold
no body torque and supply no thermal observer. The module takes both: a
`Wrench` per body, whose torque the two network Components now take on an
opt-in `torque<i>` port (crates/fs-mbd/CONTRACT.md; not yet compared with
the module), and a caller's `DepositionObserver` (sources, conductances,
heat transfers, a transition time and a step limit on the thermal clock),
which the timeline Component does not take. Not reproduced either:
`temperature::shift::CalibratedKelvinDamping`, which is a law and not an
entry point, has no adapter to the Component's `DampingLaw`; the module's
refusal behaviour and continuation from a candidate are the module's own and
are not compared. Nothing in the module is deleted.

Gas piston, `crates/fs-couple/tests/gas_piston_graph.rs` (run with
`cargo test -p fs-couple --test gas_piston_graph`): the hand-written
`fs_couple::cavity_phs::piston` against a graph with no pairing code. The
wall is a `Mass`; `common_flow_ratios` joins to its velocity the fs-gas
cavity (`fs_gas::component::GasCavity`) at the face area, a constant
exterior pressure at minus that area, and a constant force at one; heat is a
boundary on the cavity's heat port. Both sides use the gas, controls and
0.05 ms clock of the module's own tests (a two-species mixture with a
temperature-dependent heat capacity, a litre at 400 K and 1 bar, a 1 kg wall
of 0.01 m²) for 48 steps.
- `GasPiston::step_discrete_gradient`, three runs: released into vacuum;
  pushed inwards at 0.2 m/s against 1 bar by 0.5 N while the gas takes 20 J
  a step for eight steps; thrown outwards at 0.2 m/s against 1 bar.
- `DoubleActingGasPiston::step_discrete_gradient`, a second gas of another
  composition on the other face at minus its area, two runs: equal faces,
  and a smaller opposing face (0.008 m²) with a 1.5 N force; the first gas
  is heated as above.
- The module's material wall (`piston::material`), whose mass is the
  resolved body's: the fixture's plate (0.1 m square, 10 mm, densities 1000
  and 8000 kg/m³, so 0.1 and 0.8 kg) as a free fs-mbd `RigidBody` in three
  dimensions, the gas on its face through a 1 x 3 `power_map` of the area,
  its torque port held at zero, against `GasPiston` with that mass and the
  discrete-gradient step. Two runs, the gas heated as above. Wall velocity
  within 9e-12 m/s, volume 2e-16 m³, temperature 4e-11 K; the body stays on
  the axis to 9e-19 (momentum, position, spin) with nothing constraining
  it; ledger 4e-14 J. `MaterialGasPiston` itself offers only the midpoint
  step, which the fs-gas cavity does not have.
Measured at every step: wall velocity within 1e-11 m/s (allowed 1e-10, ten
times the module's impulse tolerance on the unit mass), gas volumes within
2e-16 m³ (allowed 1e-15), gas temperatures within 4e-11 K (allowed 1e-9).
The sources' accounts equal the module's exterior-pressure and force work to
2e-11 J of 1 J and 8e-15 J, and the ledger closes to 2e-14 J.
The two-gas runs give the module a volume bracket of 0.9 to 1.1 litres. With
the 0.7 to 1.5 litres of the module's own tests it refuses the first step:
it tries both ends of the bracket, and at the far end the opposing gas would
be halved, which the gas owner refuses for this mixture. The graph has no
bracket. Not reproduced by a graph: the module's time-resolution admission
(the graph accepts any step), its fixed support, impulses and impacts, and
its real-gas, reactive, nozzle-exchange, thermal-body and electrical
variants. Nothing in the module is deleted.

Friction heat, `crates/fs-couple/tests/friction_heat/graph.rs` (run with
`cargo test -p fs-couple --test friction_heat graph::`): the hand-written
`fs_couple::friction_heat`, which solves the temperature of a massless
plane between two cells inside one backward-Euler thermal step with the
friction work generated on the plane, against a graph with no pairing code.
The cells are the module tests' own network on a backward-Euler
`ThermalTimeline` with end-of-step ports and the contact link named as
carried outside (fs-thermal); the plane is `elements::HeatNode` with the
route's two half conductances; the friction is fs-tribo's
`tangential::component::SlidingInterface` under the tests' prescribed slip,
its heat port on the plane's `source`. The criterion is the module's own
verdict: the endpoint of every kernel step (each cell's energy increment and
the temperature the interface was given) is handed to
`FrictionHeatStep::evaluate` with the module tests' law callback and must be
accepted at the tests' 1e-10 J. Acceptance is the evidence: the module forms
each cell's residual from its own authored conductances, its own route
conductances and its own call of the law. The tangential state and traction
of that call equal the graph's exactly, which shows that the interface
component drives the law with the tests' input and history; the cells the
module publishes are the handed increments added back and say nothing more.
The six base cases:
- one cell a side with unequal half cells, 10 W for 0.2 s: accepted, the
  module's largest residual at the graph's endpoint 9.1e-14 J (2.2e-13 J
  before a block's differenced columns were resolved, 2026-11-14); the plane's
  temperature and the cells' energies equal the tests' own LU solve in
  temperature coordinates to 1e-12; the share of the first body is the
  tests' closed form, on the module's report and on the graph's own wires;
  the ledger has 2 J dissipated by the interface, 2 J stored in the cells,
  and the interface's entropy is the heat over the plane's temperature;
- no slip, the first cell 10 K warmer: accepted (1.6e-13 J; 6.8e-14 J
  before), the exchange is
  the tests' closed form and nothing is dissipated;
- a coefficient `0.8 - 0.01 (T - 300 K)` with films to a 300 K room, 80
  steps of 0.25 s: every step accepted (4.3e-12 J at worst; 4.1e-12 J
  before); the plane is
  2.8e-12 K at worst from the module's own run from its LU endpoints, and
  1.5e-12 K from the analytic steady temperature at the end; 153.562 J are
  in the room's sink, as the module's reports have them;
- a law that refuses during the step, and cancellation requested inside
  the step by the interface's own law once the solve has moved the plane:
  the kernel's step fails and cells, clock and interface are as they were;
  the same graph with the timeline's default mean-temperature ports is
  taken by the kernel and refused by the module (`EnergyBalance(-0.276)`
  J): end-of-step ports are what make the two one discrete rule;
- one owner for the contact: the conduction binding refuses a contact both
  conducted and carried outside, and the implicit solve refuses one left
  unsaid;
- 16, 32 and 64 cells a side, 512 steps: every step accepted (6.5e-14 J at
  worst; 1.7e-13 J before); beside each grid the module's own run from its
  tests' LU endpoints, from which the graph's plane is 8.0e-13 K at worst
  at any step (2.4e-12 K before; asserted to 1e-11 K) and each body's total
  heat 3.0e-14 J at worst (2.9e-14 J before; asserted to 1e-12 J); the
  partition and plane-temperature errors against
  the half-space solution (5.20833182e-3 and 4.31961554e-3 at 16 cells,
  3.25520296e-4 and 1.38133321e-4 at 64) fall with the grid as the module's
  test requires.
The test runs at the kernel's default solve settings. When a loop was closed
only to 1e-12 of the largest temperature in the block (3e-10 K near 300 K),
the module refused the third case's ninth step with 2.4e-10 J at those
settings, and the test set a tolerance of its own; since the kernel also
judges a temperature's residual by the heat it drives (above, "A residual is
judged a second time"), it needs none. The wire between the plane and the
interface is open by 1.1e-13 K at worst. Not reproduced by a graph: conductances
that follow the materials' conductivities at the endpoint (the plane's are
constants, as the six cases' conductivities are), a contact resistance split
about the plane, the module's entropy audit, and its area, moving-patch,
third-body, rigid and other variants. Nothing in the module is deleted.

Two source-backed cavities exchanging gas, `crates/fs-gas/tests/component.rs`
(`cargo test -p fs-gas --features component --test component`): the fs-gas
cavity with a stream port (`GasCavity::streaming`, crates/fs-gas/CONTRACT.md)
on each side of the kernel's `Restriction`, against the gas owner's own pair
transfer. Masses and species masses equal to the bit over 40 steps,
temperatures within 5.7e-13 K. This is the first physics owner on a stream
port; the cavity is added without a tangent and the kernel differences it.

Other domains, `tests/elements.rs`: the gyrator and one element of each law
in the rotational and electrical domains pass the conformance and tangent
checks at 100 generated points, and the gyrator's two port powers cancel to
rounding; a DC motor of seven elements (supply, inductance, resistance,
gyrator, rotor, bearing, one thermal case that takes both heats) matches the
monolithic midpoint rule for its two coupled equations over 20,000 steps to
1e-9, the supply's account equals stored magnetic and kinetic energy plus
heat, the gyrator's account is zero, and a rotor port joined into the current
loop is refused. Two planted gyrator defects were each caught (both efforts
with the same sign; the transpose with the wrong sign).

Streams, `tests/elements.rs`: `elements::GasVolume` (a well-mixed volume of
perfect gas with constant heat capacities: a stream port, a heat port and,
with `with_wall`, a pressure-volume port; trapezoidal, energy-exact for any
step) and `elements::Restriction` (mass rate `conductance (p_a - p_b)`, the
enthalpy and carried quantities those of the side the mass leaves). Both
supply measured mass: the volume's validated inventory and the restriction's
zero storage, so their mass balances are checked before publication.
- Both pass the conformance and tangent checks at 100 generated points, the
  restriction on both sides of zero flow; at equal pressures nothing crosses
  (the mandatory check only: the port map has a kink there). The kit's
  tangent check takes one norm over all outputs, which pins the small rows
  of a mixed-unit port only loosely, so every entry of both tangents is also
  compared on its own with a central difference (1e-6 of the entry). The
  volume's returned pressure, temperature, enthalpy and fractions are
  checked against the means of their values before and after the step.
  Steps that would empty a volume, overdraw a carried quantity or take more
  energy than it holds are shortened with the candidate untouched.
- Two vessels (5 bar at 350 K in 10 l, 1 bar at 290 K in 30 l, each holding a
  quantity the other lacks) through one restriction against an independent
  RK4 integration of the same equations: errors 1.35e-6 and 3.37e-7 at 100
  and 200 steps, a ratio of 4.0. While mass leaves the first vessel it keeps
  its composition to 1e-11 and receives none of the other's quantity, which
  no element states. Mass, each carried quantity and energy are conserved
  (4e-15 relative; ledger residual zero). Left alone the pressures meet at
  2 bar, `(p_a V_a + p_b V_b) / (V_a + V_b)`, and the temperatures do not.
- Reversal: the vessels start level; the first is heated for a second, then
  the second at five times the power. The flow turns once. Until the turn
  the first vessel holds none of the second's quantity; afterwards it does.
  Mass and both carried quantities are conserved through the turn and the
  energy rises by the 4.8 kJ supplied. Mass and energy are second order
  against RK4 (6.5e-8 and 1.5e-8, ratio 4.3). The carried quantities' error
  (1.2e-7 and 6.4e-8) is bounded by 5e-3 of the step squared and has no
  clean ratio: the step in which the supplier changes costs up to its length
  squared, depending on where in it the turn falls.
- Throttle entropy. A `Restriction::for_gas` reports `dt |m_dot|` times its
  gas's `GasVolume::specific_entropy_change_j_kg_k`, taken from the upstream
  `(h, p)` to the downstream pressure at the same `h`. The test judges it by
  the closed form `c_p ln(T2/T1) - R ln(p2/p1)`, computed in the test, with
  an allowance of 32 eps of `(want + dt |m_dot| R)`. The allowance's
  derivation is in the test.
  - Over 200 generated points, every tenth a drop of one ulp, it is within
    0.0401 of that allowance and positive at every point.
  - Its port map equals `Restriction::new`'s to the bit. With `a` and `b`
    swapped the entropy is the same to the bit.
  - At equal pressures it is +0.
  - While mass crosses, a pressure that is not positive on a side it reads
    is refused. `new` and `for_stream` report no entropy.
- Throttle entropy in the reversal run (600 steps of 5 ms, 293 one way and
  307 back), the plain and the reporting restriction in one binary:
  - The step reports, both vessels' contents, the clock and every port's
    last input and output are the same to the bit at every step.
  - The ledgers' Debug texts are equal once the pipe's entropy account is
    masked. That account is 2.390423e-2 J/K.
  - Each step's report is within 0.0113 of the allowance above.
  - The account equals, to the bit, a compensated sum of the element's
    reports at the inputs the system accepted.
  - Two vessels level at equal pressure: nothing crosses, and the account
    stays +0.
- A restriction that reports its throttle entropy with the sign turned:
  - The conformance kit refuses it (`NegativeProduction`).
  - In a system, after ten honest steps, the step's audit refuses it
    (`SystemError::Contract`, `NegativeProduction`, node `pipe`). That
    comes before the second-law sum. States, clock, port values and ledger
    are unchanged to the bit, and an honest step afterwards is taken.
- Three vessels in a line at 4, 2.5 and 1 bar, the first holding a quantity:
  it reaches the third only through the second, which receives and gives in
  one step. Against RK4 the quantity's error is 1.6e-7 and 4.1e-8 at 100 and
  200 steps and mass and energy 8.8e-7 and 2.2e-7, ratios of 4.0. A storage
  returning its fractions at the start or end of the step would be first
  order here; the two-vessel cases cannot tell, because a donor's
  composition does not change within a step.
- No step of these runs is shortened, the one in which the flow turns
  included.
- A vessel filled from a line through a prescribed stream boundary ends at
  `(m0 c_v T0 + m_in c_p T_line) / ((m0 + m_in) c_v)`, 355 K from a 300 K
  line, to 1e-12.
- A 1 kg wall on a litre of gas at 1 bar (the outside pressure a constant
  force, all on one ratio junction) rings at the adiabatic gas-spring
  frequency `w² = ratio p A² / (M V)`: the period is within 1e-6 of the
  trapezoidal rule's for that frequency and 1.2e-5 from `2 pi / w`;
  `p V^ratio` holds to 8e-15.
- Three vessels in a line: the middle one is the node two restrictions meet
  at; mass and energy are conserved and the pressures meet at their mean.
  Refused, claiming nothing: a common flow to two storages, two storages or
  two restrictions joined directly, streams carrying different numbers of
  quantities, a stream with a heat port, a ratio, a transformer and a heat
  route on a stream port.

Mass accounting, `tests/mass_accounting.rs`:

- two `GasVolume`s and a `Restriction`, heated in turn for a flow reversal,
  conserve their measured total over 300 steps; per-node inventory changes,
  port transfers, audit coverage and connection defects agree;
- a prescribed filling boundary books the independently known added mass;
  restoring a checkpoint and replaying restores every mass account exactly;
- an energy-neutral fixture creates mass in the last node: the named
  `MassBalance` refusal publishes no node, clock or ledger change, and a
  corrected retry succeeds with only the accepted audits counted;
- an owner that supplies no inventory still records boundary transport but
  has no mass account or system mass-conservation verdict;
- a finite balanced transfer of `1e308` kg succeeds with zero mass defect and
  residual even though the unscaled sum of its absolute transfers overflows.

Entropy accounting, `tests/entropy_accounting.rs`, a hot and a cold body
(backward Euler, entropy `C ln(T / T_ref)`) through a stateless wall, all
three declaring; every step is asserted taken whole, with no shortening:

- over 200 steps every member has 200 audited steps and an absolute defect
  below 1e-10 J/K, and its port, generation, change and signed defect
  accounts agree; each body's booked change is its measured
  `C ln(T_end / T_start)`; the members' generation sums to the ledger's
  `entropy_generated_j_per_k`, and the pair's stored change equals it;
- undeclared, the same pair has no entropy account and a bitwise equal
  trajectory and generation: the audit moves no state;
- a hot body whose entropy is skewed by 1e-6 is refused at its first step as
  `Contract` with `EntropyBalance`, naming it, and nothing is published;
- a cold body that takes `ln` of the quotient passes its first step from its
  reference and is refused at the second;
- a hot body that books 4e-14 J/K too little at every step, inside its
  allowance, ends with a negative signed defect account and a positive
  absolute one of about 200 times that;
- a body whose measurement refuses is named as `Component` in the
  `Measure` phase, and nothing is published;
- a checkpoint rewinds the entropy accounts and a replay restores them
  exactly; a snapshot resumed into a fresh pair carries them, and both
  continue to the same ledger;
- the shipped `Conductance` in the wall's place declares `Some(0.0)` and,
  from 400 K and from 3000 K against 300 K, gives the test wall's states and
  entropy accounts to the bit over 200 steps: 200 audited steps, no stored
  change, positive generation, and port entropy that cancels it within
  1e-12 of it.

Shared receivers, `tests/elements.rs`: two heaters routed into one body (one
of them also into a second body and a named sink) deliver the sum of their
shares, the pair is audited in the first route's row, and the ledger closes;
a receiver named twice in one call, a source used by two routes, and a
receiver that is also in a junction (either order) are refused by name, and
the refused calls change nothing.

The load line, `tests/load_line.rs`, with `elements::Slider`: ideal dry
friction between two ends under a constant normal force, with constant static
and kinetic coefficients, its loss leaving through a heat port. It is the
element the tier is exercised with and nothing more: friction that follows a
material, a surface card, a temperature or a slip history stays in `fs-tribo`.

- The slider's law of a step. The elements here use the midpoint rule, so a
  mass returns its mean velocity over the step. "Stuck means the mean slip is
  zero" would leave a thrown mass flipping the sign of its velocity for ever
  with nothing dissipated, and "stuck means the end slip is zero" alone gives
  a negative loss in a step whose drive reverses the body. The slider keeps
  the slip `s0` its last step ended with and, on the line
  `mean = r - c e`, takes one of seven rows (`elements::SliderRow`, tabulated
  in `src/elements/slider.rs`): stuck and held (`e = r / c`), stuck and broken
  away (`mu_k N` along `r`); sliding and arrested at exactly zero end slip
  (`e = (r - s0/2) / c`, when that lies between zero and `mu_k N` along the
  slip), sliding on, and three rows for a step that reverses, of which
  turned holds the mean slip at zero (`e = r / c`, inside `mu_k N`) and so
  ends stuck, as held does: a hold is stuck from whichever state it came
  (before, turned ended sliding, and a belt caught on the boundary of a
  step stayed turned at every step after). The loss
  `e * mean` is not negative in any row, the force is inside `mu_s N` stuck and
  `mu_k N` sliding, and in the sliding state the force does not decrease with
  `r`. One row trades the friction of a step for first-law exactness: where
  the drive alone reverses the body within the step, no friction acts in that
  step, because the force that would end the step at zero slip has a
  negative loss there, and no force keeps the row continuous with its
  neighbours. A force that is a quotient `r / c` is kept inside the limit
  its row was chosen by (the quotient could exceed it by a unit in the last
  place). The slip is the flow `a` receives less the flow `b` receives,
  with the signs the ports receive them in: a junction member receives minus
  its setter's flow, so a slider joined to a mass moving at `v` starts with a
  slip of `-v`, and nothing checks a starting slip against the ports; a
  pair that starts stuck (`Slip::stuck()`) needs none.
- What a port tells the slider. What an end receives is a step mean: a
  mass returns its mean velocity, and a prescribed end gets the mean of
  what is prescribed (`FlowSource` delivers `Waveform::mean`, a boundary is
  handed `(t0, dt)`). Friction needs the slip a step starts and ends with.
  For an end that answers, the end velocity is twice the mean less the
  start. For a prescribed end a mean does not say where the side ends its
  step, and the slider takes the mean for the end: the side is held over
  the step and jumps between steps. The line says which end is which: an
  end whose row of the matrix is zero is given. The slider keeps the slip
  its last step ended with, which goes to twice its mean less its start
  within a step either way, and the flow each given end was given; a given
  end that comes back with another value has jumped, and the step starts
  with that much more slip. It also keeps whether it ended stuck. Stuck
  with such a start slip, it gives the force that ends the step at zero
  slip, up to `mu_s N`, where that has no negative loss; else it holds on
  the mean, up to `mu_s N`; else it breaks away. (The first version looked
  for no jump: after a prescribed side jumped its record alternated in
  sign for ever, friction acted every other step, and every audit passed.
  Found by the independent review of 6ad46a38fd.)
- A prescribed side that is told its span. `Slider::told(SliderEnd)` gives
  an end the signal input `span_a` or `span_b` of width two: the
  velocities of what moves that end at the start and at the end of the
  step, as the mover returns them. `FlowSource::publishing` puts them
  out; a boundary can prescribe them. Nothing is inferred (three ways of
  inferring the end of a step from the means were emulated, and each
  fails a basic case). The sign is the mover's: every connection hands a
  port minus the flow of what it is joined to, so the flow the port
  receives is the mean with the other sign. The end's step then differs
  from the midpoint relation's by `2 mean - from - end` (in the port's
  sign; `from` is where its last step ended, the span's start in the
  first step), and the step's start slip is the last end slip and that,
  for `a` less that for `b`. A held end is the case `end = mean`, where
  this is the jump, in the arithmetic it had. `Slider::on_line` and its
  rows are untouched: only the start slip they are handed differs. With
  both ends given the slip ends where the two ends end.
  Refused by name (`ElementError::Span`, nothing published): a span on an
  end that answers the slider's force (its row of the matrix is not
  zero), at once; and a span of the port's sign, at its third step. A
  span on one side of zero belongs to a flow on the other. One step
  cannot say that a span on the flow's side is wrong: a side that reverses
  twice within a step has both its ends on one side and its mean on the
  other. And the slider may not decide by a step's inputs at all, because
  the kernel also proposes at inputs it has perturbed to difference it,
  and the perturbed flow of a slow side changes sign. So the slider
  counts, per told end, the steps taken in a row with the span on the
  flow's side of zero, and refuses the step after two of them. Two such
  steps in a row are a span given in the port's sign, or a motion the
  step does not resolve; both are refused, and two steps with wrong
  friction have been taken before the refusal. (b47fe4f2a8 decided by one
  step's inputs. It refused a side at 1e-9 m/s beside a body at 0.2 m/s at
  its first step, and a side that reversed twice in a step. Found by its
  review.)
- Not covered: a prescribed side that accelerates and is not told its
  span. Held is the truth for a
  side at a constant speed and for one that does jump. For a side that
  accelerates at `a` it is off by `a dt / 2` at every step, and what
  follows does not get better with a smaller step (2 kg, `mu_s N` = 9.8 N,
  `mu_k N` = 5.88 N, a side at 4 m/s² that needs 8 N):
  - A side that leaves rest at `t = 0`, through a `FlowSource` or through a
    boundary that gives the same means. Its first mean is `a dt / 2` with
    the mass at rest, a slip the slider is not told, and its record is off
    by that from the first step on. The first step gives 8 N. Then, by the
    sign of a rounding residue, the slider either goes on holding on the
    mean for a while, or sees a jump twice what there is, gives nothing in
    that step and breaks away in the next. It breaks away within 27 steps
    at every step length tried: at step index 2 at 1, 0.5 and 0.25 ms
    (powers of two apart, the same arithmetic rescaled), and at the six
    others at step index [4, 4, 4, 4, 8, 26] at 0.9, 0.7, 0.6, 0.3, 1.5
    and 1/3 ms (ts2 chain fsv, 2026-11-14; the 26 first on chain fsz).
    While the component was first asked at the released flows of the
    step before, they were [4, 2, 16, 16, 6, 6] (chain fsx). At 0.2 s
    the mass has between 0.72 and 0.76 of the side's 0.8 m/s, and between
    0.11 and 0.14 J of heat are booked for a pair that in truth never
    slips.
    Asserted: 8 N or nothing until it breaks away, at the six step lengths
    that are not powers of two apart, the motion and the heat at each, and
    then the step index it breaks away at to 26 at most, a bound read off
    chain fsz, not derived (it hangs on the sign of a rounding residue);
    the three forces and the end speed to 1e-12 at the other three.
    (74529fda7f said "at the third step at every step length", having run
    only those three. Found wrong by its review.)
  - A side that rests one step and then accelerates. Every change of the
    mean is then the jump the slider takes it for: the mass is carried with
    the right force after the first step (which gets half), but one mean
    behind, `a dt / 2`, and every step books the catch of its jump as heat,
    `m a² dt² / 2`: 3.172e-3, 1.593e-3 and 7.98e-4 J in 0.2 s, first order
    in the step, where a stuck pair loses nothing.
  - A side that swings from rest at 2 Hz and needs at most 6.32 N. The
    slider survives on holding on the mean until, at one of the side's
    reversals, it breaks away. Whether it does, and at which reversal,
    hangs on the last digit of the side's means: of the side as it is and
    39 more with every mean left or moved by a unit in the last place, 24
    break away within 4 s and 16 never do (1 ms; through a `FlowSource`
    the side as it is broke away at 0.751 s, and at 0.251 s on
    ef1631e0d2). Asserted: some of the 40 break away, each in the step of
    a reversal or the one after.
  - With a force on the mass that exceeds `m a` along the acceleration the
    slider falls back on holding on the mean, which for the step means of
    a ramp happens to be the exact answer (-7 N under a push of 15 N in a
    scratch run of the kernel, not a test). That is no claim: nothing
    tells a caller which of these cases a run is in.
  The slider cannot refuse such a side by name: it cannot tell a side that
  accelerates from one that jumps. A prescribed side that does not stay at
  one speed is to be told its span. The record of an end that answers is
  as good as the midpoint relation is for what is attached there.
- Not covered: a change of the side's acceleration inside a step, with
  the span told. The slider is told the step's two ends and its mean, no
  more. In the step of the change, where the side speeds up more, the
  force that would end the step on the side has a negative loss (the mass
  would move faster than the side on the mean while being pulled along),
  and the slider holds on the mean; the mass ends that step behind, and
  the next step makes it up with a force above the needed one by `m` times
  a quarter of the change of acceleration at most, for one step. After a
  rest of 2.5 steps and then 3.5 m/s² (7 N needed): 1.75 N, 8.75 N, then
  7 N, stuck throughout, the mass a quarter of a step's gain behind after
  the first of them and on the side from the second. With 4 m/s² (8 N of
  the 9.8 N held): 2 N, and the 10 N the next step wants break the pair
  away, though the side never asks for more than 8 N. A change of
  acceleration at the end of a step is followed exactly.
  Changes in consecutive steps add up: a step held on the mean leaves a
  slip in the record, and a change in the very next step adds its own
  before that one is made up. With changes of 1.48, 1.88 and 3.02 m/s² in
  three consecutive steps the force is 2.917 N above the 4 N needed, 0.483
  of `m` times the largest change. Over 2000 drawn sides with a change in
  each of two to five consecutive steps it is at most 0.424 of `m` times
  the largest change and 0.248 of `m` times their sum (measured; each
  change adds at most its own quarter, asserted for those draws and not
  proved). A side with changes in four consecutive steps that needs 6 N
  of the 9.8 N held breaks the pair away. The quarter of one change holds
  for changes that are at least three steps apart. (b47fe4f2a8 gave the
  quarter as the bound for any side. Found wrong by its review.)
- Not covered: a starting slip that contradicts the ports. A mass at rest
  on a belt that already moves at 0.2 m/s slips at 0.2 m/s; declared
  `Slip::stuck()`, the slider's record is off by that slip for ever and
  alternates in sign, and its force is a difference of flows of 0.2 m/s
  over `c`, which moves in steps of 1.1e-13 N. Up to ef1631e0d2 such a run
  went on with its wrong record (6000 steps of 1 ms on a spring of
  2000 N/m). Since the solve judges a flow on the line by what crosses, it
  ends not converged when the mass has come to rest on its spring and
  those steps are more than its wires may be left with (after 3569 steps
  on the build this was written on). With `Slip::sliding(0.2)` the 6000
  steps are taken. Asserted: the first is refused as not converged, the
  second is not refused.
- The table, row by row: over six start slips, three compliances and 801
  released slips each, every row is reached more than a thousand times, the
  three claims hold at every point, and 4806 of the points pass
  `check_load_line` on a matrix that is not symmetric and `check_component`
  off any line.
- A mass thrown at 1.3 m/s against 5.88 N stops, and the step it stops in ends
  at zero slip: the arrest is late against `v0 / a` = 0.4422 s by 7.8e-3,
  2.8e-3 and 3.2e-4 s at steps of 10, 5 and 2.5 ms. What is asserted is the
  bound, late by less than a step; the lateness itself is the distance to the
  next multiple of the step and is not a smooth function of it. It rests
  within a step's travel of `v0² / (2 a)`, and its kinetic energy is in the
  body that takes the slider's heat to 1e-12 (by the ledger; the body's own
  temperature agrees to the rounding of adding to 300 K). At rest after a
  slide the slider's slip is zero exactly and the mass's velocity is zero to
  rounding, not to the bit (3.5e-15 m/s): the two keep their own records and
  agree as closely as each step's loop was closed. What is left flips sign
  every step with a zero mean, so the position holds to the bit, and it costs
  one update a step.
- A mass pushed with 0.999 of the static limit does not move, to the bit, over
  500 steps, takes no update after the first, and its port is published with
  no flow; pushed with 1.001 of it, it accelerates at `(F - mu_k N) / m`
  (velocity and distance to 1e-12) and the loss is `mu_k N` times the
  distance.
- Held under 2000 pushes between 0.05 and 0.92 of the static limit: none is
  refused, and the mass is at rest to the rounding of the released flow.
  Since every input of the component is released before its first ask,
  none of them ends with a flow or a velocity that is not exactly zero
  (282 did before, 1.6e-18 m/s at worst), so that the audit is told the
  rounding rests on the planted negative. The planted negative: a net
  force planted on a held mass, with the solve told to take a wire that is
  open by four times its own size, leaves the junction `c dt F²`. (Up to
  ef1631e0d2 a tolerance of 1e-2 planted it. The solve now judges a line
  flow by what crosses, makes an update at 1e-2 and the plant is gone.)
  The audit's tolerance, read from a refusal, is 1.31e-17 m/s at the
  slider's port, 1.5 times `16 eps (|released| + |C effort|)`. It is the
  sum of two terms, asserted by name to 5 %: that one, the rounding of the
  published flow, and `16 eps |released|`, the wire's own, told where a
  step is taken with slopes; the two are `last_allowance` of the step that
  is taken, which a restore then forgets. A hundred times the tolerance
  is refused by the junction's audit by a factor between 80 and 120, a
  quarter is taken. At the solve's ordinary
  tolerance the force that was refused a hundred times over never reaches
  the audit: it is corrected in one update, and the mass is back on its
  held state to the rounding of the released flow. Against a scratch
  kernel that tells the audit a thousand times the rounding, this test
  alone fails.
- Held under a push that changes: 4.9 N rising to 4.9196 N in 0.4 s,
  nothing rising by 2 N/s, and 4 N sin(2 pi 3.14 t). All 400 steps of each
  are taken (on ef1631e0d2: 1, 222 and 35, then refused by the junction's
  audit for ever), with no update in any step (ts2, a370fb3526,
  2026-11-14; 768, 797 and 799 updates, two at most in a step, at
  0e3df343f2). At every step what the audit is told for the wire is
  between the rounding of the published flow and 1.5 times it, read back;
  the flow left at the slider's port and the mass's mean velocity are
  inside it (0.017, 0.020 and 0.020 of it at most; 0.49, 0.02 and 0.02 at
  0e3df343f2).
- Held under a slow push for as long as it lasts: 3 N rising to 3.3 N in
  20 s for 2000 steps, and 6 N + 2 N sin(2 pi 0.05 t) for 16000. On
  74529fda7f the first was refused after 1728 steps and the second after
  14939, each over the audit's allowance by less than a thousandth. The
  flow and the mean velocity stay within 0.6 of what the audit is told:
  the half the solve takes the wire at, and the rounding of forming the
  residual. They reach 0.014 and 0.020 of it (ts2, a370fb3526; 0.50 at
  0e3df343f2).
- Held beside another mass (2 kg, joined to the held one, which is under
  the first of those pushes): at 0.2 m/s on a spring of 100 N/m, at
  0.2 m/s on a damper of 5 N s/m, at rest on the spring. 2000 steps each,
  all taken, with the same per-step assertions: 0.020, 0.024 and 0.017 of
  what the audit is told at most, with 2000, 2001 and 0 updates (ts2,
  a370fb3526; 0.50 in each at 0c0643548b). On 74529fda7f: refused at the
  second step, at the first, and after 68.
- A block riding a cart: the slider between two bodies that both answer.
  0.5 kg on 2 kg, the cart pushed with 5 N against a spring of 100 N/m to
  a wall. 6000 steps, all taken (on b47fe4f2a8: 248, then refused), stuck
  in every one; the block has the cart's velocity to 1e-12 m/s and the
  friction is the block's mass times its change of velocity to 1e-9 N,
  between 1.07e-4 and 0.99999 N: it passes through zero as the two swing.
- A standing property: 160 drawn pairs, a mass on the ground through the
  slider joined by a spring, a damper or both to a second mass that moves
  or rests, pushed by a constant, a sine, a ramp or nothing, 2500 steps
  each at 1 or 0.1 ms. No step is refused. The slider holds throughout in
  144 and slides at some step in 16. (A junction's audit refused 121 of
  the 160 on ef1631e0d2 and 76 on 74529fda7f.)
- Standing properties, the attack families of the review of 0e3df343f2 and
  b47fe4f2a8; their counts and seeds are in the test for a later reviewer
  to widen. 40 drawn carts, the block riding the cart under a push of 2 to
  8 N and a spring of 50 to 400 N/m, 4000 steps of 1 ms each: no step is
  refused, stuck in every one (on b47fe4f2a8 a junction's audit refused
  33 of the 40). And 300 drawn graphs of 1500 steps: two or three masses
  with drawn velocities, the slider to the ground or between two of them,
  springs (a quarter of them differenced) and dampers between neighbours
  and to a wall, a push up to 14 N on one, steps of 3e-6 to 1e-2 s. No
  step is refused; the slider slides at some step in 208 of them; 492653
  updates in their 450000 steps on a Linux host (ts2, a block's differenced
  columns resolved, 2026-11-14; 536349 at 70e5607b6a before it; 697264 at
  3075f215f3, with
  turned recorded stuck; 697206 there before it; 695720 at 3a42d0c16d;
  695799 there before the
  columns of a line were resolved, 2026-11-12; 695778 on a Mac before, the
  two hosts' libm apart; 691915 on b47fe4f2a8, which also took them all:
  this family did not show the defect, the carts did).
- The planted negative for the rounding of the largest wire of its kind:
  three masses in a row on two springs, the first held by the slider with
  half the static limit, the other two at rest. The third mass's force
  wire is told `16 eps` of the block's largest force, the 9.8 N of terms
  of the held mass's balance (`last_allowance`, to 5 %). A net force
  planted on it, with the solve told to take a wire open by four times
  its size, moves the mass and leaves its junction that force times the
  motion; the audit allows the rounding term times the motion (the
  tolerance it reports, to 5 %). A hundred times the term is refused by a
  factor between 80 and 120 (100.0), a quarter of it is taken with no
  update. Only in a block with a line: in the twin without one, the first
  mass kept where it is by a spring to a wall under the same 9.8 N of
  terms, the third mass's wire is told nothing, and that quarter is
  refused by a factor of 1e9.
- Controls on scratch kernels, each failing the tests named and no other
  of this file: the other wires of a block with a line left on the block's
  scale (the cart; the drawn carts); no rounding of the largest of its
  kind (the 160 pairs, not converged; the planted negative); that rounding
  left out of the second test (the same two); that rounding in every
  block, and the audit told a thousand times it (the planted negative,
  each); a line flow taken at all of its allowance (held beside another;
  the slow push; the 160 pairs; the carts and graphs); a restore that
  keeps the last allowance (the held planted negative).
- Break-away beside a differenced neighbour at rest (break-away): the mass
  on a spring of 400 N/m that holds 2 N of a push rising through the 11.8 N
  the slider and the spring hold together. The spring's tangent declared or
  differenced, the mass breaks away in step 400 and the run goes on to 1000
  steps; the two end 2.7e-15 m apart (2.5e-15 m on ts2 at 70e5607b6a,
  before a block's differenced columns were resolved; this entry said
  1e-17 m). (Until 2026-11-12 a stated limit: with
  the spring differenced, step 400 was refused as not converged, the mass's
  force wire flipping by the static limit less the sliding force, 3.92 N.)
- A block's differenced column (`tests/elements.rs`, the capacity near
  rest): two bodies of 0.02 J/K joined by a conductance of 1.98 W/K, one
  4e-9 K above the other at 300 K, one step of 10 ms. Added without their
  tangents, the first difference moves each body's temperature by about
  1e-3 of its last place, and the entry is nothing; an update made with it
  is a fixed-point iteration of factor 0.99. Not resolved (70e5607b6a, ts2,
  2026-11-14), the step was refused: `NotConverged` after 25 iterations,
  the worst input `wall.a` at a residual of 1.6e-9. Resolved, the step
  takes the one iteration of the same bodies with their tangents declared
  (an affine block), no entry left unresolved, and ends at the same
  temperatures (equal on ts2; the test asserts a last place).
- Differenced columns of a line (`tests/load_line.rs`):
  - The 18 masses: a probe on a line beside a 2 kg mass added without its
    tangent, at 0, 1 and 100 m/s under 1 N, 1e-3 N, 1e-6 N, 1e-9 N, 1e-12 N
    and none, in a step of 1 ms: handed `dt / (2 m)` within 1e-6 in all 18
    (8.3e-8 at worst), no entry left unresolved. Before, 7 of the 18 were
    handed nothing and 12 were off by more than 1e-6 (1.3e-2 too much at 1
    m/s under 1e-3 N).
  - The two-port neighbour: one under 1e3 N at `a`, whose force at `a` moves
    `b` (1e-5 m/(N^2 s) times its square) at once and `a` (100 m/s and
    2.5e-4 m/(N s)) only at a step 401 times longer: each entry within 1e-6
    of its slope, where `b` taken at the longer step would be 3.0e-6 off;
    one entry counted, `a` under the force at `b`, which does not answer it.
  - The jump: the same neighbour under 1e-6 N, with no mobility at `a`
    within 1e-4 N of it and 1 m/s faster above, as much slower below: `a` is
    handed nothing, not the 1024 of the difference across the jump that the
    longest step, `2^-10` N, crosses (a jump the other way as far on the
    other side agrees); two entries counted.
  - The kink: the same, with ten times the mobility at `a` beyond 1e-4 N
    more: `a` is handed what the first step had, nothing, not the secant
    across the kink inside the longest step (2.27e-3 where the slope is
    2.5e-4); two entries counted.
  - The cube: the same, with 5e-5 m/(N^3 s) times the cube of the departure
    from the push at `a` (the two sides agree): within 1e-6 (1.6e-7, the
    cube's 1.9e-7 less the rounding of the 100 m/s) where the step aimed at
    the rounding, some 6.0e-3 N without the bound, gives 7.1e-6.
  - The aimed step: under 1e3 N, with 1e-8 m/(N^3 s) times that cube: within
    1e-6 at the aimed step, where the longest step, 0.98 N, gives 3.8e-5.
  - The side not taken: the same neighbour asking to shorten above its push
    of 1e-6 N: its longer steps taken on the side it steps, but not checked
    on the other; `a` is handed what the first step had, nothing; two
    entries counted.
  - The declared neighbour: beside a differenced mass, a mass added with its
    tangent hands the line the bits it hands it beside a declared one
    (0x3f25d867c3ece2a5).
  - Controls on scratch kernels, each failing the tests named and no other
    of these nine: no fix at all, and no longer step (all nine); every entry
    taking the last step (the two-port neighbour, the jump, the kink, the
    cube, the aimed step, the side not taken); no jump guard (the jump); no
    check on the other side, or the check made on the same side (the kink,
    the side not taken); a check that compares at the wrong sign and so
    rejects every entry (the 18 masses, the two-port neighbour, the cube,
    the aimed step, the declared neighbour, break-away); two longer steps
    (the 18 masses: one entry counted at 100 m/s under 1e-12 N); no bound,
    or a bound by the count alone (the cube); every longer step `2^26` times
    the last, not aimed (the aimed step). Not pinned by a test: the cap of
    four longer steps (with three all nine pass).
- A side that jumps: a mass at rest stuck to a side prescribed at 0 and
  then at 0.5 m/s. In the step of the jump it breaks away (0.5 m/s cannot
  be caught within a step by 9.8 N); it is dragged at exactly 5.88 N in
  each of 170 steps, 2.94 mm/s a step; the 171st ends at zero slip with the
  0.4 N that needs; then it rides with no force. The slider's slip equals
  the slip of the two sides to 1e-12 at every step. The side's work is
  `m V²` = 0.5 J, half of it the mass's kinetic energy and half the heat,
  each to 1e-12. The same with the slider's ends exchanged (the mass at
  `b`, the side prescribed at `a`): the same motion and force on the mass,
  the record with the other sign.
- A small jump under a push: a mass held against 6 N on a side that jumps
  by 1 mm/s along the push. Ending the step of the jump at zero slip would
  be a force against the push with a slip that asks for one along it, a
  negative loss; the slider holds on the mean instead (-2 N, the mass ends
  1 mm/s ahead), and catches that slip in the next step (-8 N, 4e-6 J
  lost); then -6 N. Forces to 1e-9 N, velocities to 1e-15 m/s, the heat to
  1e-12. This is the stuck slider's hold on the mean with a start slip,
  reached through the kernel.
- A side that accelerates: the numbers of "Not covered" above.
- A side told its span (2 kg, 9.8 N held, the sides of "Not covered"):
  - It leaves rest at 4 m/s² at `t = 0`; it rests one step first; and that
    under a push of 15 N on the mass from the moment it moves. At steps of
    1, 0.5, 0.25, 0.9 and 0.7 ms the force is the needed one at every step
    (8, 8 and -7 N; to 1.3e-10, 1.5e-10 and 7.5e-9 N, the last being the
    billionth of a step in which the push rises), the pair never slips,
    the slider's record is the slip the two end each step with to 1e-12
    m/s (asserted in every test of a told side that has a mass), the mass
    ends on the side to 1e-12 m/s, and the heat booked is at most 2.7e-17
    J. Taken as held: broken away within 27 steps at the step lengths
    tried.
  - It swings from rest at 2 Hz, at most 6.32 N needed, for 1 s at 1, 0.5
    and 0.25 ms: it never slips. The force is `m` times the side's
    acceleration at the middle of the step to 1.3e-2, 6.6e-3 and 3.3e-3 N
    (first order: the midpoint relation on a speed that is curved within
    the step); the mass's end speed is the side's to 6.6e-6, 1.7e-6 and
    4.1e-7 m/s (second order); the heat booked for a pair that never slips
    is 3.4e-6, 8.4e-7 and 2.1e-7 J (second order). Asserted: the first of
    each within a bound, and each halving or quartering. And the 40 sides
    whose means are moved by a unit in the last place, of which 24 break
    away when taken as held: told their span, none does in 4 s.
  - It speeds up as `(rate / 2) u²` after a step at rest: the force is
    `m rate u` at the middle of the step to 0.02 N a millisecond of step,
    and it breaks away in the step that contains `u*` or the one after
    (at step 103, 205 and 410 at 1, 0.5 and 0.25 ms; `u*` lies in steps
    103, 205 and 409).
  - It jumps, held over each step, its span that one value twice: the
    numbers of the side that jumps above, at either end of the slider.
  - A standing property: 60 drawn sides, piecewise linear with one to five
    changes of slope in 0.4 s, anywhere within their steps but at least
    three steps apart, and slopes within 2 m/s², under drawn pushes, the
    needed friction within 0.7 of the static limit. None breaks away; the
    force is above the largest needed by at most `m` times a quarter of
    the largest change of acceleration (0.250 of it at worst); from three
    steps after the side's last change the mass is on the side to 1e-12
    m/s; the heat is at most 3.5e-6 J.
  - A change of slope inside a step, and changes in consecutive steps: the
    numbers of "Not covered" above, the forces of the two cases asserted
    to 1e-3 N (-0.558, 1.194, -2.752, 6.917, then 4 N; and -4.637, -2.89,
    -2.168, -1.878, 3.507 N, then broken away), and the 2000 drawn sides.
  - The sign of the span, pinned: the same side told its span by its
    `FlowSource` and by boundaries (minus the mean on the port, the
    waveform's values on the span), a side that swings and one that
    accelerates under a push. What the slider receives, the flow and the
    two velocities of the span, is the same to the bit at every one of
    1000 steps. The forces of the side that swings differ in 770 of the
    1000 steps, by 5.2e-12 N at most, and the mass's velocity by 1.4e-15
    m/s at most; those of the side under a push differ in none (asserted
    to 1e-9 N and 1e-12 m/s). Not to the bit throughout: the source is a
    component of the block, the solve has other unknowns, and a held
    slider's force is the mass's velocity change over `c`. The bits of
    the flow at the slider's `b`, an unknown of the block when the source
    tells the span, rest on the LU residue of its Newton update staying
    below half a unit in the last place (at most 3.9e-5 of it
    over the updates traced on ts2 chain fsz, 2026-11-14), so they are not
    guaranteed by construction.
  - The two refusals, planted: a slider told the span of the end its mass
    is on is refused at the first step by name and nothing is published. A
    side at 0.2 m/s and faster whose span is given with the port's sign is
    taken for two steps and refused by name at the third and at every one
    after it; with the mover's sign every step is taken.
  - What the first sign check refused, taken: a side at 1e-9 m/s, and one
    leaving rest at 1e-6 m/s², each beside a second mass at 0.2 m/s on a
    spring, 300 steps, stuck in every one; a side that reverses twice in
    its first step (ends 0.1 and 0.5 mm/s, mean -0.18 mm/s), 20 steps;
    the same with a second such step five steps later, 20 steps. Two such
    steps in a row: the third step is refused by name.
  - Both ends prescribed and told (no line): the slider slides, its loss
    is its sliding force times the mean slip to 1e-12, and its record
    after each step is the slip the two sides end the step with, to 1e-15.
    The slip it is declared to start with is not used: declared sliding
    with the slip the sides start with, sliding the other way, or stuck,
    the forces and the records are the same to the bit.
  - Controls on scratch sliders, each failing the tests named and no other
    of this file: the span ignored and the end taken as held (eight of
    the tests of a told side: all but the refusals, the slow and reversing
    sides, and the jump, which is held); the span taken in the port's sign
    (eight: all but the refusals, the slow and reversing sides, and the
    test of the two ways, which agree with each other in any sign); a
    span on an answering end not refused (the refusals); the sign decided
    by the step's own inputs (the refusals; the slow and reversing sides;
    the changes in consecutive steps); a span of the port's sign never
    refused (the refusals; the slow and reversing sides); the count not
    started anew by a step in order (the slow and reversing sides); the
    first step started from the mean and not from the span's start (the
    side that accelerates; the side that swings; the changes in
    consecutive steps; the 60 drawn sides).
- A side that speeds up: its speed `(rate / 2) u²` given as step means, a
  staircase. Stuck, the slider gives `m (V_n - V_(n-1)) / dt`, which is
  `m rate u` at the start of the step, to 1e-9 N, and the mass has the
  side's velocity to 1e-12, until that force exceeds 9.8 N: it breaks away
  in the first step that starts after `u* = mu_s N / (m rate)`, step 104,
  having been carried with more than 1.6 times the sliding force.
- Stuck to a side that moved, the three rows over 9612 points: carried
  2410, held 2304, broken away 4898; the loss is not negative and the force
  is inside `mu_s N` at every point; a carried force is the one that ends
  the step at zero slip, to the bit, and reaches above 1.6 times the
  sliding force; a held one is `r / c` inside the limit, to the bit. At a
  threshold: `on_line` returns exactly the limit where the quotient would
  be 3.0000000000000004; with no compliance the slip is the given one and
  the force the sliding one. With unequal coefficients (3 N static, 1.5 N
  sliding): held with exactly 3 N at `c mu_s N`, the sliding 1.5 N a hair
  beyond; sliding, turned with exactly -1.5 N at `-c mu_k N` and reversed a
  hair beyond; a start slip that 2 N end is carried when stuck and slid
  over when sliding.
- A mass on a spring with friction (1 kg, 100 N/m, 3 N sliding, 5 N
  static, thrown at sqrt(22) m/s): the turning points 0.44, -0.38, 0.32,
  -0.26, 0.20, -0.14, 0.08 fall by `2 mu_k N / k` a half cycle, and it
  stops at -0.02 inside the dead band. Off the closed form by 4.3e-6 m at a
  step of 1 ms and 1.1e-6 m at 0.5 ms at most; asserted: within
  `8 (mu_k N / m) dt²`. What it was thrown with is in the spring and in the
  heat to 1e-12. Every one of its steps is a point of `Slider::on_line`:
  the force to 1e-9 N and the slip to 1e-12 m/s for the released slip
  rebuilt from the two on the line the mass and the spring present
  (`c = g / (1 + g k dt / 2)`, `g = dt / (2 m)`). That checks the rows and
  their forces, the reversals among them; it does not check `c`, since a
  released slip rebuilt from a held step satisfies the law for any
  compliance (the probe pins `c`). The rows reached are, at 1 ms, slid 2342,
  held 1648, released 4, reversed 2, turned 1, arrested 2, broken away 1,
  and at 0.5 ms 4693, 3298, 4, 3, 1, 1, 0: the three reversal rows are
  asserted to be among them. (The same rows and turning points to the
  digit whether turned ends stuck or sliding.) Each released step is
  asserted to end sliding, as a step with no friction in it does; recorded
  stuck, the step after each one would be judged from stuck and break away
  (at 1 ms: broken away 5, slid 2338).
- fs-cli's friction belt, built as its command builds it (2 kg, 10 N held,
  4 N sliding, a belt tabulated to its last speed at 0.5 s, steps of 5 ms):
  held through step 50, broken away in 51, sliding at 4 N. At a last speed
  of 3 m/s the belt catches the block at 1.25 s, the end of step 250, where
  the force that ends the step at no slip is the sliding force in exact
  arithmetic. On ts2 at 3075f215f3 step 250 slides (the slider's record
  4.26e-14 m/s, the block's slip 2.09e-14), 251 is released, and from 252
  each step holds the mean slip at zero, at a force of 1.67e-11 N growing
  to 5.01e-11 N by step 300, the two records on opposite sides of zero:
  turned at 252 and held after it, stuck in every one of those 49 steps.
  While turned ended sliding, all 49 were turned and none stuck. At
  3.005 m/s the catch is inside step 251, arrested there at exactly zero
  by 2.0 N; at 3.000000002 m/s step 250 slides and 251 is arrested at
  exactly zero by 8.0e-7 N. Asserted: stuck from step 252 (catch + 2) to
  1.5 s in the first; the arrest and stuck after it in the other two.
- The matrix a component is handed, read by a probe that returns no effort:
  `dt / (2 m)` per mass with nothing between two masses; the closed form
  with a damper between them (symmetric, off-diagonal entries not small);
  with a gyrator between them (`[[g1, -g1 g2 G], [g1 g2 G, g2]] /
  (1 + g1 g2 G²)`, not symmetric); a zero row and column for a prescribed
  side; `g1 / (1 + g1 k dt / 2)` for a mass on a spring to a wall (with
  3e6 N/m, which takes more than a fifth off). Each entry to 1e-12.
- The slider's forces through such a matrix: the two masses (2 kg thrown at
  1 m/s, 0.5 kg at rest) rub on each other beside the damper and beside
  the gyrator, so every flow of the line is the released flow less what
  both forces take away, one through an off-diagonal entry. Beside the
  damper, with 30 N on the first mass: the momentum is the 2 kg m/s and
  the push's impulse to 1e-12 at every step; the second mass needs 6 N to
  go along, more than the sliding 5.88 N, so the slider slides at every
  one of the 80 steps, at the difference at which the damper gives the
  missing 0.12 N (0.171 mm/s, to 1e-12); the push's work is what the two
  gain and the slider and the damper lose, to 1e-12. Beside the gyrator
  (not symmetric): the two lose exactly the slider's heat at every one of
  6000 steps of 0.1 ms (to 1e-12 J of 1 J); they rub in 2343 steps and are
  stuck in 3625; they end moving as one at 5.98 mm/s, held with the
  5.38 N the gyrator pushes them apart with (`G v`, to 1e-9 N). (On ts2 at
  3075f215f3. While turned ended sliding, the step after a reversal held
  at no mean slip could take only `mu_k N`: 2346, 3616, 5.96 mm/s and
  5.36 N.) Against a
  scratch kernel that publishes the line's flows through the transposed
  matrix this test alone fails, at the gyrator's first step
  (1.55e-8 J against 1.74e-11 J).
- A refused step leaves nothing of its line behind: with a spring that the
  kernel differences beside the slider, a step asked with no iteration
  allowed is refused after the line was found, and a step asked with one
  iteration allowed is refused after an assembly and an update (the slider
  differenced, proposed at a perturbed flow), three times each in 60 steps;
  every step that follows equals, in iterations and bits, the step of a
  system never refused.
- The row without friction, planted: a mass sliding at 1 mm/s pushed back with
  `1.5 m v / dt` reverses in one step with no force from the slider and no
  heat, and enters the next step sliding the other way; friction then brings
  it to rest.
- A mass thrown onto a free slab (a line over both ends): the momentum holds
  to 1e-12 at every step, the two join at the velocity of their centre of
  mass late by less than a step, and the kinetic energy of the relative
  motion is the heat to 1e-12.
- A block on a belt held by a spring (the belt is the slider's other end,
  prescribed): against the closed-form stick-slip period
  `2 d / V + 2 (pi - atan(d w / V)) / w` = 0.80689 s the measured period is
  off by 8.9e-4, 1.1e-4, 1.1e-4, 1.1e-4 and 1.4e-5 s at a step of 2 ms and
  four halvings. Asserted: off by less than a step. The discrete motion
  repeats after a whole number of steps, so the error moves in jumps.
- A slider between prescribed velocities whose heat warms a body through a
  junction: no line, a loop through the heat port, the body takes force times
  slip to 1e-12.
- A run with a line and its replay through a step of another length and a
  restore take the same iterations and reach the same bits over 80 steps:
  sliding, at the arrest and at rest.
- Refused: a line on a heat port, over no port, out of order, over a port
  twice or over a port the component does not have (at `add_on_line`); a neighbour of negative
  inertia, which gives a line of negative compliance, and one of no
  inertia, whose answer is not finite (the first step, the slider named
  with the reason, nothing published). A slider with both ends prescribed
  and no loop is stepped by its mandatory tier with no iteration.
- Multiple owners, `tests/multiple_load_lines.rs`: two brakes on one moving
  body against constant-force motion and heat, two spring-coupled bodies with
  one brake sticking and one sliding against an independent midpoint
  oscillator update, iteration-budget refusal and cancellation without
  publication, and bitwise state/wire/account replay across different step
  lengths. Two loaded sticking brakes on the same body exercise the named
  singular-block refusal; two unforced brakes at rest exercise acceptance of
  an already consistent fixed point in `tests/load_line.rs`.
- Requests to shorten inside a loop, with a test component that cuts the
  step at a stop (`CutsAtAStop`: one constant force to the ground, the stop
  located from the line as for a body under a held force; its heat warms a
  body whose temperature it is handed back, so every step takes a second
  update). A mass of 2 kg thrown at 1 m/s onto 3 N stops after 2/3 s and
  1/3 m.
  - With no starting flow given to the port, the first step of 0.05 s is
    taken whole, with no shortening, and leaves the mass at 0.925 m/s (to
    1e-15). Before: `StepUnderflow { dt_s: 1.26e-21 }`.
  - Through the stop, with the starting flow given and without: 13 steps
    of 0.05 s, then the step from 0.65 s cut once, ending at
    0.6666666666666669 s at rest (slip exactly zero, the mass at
    2.2e-16 m/s and 0.3333333333333335 m), then six whole steps at rest;
    one shortening in the run. Before, with the starting flow given: that
    step cut twice, to 0.016666666169961288 s, and not at rest.
  - No load line: a damper that the kernel differences, between a pushed
    mass and one it drags, asks to shorten where the flow at one end is
    within a micrometre per second of zero and not zero, which no iterate
    is and the kernel's perturbed proposal about a flow of zero is. Asking
    above zero only, the step of 1 s is taken whole in one update, and the
    two masses end where the midpoint rule puts them (to 1e-15): the
    column from the other side is the damper's slope with its sign.
    Before: halved three times, to 0.125 s. Asking on both sides, the step
    is halved three times, until the damper can be differenced.
  Planted in a scratch kernel, each failing these tests and no other of
  the file: a perturbed request shortening the step at once (the second
  and third test); the line not released before its component is asked
  again (the first and second); a request on both sides ignored (the
  third); the other side's column divided by the forward step (the third).
- Not exercised by a test here: a neighbour that is not linear (the line is
  then its slope at the step's start). `crates/fs-couple/tests/legacy/tribothermal/
  dynamics/collinear/graph.rs` has a friction component beside springs
  with a cubic term.

The first product on the tier, `crates/fs-couple/tests/legacy/tribothermal/
dynamics/collinear/graph.rs` (run with `cargo test -p fs-couple --test
tribothermal collinear::graph`): the hand-written
`fs_couple::tribothermal::dynamics::collinear`, which advances two
translating bodies under held forces through a stop (each sliding segment
on the midpoint relation, a segment that would reverse cut at zero relative
velocity, from rest the friction law's own static map), against a graph
with no pairing code. The bodies are `Mass`es under `Force`s, one pair
along the tangent axis of the case and one across it; between the pair of
the axis is fs-tribo's `tangential::component::FrictionInterface` on the
load line, with the Coulomb coefficients the tests' interface card holds;
the finite bodies are fs-material's `UniformEnthalpyBody` on the tests' own
caloric charts; a heat route divides the friction heat by the tests'
partition. The module is the reference and stays.

- Three of the module's own cases (a stop that lasts, a stop left the
  other way under the drive, a pair that never slides), each at steps of
  0.05, 0.07, 0.013 and 0.2 s: the instant of the stop, the velocities and
  displacements of both bodies on both axes, the kinetic energy, the
  external work, the heat each body took and its temperature are within
  1.3e-11 of the module's report, each in its own unit (asserted: 1e-10,
  the module's own event tolerance); the graph's first law over the run is
  off by at most 2.2e-12 J (the module's own residual: 2.0e-11 J). A stop
  inside a step cuts that step once; a stop at a step's end cuts nothing;
  nothing else is cut.
- Which element gives which arrest. `elements::Slider` with the same force
  in the same graph takes the stop inside its step, and on the second case
  it ends body `a` at -1.42 m/s at 0.07 s and at -1.396 m/s at 0.013 s where
  the module has -1.4 (heat 0.4366 J and 0.4131 J against 0.42 J): first
  order in the step, with its own first law exact. Those are the numbers
  an emulation of the two arrests gave before either graph was run. The
  interface, which asks for the step that ends at the stop, is within
  1e-10 at both.
- Closed forms through the interface (not the module's): a mass on a
  spring loses two friction offsets a half cycle and its eight turning
  points are the closed form's to 5.3e-16 m at steps of 1 ms and 0.7 ms,
  with all it was thrown with, less what the spring keeps, as heat; a
  block on a belt sticks and slips with both events located, at belt
  speeds of 0.0937 and 0.137 m/s and steps of 2 ms and three halvings:
  every break-away within 5e-11 s of the instant the spring reaches the
  limit, with one cut, every period 3.35 and 3.62 dt^2 off the closed
  form's (asserted: under 4 dt^2), the same from cycle to cycle to
  5.3e-11 s, the work of the two ends being the heat; a push that rises
  at 10 N/s through the 5 N held lets go at 0.5 s to 7e-16 s with one cut
  at steps of 0.03, 0.007 and 0.0011 s, the velocity after it the closed
  form's to 1.8e-15 m/s and the place a quarter of dt^2 behind; a push of
  8 N sin(pi t) lets go within 1.8e-6 s of the instant after three to five
  cuts (asserted: within the impulse tolerance's 3.2e-6 s); a side that
  jumps from 0.1 to 0.3 m/s drags the block under the kinetic force and
  catches it 0.2 / 3 s later, with one cut; the time mean of the interface
  temperature between two lumped bodies behind two conductances, and the
  bodies' end temperatures, are the closed form's to 1.1e-7 K and 3.6e-7 K
  at 0.05 s, second order in the step.
- Which element locates a break-away. `elements::Slider` with the same
  forces under the same pushes decides each step by its mean: it lets go
  0.33, 0.43 and 0.45 of a step off the instant on the rising push (its
  velocity then off by 2.1e-2, 6.0e-3 and 1.0e-3 m/s) and 0.16 to 0.35 of
  a step off on the sine, with no cut. The interface with its location
  planted out does the same.
- Beside a spring the stop is not located in one cut (fs-tribo's contract,
  "Where the stop is found"): two or three cuts in the step that holds it
  beside a linear spring; 2 to 5 beside a spring with a cubic term; with
  two bodies and that spring at 0.07 s the cut step ends short of the stop
  and three steps of 3.6e-7, 9.1e-9 and 2.3e-10 s follow. At rest in every
  case.
- Outside the claim, asserted: a prescribed side that accelerates. Under
  a side that speeds up at 1 m/s2, which the block could ride with 1 N of
  the 5 N held, every step asked for is followed by a train of steps, each
  a third of the one before (the first eight to 1e-9), until what is left
  is inside the impulse tolerance: 235 steps to reach 0.2 s at a step of
  0.01 s and 2129 at 0.001 s, with 1.045 and 1.117 times m (a dt)^2 / 2 of
  heat for each whole step where a carried block makes none, and the block
  a step behind the side. `elements::Slider` can be told the span of such
  a side.
- Two things in 1441d587a2's record of this graph were wrong, and are
  corrected above. The belt ran at 0.1 m/s, where the ride of 0.4 s is a
  whole number of steps at each step length: every break-away fell on a
  step's end, and the period came out at 1.4e-5 to 2.1e-7 s, second order,
  an accident of the numbers that the record gave as the result. (At the
  speeds above, before the break-away was located, it was taken from 0.48
  of a step before the instant to 0.50 after and a period was up to 0.54
  of a step off, different from cycle to cycle: 9f4fd342b2.) And the
  interface returned, at a stop, the force that ends the step at exactly
  no slip, in a row narrower than the kernel's perturbation (next
  paragraph): the belt at 0.137 m/s in steps of 0.25 ms was refused at its
  third catch after 9009 steps, and the accelerating side above ended in a
  step that did not converge. Both are tests now.

A rule for a component on a load line follows from the second, and is the
component's to keep, not the kernel's to check. The kernel differences
the component at a flow moved by `sqrt(eps)` of its scale. A row of the
component's relation that is narrower than that is not seen: the column
is the neighbouring rows'. So in such a row the component must return an
effort with the slope of its neighbours. `elements::Slider` meets the rule
(its arrest row is as wide as its kinetic force). fs-tribo's interface as
first landed did not: at a stop it returned a force of slope `1 / c` in a
band of 5e-11 m/s between rows of slope zero, and a step that needed one
more update at such an iterate ran Newton with a slope of zero against a
row of slope 8000 N per m/s and never closed its force wire below
1e-11 N. A tangent the component declares for the row it takes would
remove the rule; the load-line tier has none yet.

Not in the graph: coefficients that follow a material card at the trial
temperatures and speed (the module's third test), the thermal step as a
backward-Euler solve of a cell network with other sources, two tangent axes
solved together. Nothing in the module is deleted.

## No-claim boundaries

- No accuracy claim for any component; the kernel checks accounts, not physics.
- The end of a flow that is held at exactly zero. A wire is solved relative
  to its size and to the rounding of the terms its driver sums. An end held
  at nothing has no size, and its rounding is not in those terms: a modal
  neighbour forms it as a sum over its modes, which cancel at the port and
  do not vanish in themselves. The kernel is not told that sum, and with no
  absolute tolerance the block is refused within a few steps (asserted in
  `tests/load_line.rs`; on ts2 at 3a42d0c16d, step 33 of the test's graph, a
  residual of 1.4e-17 m/s). The solve's absolute tolerance has to cover it.
  An end carried at a speed that is not zero (a bow) has a size and is
  taken.
- A component that keeps the default promise of `end_of_step` and does not
  have it. Nothing can check the promise from outside: the kernel sees the
  mean. Beside the slider such a neighbour is held at no mean slip and
  never brought to rest (asserted in `tests/load_line.rs`).
- The rows of a line found by differencing are known to what a forward
  difference leaves. Beside two exact oscillators the mean and the end row
  were within 1.4e-6 and 1.3e-6 of the closed form at a step of 1 ms and
  within 2.1e-7 and 4.5e-7 at 4 ms (measured there on ts2 at 3a42d0c16d,
  each column resolved against the rounding of its outputs; asserted within
  1e-4).
- The elements are ideal and their parameters constant. A `HeatCapacity`
  stores its absolute temperature, so a small heat increment is rounded to an
  ulp of that temperature each step; its ledger account, a compensated sum of
  receipts, does not lose it.
- Streams: the kernel audits enthalpy, connection mass cancellation and
  each measured owner's total mass balance. An owner returning `None` has no
  inventory verdict, and the system mass residual stays unavailable without
  complete measured coverage. Carried species or constituent inventories are
  audited only when the owner exposes their inventories, as described under
  "Carried stream inventory conservation" below. Entropy carried by streams
  and made by mixing is not accounted. Entropy made by throttling is
  reported only by a `Restriction::for_gas`; `new` and `for_stream` report
  none. Nothing checks that joined volumes hold the same gas. A
  `GasVolume` shortens a step that passes more than about twice its mass
  through it, so a small volume between restrictions sets the step; a
  prescribed draw of a quantity a volume holds none of ends in step
  underflow. Stream inputs inside a loop start from zero on the first step.
- Field metrics are fixed, dense and symmetric positive definite; indefinite,
  nonsymmetric and state-varying pairings are not admitted. No conditioning
  bound or continuum trace-map accuracy is inferred from metric admission.
  No gyrator as a connection (the `Gyrator` element is a component);
  the non-smooth tier solves multiple load-line owners together, each with its
  self-response frozen at the step's start, and locates no event itself (a component
  on the line may ask for the step that ends at one; `Slider` does not, and
  a step in which its friction switches is first order). No general
  multi-contact complementarity algorithm or force-sharing selection is
  supplied: redundant active constraints may leave a singular block, which
  is refused with its component names; convergence and uniqueness are not
  guaranteed. No multirate stepping; no
  monolithic solve through `StageResidual`; no embedded error estimator or
  PI controller (step doubling and proportional control are available).
  All `elements` and the optional nonlinear electrical storage adapters have
  native state codecs; additional physics-owner components need their own
  codecs before the whole system can be snapshotted.
  No fork lineage or SoA buffers (bead `.5.4`); no adjoint. These are beads
  `.5.2` to `.5.10` of the epic.
- An entry of a differenced column, a line's or a block's, is resolved
  against the rounding of its output, `eps` of its size. An output that is a small sum
  of large terms (a modal structure's velocity at a port, a sum over its
  modes) rounds by the terms: its entries can be further off than 1e-6, and
  its two sides can differ by more than `2^-10` with no kink. On a bowed
  string of 64 modes at 48 kHz, a fixture outside this repository, run to
  its refusal at step 784 (it is refused there with the columns resolved or
  not, the step shortened too often): 998 columns differenced again in its
  784 steps, each at one longer step and once on the other side; 1053
  entries taken, their two sides within 9.7e-4 of the entry (657 of them
  within 7.6e-6); 60 left unresolved, 47 of them because the sides differed
  by 1.005e-3 of the entry to more than half of it. The line's mobility as
  the friction law had it: off the string's own by 1e-5 to 1e-4 in 636 of
  its calls at the full step at worst, where before 27 calls had it 13 per
  cent off (2.976e-1 for 3.428e-1; 3.7 per cent on a Mac). A component that
  declares its tangent is handed to the line as declared, and none of this
  applies to it.
- A step is slow against hand-written code: 21 to 160 times in the cases
  measured. Measured by `tests/step_cost.rs` (ignored by default; `cargo
  test --release -p fs-system --test step_cost -- --ignored --nocapture`)
  on a line of masses joined by springs, with or without dampers, against a
  hand-written loop for the same midpoint equations (one tridiagonal solve
  a step), and on modes that meet at one contact against a hand-written
  loop in the one scalar (two passes over the modes a step); the graphs
  reproduce both to 1e-12. On an Apple M5 Ultra under other load,
  2026-11-11, nanoseconds a step (the first two columns at load average 19
  to 25, the rest at 23 to 29):

  | masses | dampers | first measurement | dense, kept factors | now | hand | now / hand |
  |---|---|---|---|---|---|---|
  | 1 | yes | 1,145 | 842 | 830 | 8.3 | 100 |
  | 4 | yes | 9,401 | 4,581 | 3,129 | 38 | 82 |
  | 16 | yes | 112,718 | 55,880 | 10,359 | 152 | 68 |
  | 48 | yes | 1,109,593 | 633,621 | 37,630 | 508 | 74 |
  | 1 | no | | 380 | 322 | 9.0 | 36 |
  | 4 | no | | 1,193 | 928 | 37 | 25 |
  | 16 | no | | 5,317 | 3,247 | 146 | 22 |
  | 48 | no | | 23,695 | 10,489 | 509 | 21 |

  | modes at one contact | now | hand | now / hand |
  |---|---|---|---|
  | 1 | 258 | 9.9 | 26 |
  | 16 | 1,374 | 21.5 | 64 |
  | 200 | 27,267 | 172 | 159 |

  The first measurement re-derived the Jacobian by asking every driver about
  every unknown, factored it every iteration, and allocated 12 to 19 times a
  step. The second column kept workspaces between steps (an accepted step
  allocates nothing, asserted), assembled the Jacobian from a list made at
  sealing of which residual rows read which component's outputs, and kept
  the factors of a block that is affine throughout (every component
  `Class::ExactLinear`, every driver linear) for as long as the step length
  stays the same; its results were those of the first measurement's code to
  the bit. "Now" substitutes the inputs that are copies and eliminates only
  the ones that close the loops (the invariant above): 49 of 288 inputs for
  48 damped masses, one of 201 for 200 modes. Results moved at the rounding
  of the elimination and every suite below passes unchanged: of the reports
  of ten example scenes stored before the change, three are the same bytes
  (no loop, or one component) and seven moved, the fixed-step ones by at
  most 8.9e-15 in any number, the adaptive one by 2.8e-13 s in the times its
  controller chose and 9.3e-12 in its last normalized error; no step or
  iteration count changed. The cost is now close to linear in the number of
  components in these cases: 81 to 136 ns a component a step with one
  iteration, 166 to 258 ns with two.

  What is left, and what bead `.5.7` still has to remove: the kept inputs
  are eliminated densely (49 by 49 twice a step for 48 damped masses; the
  elimination skips zero multipliers, so for a line, whose kept matrix is
  tridiagonal plus the heat sum's row and column, that part grows with the
  square, and with the cube where the kept inputs all depend on each
  other; not measured apart); a damper makes a block not
  affine, because its heat is quadratic in its rate, so a damped line
  assembles and factors twice a step although its mechanics are linear;
  every component is proposed once more than it is iterated, audited and
  published through dynamic calls whatever its class; the plan is not
  explainable beyond `plan_shape` and `plan_loops`. A graph whose modes are
  each a component pays the per-component cost 200 times where the
  hand-written loop does two passes over an array; a modal bank that is one
  component would not.
- Entropy is accounted as generation and as flow through heat ports. Stored
  entropy of components is not tracked, so the second-law check is narrow: it
  refuses a step in which the entropy components say they generated, plus the
  entropy that appears or vanishes across connections, is negative. It does
  not check the entropy balance of each component.
  Of the stream elements, only a `Restriction::for_gas` reports
  generation: that of its throttling. A `GasVolume` does not report what
  mixing a stream into it makes, and that is never negative. So the entropy
  a stream system's ledger shows generated is a lower bound of what it made.
- The audit trusts the tolerance a component claims for its own energy
  balance. A component that claims a wide tolerance passes the audit with a
  leak; the leak is visible in the ledger's `balance_defect_j` and
  `energy_residual_j`, not refused.
- `common_flow` on heat ports is admitted for one member only (a bond). Heat
  shared between several bodies goes through `common_effort` or `heat_route`.
- The kernel does not check the entropy balance of a component, and a
  `ConductanceMatrix` does not close its own. It reports the entropy
  generated where the heat crosses, at the interpolated temperatures of its
  interface points. The kernel books the entropy through a field heat port
  by the nodal rule, `dt sum Q_i / T_i`, at the temperatures of the trace's
  points. Minus the nodal sum over the element's two ports is therefore not
  what the element generates: the two differ by the difference between a
  trace's temperature at an interface point and at the points it is
  interpolated from, with either sign, and the nodal sum can be negative
  where the generation is not. `tests/conductance_matrix.rs` pins two such
  fields, both taken: opposed jumps between points at 1200 K and 300 K
  (generated 3.8e-5 J/K, nodal sum -2.8e-6 J/K), and one point at 314.5 K
  fed in equal parts from 300 K and 330 K (2.5e-6 against -1.1e-6 J/K).
  The ledger shows the gap between a contact's generation and its ports'
  nodal sum; no order is claimed for it.
- A field heat port names no trace. Nothing checks that two joined field
  heat ports list the same points in the same order, or that a
  `ConductanceMatrix` was built for the two traces it is bonded to.
- Momentum and charge outside carried stream coordinates are not in the ledger
  yet (bead `.5.6`). Coordinate accounting alone supplies no chemical identity,
  reaction balance or mechanical momentum theorem.

## Carried stream inventory conservation

`Component::stored_carried_quantity(state, carried)` opts an owner into an
independent inventory audit for one zero-based coordinate after mass and
enthalpy. Both endpoints must be available or both unavailable. Inventories
are finite and may be signed; their declared units and coordinate meaning
belong to the owners. Every stream of an owner uses the same coordinate order.
An available inventory declares conservation under transport: internal
reactions or sources are not silently booked as incoming material.

The kernel audits each coordinate separately with scaled, compensated sums
and a roundoff allowance based on its own endpoint inventories and transfers.
Total mass or energy closure cannot hide an unreported species conversion.
Each connection also closes its carried transfers within the interface solve's
coordinate allowance. Failures identify the owner or connection and ordinal,
and refuse the entire step before any state, clock or account publishes.

`NodeAccount::carried` keeps optional measured accounts; boundary `carried_in`
and connection `carried_defect` retain the transported amounts independently.
`Ledger::carried_residual(i)` is available only after a step and complete
measurement coverage of every owner carrying coordinate `i` over the window.
Later measurements cannot fill an earlier unaudited step. Checkpoint restore
rewinds these accounts and coverage counts with the other accounts.
All added per-port, per-owner and per-connection workspaces are sized at seal.
`GasVolume` measures its `carried_kg`; `Restriction` declares zero storage.
An attempted draw past a constituent's inventory shortens the step using that
constituent's own rounding scale, so bulk gas mass cannot erase a trace species.

G0/G3 tests in `tests/carried_accounting.rs` check distinct gas mixtures through
flow reversal, boundary delivery, exact checkpoint replay, an unreported
conversion that preserves mass and energy, atomic refusal and missing coverage.
They also drain a trace constituent at a mass fraction of `1e-30`: its inventory
limits the accepted interval even when bulk mass and energy barely change.
No species-name matching, equation of state compatibility, reaction chemistry
or entropy of mixing is inferred from this ordinal-coordinate accounting.

## Original port-Hamiltonian Component

The optional `phs` feature exposes `phs::PortHamiltonianComponent` (bead
`frankensim-rc-1109-root-3gc12.5.12`). It owns an existing `fs-phs`
`PortHamiltonian`, retaining its complete operators, source storage and state
order. Each native input/output pair has one caller-declared scalar power port;
`EffortIn` and `FlowIn` specify which conjugate the native input represents.
The adapter validates count, scalar ordinary-power shape and unique identity;
`heat` is reserved. It infers no coordinate, orientation or unit conversion.

Every proposal delegates its actual negotiated duration to
`StepWorkspace::step_into_controlled`: the original variable-step
discrete-gradient equation, native Newton admission and original source laws.
The final mandatory heat port receives positive temperature and emits the
native resistance work as `-D/dt`. Physical heat is never reconstructed from
the energy deficit. The receipt retains that heat and its entropy `D/T`, the
native stored change, and the signed native work-minus-heat-minus-storage
defect as numerical loss. The latter must fit the accepted discrete effort's
L1 norm times the equation residual, plus an arithmetic allowance from endpoint
energies and the original absolute port-work contractions. Negative heat or
unrepresentable accounts refuse. The native effective residual budget is kept;
an explicit `set_solve_budget` retains below-rounding-floor refusals.

Private `RefCell` scratch contains no accepted physical state or clock. The
native solver refreshes its inputs, trial state, structure and iterations on
every proposal; a reentrant call refuses instead of panicking. Adapter outputs
publish only after native and physical admission and the final cancellation
check. Warm calls reuse allocated buffers when the candidate has the original
shape, conditional on the source storage's allocation behavior. Cancellation
polls occur inside native Newton/Jacobian work; individual storage callbacks
and the native dense solve are not preemptible.

G0/G1/G4 tests in `tests/phs.rs` cover multiport native equivalence and
conformance across different states, inputs and durations; nonlinear source
storage; a driven storage depositing its actual loss in a thermal body against
independent midpoint and continuous exponential references; adaptive
coarse/half trials and rollback/retry; cancellation inside a native storage
callback; typed native refusals; and refusal of energy created by an invalid
raw structure. This is not an exact-exponential adapter or a `PortTangent`:
the complete heat-bearing map is `Class::Smooth`. The caller must supply
differentiable source storage on its admitted support and a smooth port map;
switching or contact laws that require their original event treatment are
outside this adapter. Temperature sets the entropy account only; it does not
rebind the owner's fixed resistance or infer a material law. Numerical
conformance does not establish physical validity of the supplied model. The
exact-quadratic companion is described below; other adapter work in bead
`.5.12` stays open.

### Source-owned endpoint evaluation allowance

The PHS adapter reads the original storage's optional
`energy_evaluation_error_bound_j` at both endpoints. A declared bound must be
finite and nonnegative; their receipt allowance must remain finite. `None`
retains the existing ordinary arithmetic allowance. An owner whose Hamiltonian
evaluation cancels large internal terms or an arbitrary energy datum must
declare its complete absolute evaluation error from the source arithmetic.
The adapter adds both endpoint bounds to `Receipt::energy_tolerance_j` only;
its independent native work/heat/storage check retains the original arithmetic
and equation-work allowance. Thus endpoint measurement roundoff cannot grant
permission for unphysical work or be converted into physical heat.

The regression sources in `tests/phs.rs` retain the native 100000000.5 J stored
change when a shifted quadratic's rounded endpoint subtraction gives
100000000 J, exercise the source bound after cold sharing, refuse negative or
nonfinite bounds at either endpoint and an overflowing combined allowance
without changing caller outputs, and retain the invalid-structure refusal even
when the source declares a much larger endpoint bound.

## Exact port-Hamiltonian Component and complete tangent

The same optional `phs` feature exposes
`phs_exact::ExactPortHamiltonianComponent`. It accepts the original owner,
the same ordered scalar power-port declarations, and a native
`LinearSplitBudget`. Native zero-clock preparation admits the declared exact
quadratic storage, structure and dimensions during construction. Each positive
negotiated duration uses `PreparedExactZoh::new_owned_controlled` over
`share_storage`; the original operators, source law and state order are retained.
The adapter retains at most two duration caches, keyed by exact IEEE bits,
evicting the least recently used. Cold misses prepare before eviction, so failed
preparation preserves the existing cache. Caches and reusable scratch carry no accepted state or clock;
adaptive coarse/half trials, refusal, restore and retry require no cache rollback.

Every proposal uses the original controlled exact step and its independent
resistance heat integral and physical balance gate. The mandatory `heat` outlet
returns `-D/dt`, the receipt records `D/T` at the supplied positive temperature,
and the signed native supplied-minus-heat-minus-storage defect remains numerical
loss. The arithmetic allowance includes the original retained two-word port-work
magnitudes before opposing ports cancel, read immediately after that same native
step. The source's checked endpoint evaluation bounds are added only for the
Component's energy remeasurement, as above; neither allowance weakens native
physical admission or turns numerical error into heat.

`PortTangent` uses the cached original mean-input response rows and the native
analytic gradient of the independently integrated resistance heat. Its JVP
includes the full `-dD/du / dt` heat row; its VJP includes that row's output
weight. The temperature column is zero, since temperature sets entropy only.
The complete heat-bearing map stays `Class::Smooth`. Tangent calls first run
the same complete native trial and source-support admission as proposals;
derivatives cannot succeed outside the original physical map's support.
There are no finite-difference derivatives or substitute integration formulas.

Snapshots. The owner has the state codec of fs-ports' `Component` when its
system can name itself. The state is the coordinate count and each
coordinate's bits. It is written, and read back, only as a state the owner
holds: the system's count, finite coordinates and a finite energy. The
identity is the system's own (`PortHamiltonian::identity_words` of fs-phs:
its dimensions, `J`, `R`, `G`, its storage's words and its resistance
witness), then the two limits of the split budget and whether ends are told.
The ports are in a system's own description, and the scratch (a handle on
the same storage, at most two prepared durations, buffers) holds nothing
that is not derived. A system whose storage cannot name itself has no
identity, and the owner then has no codec: `encode_state` writes nothing,
`decode_state` refuses and `identity_words` appends nothing. Today that is
every storage but fs-phs's `QuadraticStorage`. An exact quadratic part is
not enough: the tests' shifted quadratic has one, and its energy and its
evaluation error bound, both of which this owner reads, are its own. The
owner shares its system's storage when it is built, so a plain quadratic
storage is named through fs-phs's shared handle. Three tests in
`tests/phs_exact.rs` (claims written before the first run; on that run all
three failed, because the handle did not yet forward the storage's words):
- the state is read back; no words, a word short, a word over, another count
  and a coordinate that is not finite are refused and not written; an owner
  on the shifted quadratic writes nothing, reads nothing and has no identity
  words;
- the identity is equal for a rebuild and for an owner whose scratch two
  proposals at other durations have used, and a used owner and a fresh one
  read the same words and step to the same bits; the identity differs for
  another entry of the system, each of the two limits of the budget, and
  ends that are told;
- through `System::snapshot`: a system whose owner is on a plain quadratic
  storage is written to a snapshot, and the same system with another `Q` is
  another problem with other words; a system whose owner is on the shifted
  quadratic is refused with `NoCodec { node: "owner" }`.
Not tested here: a resume of such a snapshot (fs-cli's scene tests save and
resume through the binary).

The end of a step. The step is exact for a held input, not the midpoint rule,
so the default `Component::end_of_step` promise (the midpoint relation) is not
this owner's: its ports answer `EndOfStep::Untold`, and a system does not join
one to a load line. `publishing_ends()` adds one signal output per power port
after `heat`, `<id>.end`, of width one, carrying that port's output at the end
state of the same native step, `y1 = G^T grad H(x1)`; the port then answers
`EndOfStep::Told` with it, and a load line that reads ends (fs-dcontact's
contact, fs-tribo's friction) may be joined to it. An end whose name repeats a
declared port (`a.end` beside `a`) is refused. Publishing changes neither the
step, nor the mean ports, nor the heat, bit for bit. The tangent's end rows are
`d y1 / d u` in the input direction at the call's fixed start state, not
derivatives in the start state: the step is affine in its held inputs and the
output is linear in the state, so column `j` is the end output of a step of
the same duration from rest under a unit input on port `j`, and each tangent
call takes those `m` extra native steps. The temperature's column is zero.

Caller state, outputs and derivative buffers publish only after complete
admission and the final cancellation check. Reentrant workspace use refuses.
Native controlled preparation and evaluation retain their bounded polling
boundaries; a storage callback and each documented indivisible dense/heat block
remain non-preemptible. The adapter retains the original dense limits; it
neither supplies a sparse path nor admits nonlinear storage by approximation.

`phs_exact` re-exports the fs-phs types its constructor takes,
`PortHamiltonian`, `QuadraticStorage` and `LinearSplitBudget`, so a caller
such as fs-cli's `port_hamiltonian` scene element assembles the owner without
a direct fs-phs dependency. The re-exports change no behaviour.

`tests/phs_exact.rs` covers bitwise agreement with the original exact owner,
Component and tangent conformance including transpose duality, independent
driven scalar exponentials and heat/heat-derivative integrals, variable-duration
thermal composition, exact-bit cache identity and bounded eviction, atomic
source/cancellation refusals, adaptive discard/restore/retry, native structural
and metadata refusals, and shifted-Hamiltonian endpoint evaluation allowances.
It also covers the published ends: their names and `Untold`/`Told` answers,
their equality with `G^T grad H(x1)` with the step and the means unchanged,
tangent conformance with the end rows, a scalar end and its input derivative
against the closed exponential, the name-collision refusal, and the refusal
of an exact owner beside fs-system's slider, which assumes the midpoint
relation: untold for that, told because the slider takes no end input.
fs-couple's `tests/drum_strike_graph.rs` joins a published end to fs-dcontact's
contact on the load line: two drum heads and the air between them as one exact
owner, struck by a stick, against fs-couple's prepared linear impact image with
its volume spring, each coordinate within -120 dB of its full scale.
These numerical checks do not establish physical validity of an authored model;
the other numerical adapters and sparse work in bead `.5.12` remain open.

## Kirchhoff circuit Component with heat ports

The same optional `phs` feature exposes `circuit::KirchhoffCircuit` (epic
`frankensim-rc-1109-root-3gc12`, the circuit/enthalpy coupling of
`fs-couple`). It owns an `fs_phs::circuit::CircuitGraph` and its signed mutual
inductances: any netlist that owner admits (resistors, inductors, capacitors,
voltage and current sources, ideal transformers) is one Component. Ports, in
order: `source<p>` for each external port, `EffortIn` when the port drives a
voltage source and `FlowIn` when it drives a current source, returning the
owner's step-mean conjugate; then `heat<b>` for each resistor, by branch
index, receiving temperature and returning minus the resistor's Joule power.
An external port must drive exactly one source branch, or construction
refuses. State is the owner's differential history only: inductor flux
linkages, then capacitor charges. Node potentials and constraint currents
start each step at zero, as `fs_couple::material_circuit` starts them.

A proposal assembles the circuit at the resistances of the step and runs the
owner's `CircuitDae::step_audited`; nothing is integrated here. A resistor
given a `ResistanceLaw` with `following(branch, law)` takes `law.ohms(T)` at
the temperature `T` its heat port receives; a law with no answer, or with a
resistance that is not positive and finite, refuses the step with the branch
index (`CircuitComponentError::Law`). It is never clamped. Each resistor's
heat is `dt v^2 / R`, `v` the voltage across it at the algebraic midpoint of
the step. Its entropy is that heat over the temperature its port received.
Supplied work minus stored change minus the resistors' heat is booked as
numerical loss and refused (`Balance`) beyond the arithmetic rounding floor
plus the caller's `energy_tolerance_j`. The receipt's `energy_tolerance_j` is
that rounding floor, as in `phs::PortHamiltonianComponent`; its residual and
tolerance are those of the owner's solve certificate, and a step without one
is refused. The circuit is `Class::Smooth` and declares no tangent. Identity
words bind the netlist, transformers, mutual couplings, each law's branch and
its words, counted, the tolerance, and the operating-point iteration's
controls when the circuit has one. The codec is one word of `f64` bits per
coordinate.

Operating point. A law may also follow the voltage across its resistor and the
current through it (`ResistanceLaw::follows_operating_point`, answered by
`ohms_at(T, v, i)`), as a conductor whose resistivity follows the field or the
current density does. Such a law needs the circuit's `OperatingPointIteration`,
given by `iterating(iteration)` before the law; `following` refuses it
otherwise, and the battery refuses it always. The circuit then solves the step
by the constitutive iteration of `fs_couple::material_circuit::conductivity`:
- every trial steps the owner from the same start, each such resistor at its
  law's resistance at an applied operating point, zero in the first trial;
- the point a trial solves is the step's midpoint voltage across the resistor
  and that voltage over the applied resistance;
- a trial is accepted when every such law, asked at the point it solved,
  answers within `absolute_resistance_tolerance_ohm` plus
  `relative_resistance_tolerance` times the larger resistance of the applied
  one, and the trial's power `v i` is within `power_tolerance_w` of the law's
  `R i^2`. The heat is the accepted trial's;
- otherwise the next trial applies the solved point, or with a `relaxation`
  below one the point that fraction of the way from the applied one;
- a circuit that has not met its tolerances by its `maximum_trials`-th trial
  refuses with `CircuitComponentError::OperatingPoint`.
A circuit with no such law takes one trial, whatever its iteration.

Which temperature a law sees is a property of the graph, not of this
Component. Joined by `common_effort` to an `elements::HeatCapacity`, it is the
body's step-mean temperature. Joined by a `heat_route` with
`RouteTemperature::SinkEquivalent` to a backward-Euler
`fs_thermal::component::ThermalTimeline` built `returning_end_temperatures`,
it is the cell's end-of-step temperature. That is the protocol of
`fs_couple::circuit_enthalpy`: backward Euler for the cells, circuit midpoint
currents, endpoint material resistance.

`tests/circuit.rs` (15 tests):
- the ports for a netlist with both kinds of source;
- a series RC discharging into two heat capacities: charge equal to the
  midpoint rule's `q0 ((1 - a) / (1 + a))^n`, heats split as the
  resistances, ledger closed;
- a voltage source driving R and L: current on the midpoint recurrence,
  source work equal to field energy plus body heat;
- an RC whose resistor follows the temperature of the body it heats, against
  a fixed point solved by bisection in the test, with the resistance moved
  by more than 5 %;
- a current source into a resistor: the port returns `I R` and the resistor
  takes `I^2 R`;
- `fs-phs`'s own reflected-load case (12 V, 100 ohm, a 3:1 ideal transformer,
  50 ohm), a circuit with no stored coordinate: primary current
  `U / (R1 + N^2 R2)`, each resistor's heat rate, the step's heat;
- `check_component` conformance with a source and a law;
- refusals: a law without an answer, by branch; a cold port; a port
  driving two sources; a port driving none; a law given to a capacitor or
  given twice; a zero tolerance;
- the receipt's heat and entropy at the received temperatures;
- codec round trip, a nonfinite word refused, and an identity that tells
  two laws apart by value and by branch and keeps apart laws whose words run
  together;
- a resistor following its field (`R = 1 + 0.05 |v|` in an RC), unrelaxed and
  relaxed by half, against the midpoint step solved by bisection in the test:
  within 2.1e-13 C and 1.5e-11 W on the export of 1e611fdf33 with series 1-8,
  under bounds of 1e-11 C and 1e-9 W. With the resistance tolerance at half
  the resistance, the power tolerance alone holds the step to 4.8e-11 C of
  it (bound 2.5e-9 C). One trial is refused;
- a steep current law (`R = 0.1 + 1000 |i|`): unrelaxed it is refused after
  64 trials, and relaxed by half it meets the closed form (charge to the bit,
  heat rates within 2e-14 W, bounds 1e-12);
- the iteration's admission: a following law before it, the iteration given
  twice, and each control out of its range are refused;
- a circuit whose law follows temperature alone steps bit for bit the same
  with an iteration as without; the iteration and its relaxation each change
  the identity; the iterated circuit passes `check_component`.
Ten defects planted one at a time in the iteration each failed at least one
of the last four tests: the relaxation ignored or reversed, the law asked at
the applied point, the point taken at the step's end, the power check
dropped, the first trial accepted, one trial past the cap, the iteration left
out of the identity, a following law admitted before the iteration, a
relaxation up to two admitted.
Fourteen defects planted one at a time in the Component each failed at least
one of these tests: heat at the end potentials, the law asked at 300 K, the
heat-rate sign, the source-output sign, a state not advanced, entropy at
300 K, every resistor at the first resistance, the stored-change sign, a law
left out of the identity, temperatures given to the wrong resistors, a port
driving no source admitted, a port driving two sources admitted, the law
words not counted, one resistor's heat rate left unwritten.

`fs-couple`'s `tests/circuit_enthalpy/graph.rs` builds the module's own
endpoint-feedback fixture (two material cards, one resistor rising 0.5 % per
kelvin, a 1 W/K link, eight intervals of 0.25 s) as such a graph. The
resistors follow the cards through the same `MaterialOhmicBody` resolution as
the module, and the comparison is against `circuit_enthalpy::step_candidate`.
On the export of 1e611fdf33 the cell temperatures and the charge agree bit
for bit over the first four intervals and the Joule heat within 9e-16 J;
after that all agree within 2e-9 K, 1.4e-12 C and 4e-12 J. The test's bounds are
1e-7 K, 1e-10 C and 2e-10 J. Six defects planted in the graph each failed
it: heat routed to the other cell, one resistor left at its placeholder
resistance, resistors at a fixed 300 K, cells returning step-mean
temperatures, the link conductance halved, a resistor's cell volume doubled
(the card then refuses the conductor's geometry). The same file declares the
module's mutual-coil fixture (2 H and 3 H coupled by 1 H, each shunted by a
card resistor rising 0.5 % and 0.2 % per kelvin, linkages 3.5 and 0.5 Wb, six
intervals of 0.25 s) with the same graph and only the netlist changed: it
agrees within 2e-9 K, 1.1e-11 Wb and 7.1e-12 J at every interval, under
bounds of 1e-7 K, 5e-10 Wb and 4e-10 J. The six planted defects above fail
it as well, and three more fail it alone: the coupling left out, its sign
flipped, the starting linkages swapped.

A third fixture is the module's independent-sources case: constant cards, a
2 W/K link, 2 W into cell 30 and 1 W out of cell 10. The graph needs no new
code for it. The sources become the cells' authored boundary sources, which
the timeline books on its `environment` port and the kernel routes to the
300 K sink. Over three intervals of 1 s the charge agrees bit for bit, the
temperatures within 3e-12 K and the Joule heat within 4e-15 J. The heat from
outside the circuit agrees exactly: 1 J an interval on both sides, the
module's `non_electrical_heat_j` against what the sink gave. The bounds are
1.5e-10 K, 1e-10 C and 2e-10 J. Four planted defects fail this fixture alone:
the graph's cells without their sources, the sources swapped between the
cells, the outside heat read with the sink's sign, and the module side without
its sources. Four of the six graph defects above also fail it; resistors at
300 K and step-mean cell temperatures do not, since its constant cards make
the resistance independent of temperature.

Entropy. Each interval of the three fixtures also holds the module's own
entropy account (`circuit_enthalpy::entropy::step_candidate`, which returns
the same candidate) to the kernel's ledger:
- its `electrical_loss_generation_j_k` to the circuit's entropy, which is
  each resistor's heat over the temperature its port received;
- its `conduction_generation_j_k` to the timeline's entropy, which is the
  owner's own now that fs-thermal's implicit timeline books it.

On the export of 1e611fdf33 the two sides agree within 1.7e-14, 6.6e-15 and
3.5e-16 J/K; the bounds are 1e-12, 3e-13 and 2e-14 J/K. The heat routes
carry only the interface iteration's residue: each route's heat times the
difference of the reciprocal temperatures that the port and the cell last
saw. The residue is at most 5.4e-16 J/K, 4.3e-15 of the total, and the
bound is 2e-13 of the total. Three defects planted in the comparison each
failed it: the module's two terms swapped, the kernel's two nodes swapped,
and the environment route at 290 K against the 300 K sink (the sources
fixture alone). Of the timeline-wiring defects listed with the timeline
slice, the three that change what the implicit step books also fail it: zero,
the storage change, and the receipt ignoring the solve. The one that adds the
reservoir input does not, since these fixtures have no film or reservoir.

Conductivity. `fs-couple`'s `tests/circuit_enthalpy_conductivity/graph.rs`
declares the fixtures of `circuit_enthalpy::conductivity` the same way. Each
resistor follows its card through a `ResistanceLaw` that resolves it as the
module does: the seed conductor at the phase state of the temperature its
port receives, as the fixtures' `coefficients` resolve it, then at the field
`V / L` or current density `I / A` of its operating point on each axis the
card binds, as the material circuit refreshes it. The circuit iterates with
the module's conductivity controls (64 trials, 1e-12 ohm, 1e-12 of the
resistance, the circuit's 1e-8 W). Both run on the backward-Euler timeline
returning end temperatures.
- The module's field fixture: resistivity `1 + 0.05 E` on resistor 1, resistor
  2 rising 0.5 % per kelvin, the 1 F capacitor at 10 C, a 1 W/K link, eight
  intervals of 0.2 s, unrelaxed. The two sides agree within 3.5e-9 K,
  9.2e-12 C and 9.9e-12 J at every interval, under bounds of 2e-7 K, 5e-10 C
  and 5e-10 J. Cell 10 warms 9 K.
- The module's steep current-density card, `0.1 + 1000 J`, with a constant
  heat capacity and resistor 2 as above, relaxed by half, six intervals. The
  charge agrees bit for bit, the temperatures within 1.8e-12 K and the Joule
  heat within 1.7e-16 J, under bounds of 1e-10 K, 5e-10 C and 1e-14 J.
  Unrelaxed, both sides refuse the first interval after 64 trials: the
  module with `ConductivityConvergence`, the kernel's circuit with
  `OperatingPoint`.
Of the ten iteration defects above, six fail these comparisons too: the
relaxation ignored or reversed, the law asked at the applied point, the end
point, the first trial accepted, one trial past the cap. Two defects planted
in the graph also fail them: the card's law at temperature alone, and the
kernel unrelaxed beside the relaxed module.

Not claimed: a tangent or `PortTangent`; second-order accuracy of the
coupled step, which is first order with endpoint coefficients as the module
documents; the module's storage, reservoir-input, temporal and algebraic
entropy terms, since the ledger books generation, not storage; the module's
phase-dependent resistance at a latent plateau, which the fixture's
single-phase cards do not reach; sources, conductances or reservoirs that the
module's coefficients resolve from the state of each trial (the comparison
uses constant sources only). For conductivity, also not claimed:
- the module's adaptive relaxation (`with_adaptive_relaxation`);
- series electrical contacts and their own contact laws;
- the module's relaxed fixture with a heat capacity that varies with
  temperature, which the comparison's timeline of one caloric segment does
  not hold;
- relaxation of a magnitude coordinate. The kernel relaxes the signed voltage
  and current, the module each bound coordinate, a magnitude on a magnitude
  axis. The two agree while a resistor's voltage keeps its sign from trial to
  trial, as in these fixtures.
The module is kept as the regression reference, and no `fs-couple` source
changes. Its test fixture gains `sourced`, which `initial` now calls with no
sources.

## Equivalent-circuit battery Component with heat ports

The optional `phs` feature also exposes `battery::EquivalentCircuitBattery`
(bead `frankensim-rc-1109-root-3gc12.9.2`, the battery coupling of
`fs_couple::circuit_enthalpy::battery`). It owns the battery of
`fs_phs::circuit::battery`: an OCV/SOC curve that owns the chemical energy, a
series resistance, one polarization RC pair and a bus capacitor. Ports, in
order:
- `terminal`, electrical and `FlowIn`: it receives the current entering the
  bus and returns the owner's discrete bus voltage;
- `heat_series` and `heat_polarization`, each receiving the temperature where
  that resistance's heat leaves and returning minus its heat over the step
  divided by the step.

State is the owner's `[chemical, polarization, bus]` charge in coulombs.
Stored energy is the owner's Hamiltonian, which no resistance changes.

A proposal builds the owner at the resistances of the step and runs its own
`PreparedBatterySource::step_controlled`, the componentwise discrete gradient,
exact across OCV knots; nothing is integrated here. A resistance given a
`circuit::ResistanceLaw` with `following(loss, law)` takes `law.ohms(T)` at
the temperature its heat port receives. A law with no answer, or with an
answer that is not positive and finite, refuses the step with the loss
(`BatteryComponentError::Law`); it never clamps. A law that follows the
operating point is refused when it is given: the battery asks its laws at
temperature alone.

**Heat and balance.** The heats and the terminal work are the owner's own
interval accounts. Work minus stored change minus heat is booked as numerical
loss, and that loss and the differences between these accounts and the owner's
native step record are each refused (`Balance`) beyond the arithmetic rounding
floor plus the caller's `energy_tolerance_j`. A refused step writes neither
outputs nor state. The receipt's residual and tolerance come from the owner's
solve certificate, and a step without one is refused.

**Class and identity.** The battery is `Class::Smooth`: the chemical energy
integrates a continuous piecewise-linear OCV, so it and its discrete gradient
are continuously differentiable. No tangent is declared. Identity words bind
the capacity, the knots, the parameters, the solve budget and its halvings,
each law's words (counted), and the tolerance. The codec is one word of `f64`
bits per charge.

`tests/battery.rs` (8 tests):
- the port order and kinds;
- a flat curve's stationary state under discharge and charge: only the
  chemical charge moves, each resistance takes `R I^2`, the terminal returns
  `U + I (R_s + R_p)`;
- both resistances following their temperatures, against the bus and
  polarization midpoint equations eliminated by hand (the module's own
  oracle);
- a step across an OCV knot whose stored change is the curve's closed-form
  chemical energy plus both capacitors';
- five steps on the kernel with a current source and two heat capacities:
  charges, both temperatures and the source's energy in closed form, and the
  ledger closed;
- `check_component` conformance with both laws;
- refusals that publish nothing: a law with no answer and one answering a
  negative resistance, each by its loss; a cold port; a zero step; a step
  that would empty the cell; a law given twice; a law that follows the
  operating point; a zero tolerance; a curve the owner refuses;
- the codec round trip, nonfinite and short words refused, and identities
  told apart by law, by loss and by curve.

Sixteen defects planted one at a time in the Component each failed at least
one of these tests:
- the law asked at 300 K;
- the laws given to the other loss;
- the port temperatures swapped;
- the heat-rate sign;
- the heat rates on the other port;
- the terminal returning the OCV;
- the state not advanced;
- the entropy at 300 K;
- the terminal current's sign;
- the law's answer ignored;
- a negative answer passed to the owner;
- the laws, and then the knots, left out of the identity;
- a law given twice admitted;
- a nonfinite state decoded;
- the heat booked as stored.
A seventeenth, the refusal of a law that follows the operating point removed,
fails the refusals test.

`fs-couple`'s `tests/battery_enthalpy/graph.rs` builds the module's three
stepping fixtures as graphs. Each graph is the battery, an
`elements::FlowSource` drawing the terminal current, the backward-Euler
timeline of `tests/circuit_enthalpy/kernel.rs`, and one `SinkEquivalent` heat
route per loss into the module's receiving cell. A cell that receives no loss
gets a zero heat boundary. Resistances that the module resolves from the
receiving cell's card follow that card through the same state-point
resolution. On the export of 1e611fdf33:
- **The knot fixture** (no current, one step of 30 ms across the 0.5 SOC knot)
  agrees bit for bit.
- **The card fixture** (both resistances rising 0.4 % per kelvin of their
  cell's end temperature, 1 A drawn for four steps of 0.2 s, cell 30 warming
  by more than 15 K): the first step agrees bit for bit, and every step
  within 3.9e-12 K, 4.7e-15 C and 5.6e-16 J.
- **The common-receiver fixture** (both losses into cell 30, both signs of the
  current) agrees bit for bit except one rounding of 5.7e-14 K.

The bounds are 2e-10 K, 2.5e-13 C and 3e-14 J for the card fixture and
1e-11 K, 1e-13 C and 1e-13 J for the others. The kernel's ledger closes, and
the cells hold the heat the battery released. Eight defects planted in the
graph each failed it:
- the losses routed to the other cell;
- a law from the other cell's card;
- the resistances not following their cards;
- cells returning step-mean temperatures;
- the load current's sign;
- the card law asked at 300 K;
- the unused cell heated by 1 W;
- the module side at fixed resistances.

The kernel side shared by both comparisons is in
`tests/circuit_enthalpy/kernel.rs`: the cell timeline, its conduction law and
the execution context, moved there from `graph.rs` unchanged.

Not claimed:
- OCV that depends on temperature, reversible entropic heat, capacitance
  drift, ageing, or protection events;
- second-order accuracy of the coupled step, which is first order with
  endpoint resistances, as the module documents.

The module stays as the regression reference, and no `fs-couple` source
changes. The module's fourth test,
`unsupported_sources_exhaustion_events_and_cancellation_publish_neither_history`,
fails on the unmodified export of 1e611fdf33 (a card query outside its
claim's domain at line 530) and is not touched here.

## Signal-controlled power sources

`elements::ControlledSource::effort(kind)` and `flow(kind)` connect controller
signals to an ideal physical source. Port `p` receives the conjugate value and
returns the scalar `command` signal in that power domain's SI effort or flow
unit. The command is held over the step. Both signs are accepted, with every
physical port value positive into its owner: a negative electrical flow command
supplies positive current to a connected capacitor. A command signal carries
zero power; all exchanged work enters the existing `SourceAccount`, negative
when the source delivers energy. Its compensated clock advances with the step.

The complete port derivative is `[0, 1]` in input order `[conjugate, command]`;
JVP and VJP preserve that signal-to-power response. Wrong extents, nonfinite
inputs, invalid steps and unrepresentable work/accounts refuse before publishing
state or outputs. Proposals and both tangent actions poll cancellation before
work and before publication, allocate nothing, and retain deterministic order.

`tests/controlled_sources.rs` checks Component/tangent conformance in four power
domains, signed capacitor charging and recovery against independent charge/work
equations, and an actual Mass driven and braked by a signal feedback component
against its independent discrete recurrence and the full energy ledger. Invalid
inputs and cancellation preserve caller state, port and derivative buffers.
Sampling, event scheduling, control policy, gains, limits, delay, source
impedance and finite supply physics belong to the caller or their own
Components. This ideal source makes no control stability or hardware claim.

## Reusable signal feedback

`elements::AffineSignal::new(n, m, matrix, bias)` is the stateless authored
map `y = bias + matrix * input`, with scalar signal ports `input0` through
`input{n-1}` and one `m`-wide `output`. Both dimensions are positive, matrix
and bias extents are exact, and every coefficient is finite. It is
`ExactLinear`, with analytic JVP `A` and VJP `A^T`, zero physical energy and
dissipation, and no per-step allocation. Its native stateless codec binds
dimensions and every coefficient/bias bit into the problem identity.
Arithmetic or shape errors and cancellation refuse; the kernel publishes
none of a failed coupled trial.

`Mass::publishing()` adds accepted, left-end `position` and `velocity`
signals without changing motion, power outputs or receipts. Their port-input
derivatives are zero at fixed accepted state. The same observations apply to
the inertance's chosen domain: position/velocity, angle/angular velocity,
or passed charge/current. Repeated publication requests add no duplicate ports.

`tests/affine_signal.rs` checks rectangular maps and transpose actions,
unchanged Mass stepping, independent PD feedback dynamics and source work,
rejected/cancelled steps and bitwise checkpoint replay. The existing durable
snapshot tests reject changed gains, bias or published-port shape. Gains and
signal units are authored; this map infers no saturation, delay, sample period,
dynamic controller state or stability. The `feedback_control` example builds
a complete body/controller/actuator graph from these production components.

## Sampled control and replayable references

`SampledController` implements the authored discrete state-space law
`x_next = A*x + B*u`, `y = C*x + D*u`. It samples at time zero and at integer
multiples of its positive period. The output uses the pre-update memory, then
is held through physical substeps until the next sample. Scalar `input0`,
`input1`, etc. and a vector `output` use the same wiring as `AffineSignal`;
zero memory coordinates admit a sampled gain. Matrices are row-major and their
exact finite extents are checked at construction. State, held output, sample
count and the compensated clock are accepted together and are restored on
rejection or checkpoint replay. No state update is committed by a trial alone.

Sample boundaries constrain the kernel's event schedule. Indexed sample times
avoid repeatedly adding the period; loss of representable time resolution is
an explicit refusal. The input tangent is `D` when sampling and zero while
holding, with the matching transpose; the component is `Smooth`, so the
kernel does not cache a state-independent `ExactLinear` tangent across phases.
It stores no physical energy: actuator work belongs to `ControlledSource` and
the connected plant. Scratch is allocated at construction, separate from the
native replay state and problem identity.

`SignalSource` produces a nonempty vector of authored `Waveform` point values
at the accepted start time, independent of proposed step duration. This makes
a sampled controller's reference independent of physical subcycling. Table
knots are declared events; sine references leave temporal resolution to the
authored schedule. Its compensated `SignalClock` has a native codec, and every
waveform is bound into the problem identity. These are zero-power signals;
gains, units, filtering, saturation and stability remain authored choices.

`tests/sampled_controller.rs` compares controller recurrence, sample-and-hold
behavior and physical capacitor work against independent calculations, checks
analytic transpose actions, and exercises refusal, cancellation and replay.
`tests/signal_source.rs` checks point sampling, table event clipping, codec
identity and refusal without publication. It also checks held-command jumps,
interval means, distinct interpolation identities and a switched force's
trajectory and physical work against piecewise constant-acceleration formulas,
including checkpoint replay. `feedback_control --sampled-pi`
composes a replayable reference, a 50 ms PI controller, a measured inertance
and a controlled actuator while stepping the physical system at most 10 ms.

## Executable general Dirac junctions

`elements::DiracJunction::new(schema)` executes an admitted scalar or
Euclidean-vector `fs_ports::JunctionLaw::Dirac` relation as a stateless,
exact-linear Component. Its ports retain the schema's identifiers and order,
receive effort, and return `J * effort`. For a common vector width `d`, the
relation is `J ⊗ I_d`: each coordinate is coupled independently across the
ports. Runtime arrays are port-major, with a port's `d` coordinates contiguous.
The schema endpoints denote attached owners: their outward flows enter the
junction, and ordinary bonds negate those flows into the owners. No matrix
inversion, storage, dissipation or integration rule is introduced;
odd-dimensional, singular and zero skew matrices are valid.

The existing `JunctionSchema` owns structural admission, including equal
endpoint shape, basis, frame, orientation, clock and power pairing. The adapter
refuses other junction laws, tensor and field endpoints, and an unrepresentable
packed extent. JVP applies `J ⊗ I_d`; VJP applies `-J ⊗ I_d`. Finite vector
extents, positive steps and representable arithmetic are checked, and
evaluation polls cancellation. The energy allowance counts the magnitudes of
individual work terms across every coordinate before cancellation. Scalar
evaluation retains its original arithmetic order and problem identity; vector
shape and coordinates are already bound by the retained schema.

`tests/dirac_junction.rs` covers exact scalar and vector relations, power
cancellation, tangent duality, and closed storage networks against independent
midpoint-rotation oracles. It also covers malformed source laws, unsupported
shapes, invalid calls and cancellation. The retained schema exposes the common
coordinate/frame/clock declarations, but runtime `PortSpec` power ports do not
enforce those labels: attached owners must already use the same coordinates.
This relation supplies no field metric, coordinate conversion, material law or
hardware validity claim.

## Signal-modulated lossless transformer

`elements::ModulatedTransformer::new(kind_a, kind_b)` adds a scalar transformer
whose ratio is a signal from a controller, geometry or constitutive owner.
Port `a` receives flow and returns effort; `b` receives effort and returns
flow. Inputs are `[f_a, e_b, r]`, outputs `[r e_b, -r f_a]`, with all flows
positive into the transformer. The same held ratio multiplies both sides, so
power cancels to arithmetic roundoff even as it changes sign or passes through
zero. The ratio signal carries no energy. A constant nonzero ratio reproduces
the existing fixed transformer between the attached owners. Domains may differ;
the caller supplies the ratio in the required SI conversion units.

The complete map is bilinear (`Class::Smooth`); its analytic JVP and VJP include
the ratio column as well as both physical columns. The component stores no
energy and dissipates none. Invalid extents, nonfinite inputs, nonpositive or
nonfinite steps, and unrepresentable outputs or work refuse before publishing
caller buffers. Proposals and tangent actions allocate nothing and poll
cancellation before work and before publication. Work roundoff is bounded from
the two exchange magnitudes before cancellation.

`tests/modulated_transformer.rs` checks Component and tangent conformance,
signed and zero ratios, agreement with the fixed transformer, a varying-ratio
mass/capacitor graph against its independent midpoint rotation, and a ratio
fed back from accepted storage state against an independently eliminated
recurrence. Both graphs close their energy ledger. Malformed, overflowing and
cancelled calls leave caller buffers unchanged.

The ratio is held over each step; its source owns sampling and event resolution.
No crank/cam geometry, displacement constraint, ratio-actuator work, physical
material validity, or control stability is inferred. Two integral-causality
inertias still need a constraint solve or lowering. Only a mean flow is returned
at `b` (`EndOfStep::Untold`), so endpoint-dependent load lines are not admitted
there. This component advances the variable-ratio part of bead `.5.3`; the
other connection and kinematics work remains open.

## Scheduled boundary events

`SystemBuilder::scheduled_boundary(port, value, event_times_s)` lets any
prescribed port declare the absolute times [s] at which its input law may
change. A schedule is finite and strictly increasing; an invalid entry is
`AdmissionError::BoundarySchedule`, naming the port and entry, and the builder
claims no port. Empty schedules retain the behavior of `boundary(port, value)`.
The callback still fills the values for the admitted interval, so this works
for power, heat, signal, field and stream inputs and for endpoint spans.

The ordinary step driver finds the earliest future boundary event after the
component limits and events have been checked. Events at or before the accepted
clock are ignored. Its interval retains the clock's compensation and uses the
same near-knot cap policy as `Force`, `FlowSource` and `Heater`: a strict cap
that would strand only clock rounding leaves half the remaining interval.
Distinct adjacent representable event times remain distinct. Lookup allocates
nothing and advances no cursor. Failed proposals, checkpoint restore and
adaptive coarse/half-step trials therefore see the same immutable schedule.
Scheduled times participate in the durable problem identity; existing graphs
with no scheduled boundaries retain their previous identities. Boundary
callback captures remain external to snapshots and must reproduce on replay.

`Waveform::event_times_s` exposes all knots of a table with at least two rows,
including its transitions to and from held endpoint values. A held table
exposes every knot after the first, where its value may jump. Constants,
sines and one-knot tables expose none. A vector input uses the sorted, deduplicated
union of its waveforms' event times. Raw unscheduled callbacks remain supported;
the kernel does not infer their events, and knowing event times supplies no
accuracy bound between them or guarantee that the declared schedule is complete.

`tests/system_g0.rs` exercises coarse requests across two boundary schedules,
independent ramp-force integrals and input work, refusal/retry and exact
checkpoint replay, compensated near-knot steps and adjacent times, atomic
schedule admission and problem identity. `tests/adaptive.rs` exercises the same
negotiation through step doubling, an invalid-estimate rewind and independently
integrated heat input. The waveform-event assertions remain in `tests/elements.rs`.

## Ideal electrical diode and commanded switch

`elements::IdealDiode` and `elements::IdealSwitch` are stateless, lossless
`LoadLine` components, added with `SystemBuilder::add_on_line`. Their scalar
electrical port `p` receives current [A] and returns voltage [V], with positive
power entering the device. For the diode, current runs from anode to cathode
and voltage is anode minus cathode: the ideal relation is `i >= 0`, `v <= 0`,
`i * v = 0`. Forward conduction has zero voltage; reverse bias has zero current.
The switch's scalar signal `closed` must be exactly `0` or `1`: open enforces
zero current, closed enforces zero voltage and admits either current direction.
Fractional commands refuse. Prescribe changes with a scheduled boundary or
`SignalSource` using `Waveform::held_table`, so steps respect authored events.

The kernel supplies `i = released - C * v`, with finite nonnegative compliance
`C`. An open device with `C = 0` and nonzero prescribed current refuses; with
both zero it chooses zero voltage without asserting uniqueness. Neither device
stores energy, dissipates heat or substitutes a large resistance. Both have
native empty-state codecs with distinct identities. Invalid inputs, impossible
current constraints and cancellation refuse without publishing outputs.

These are ideal relations between a step's representative values. They do not
locate circuit-induced zero crossings, model semiconductor forward drops,
recovery or breakdown, or resolve arcing and ideal-storage switching impulses.
The caller owns temporal resolution and event schedules; a binary command is
not a continuous control law for instantaneous feedback inside a nonlinear solve.

Run `cargo run -p fs-system --example electrical_switching` for two RC branches
whose resistors heat a shared thermal capacity. One diode blocks discharge
when its supply reverses; a commanded switch interrupts the other branch's
charging. The example emits capacitor voltages, branch currents, temperature,
source work and the energy residual from the composition kernel. The seven
regressions in `tests/electrical_switches.rs` cover constitutive and energy
conformance, analytic resistive circuits, event boundaries, refusal and replay,
native codecs, and independent midpoint RC charge and thermal-work balances.

## Signal-controlled stream restriction

`Restriction::modulated()` adds a scalar, zero-power `opening` signal to
any existing anonymous, identified or gas-aware restriction. The two stream
ports and their exact material/source, caloric-reference and carried-quantity
identities are retained. Inputs are the two native stream value vectors followed
by opening; outputs remain the two native mass/enthalpy/carried-rate vectors.
A finite opening in [0, 1] multiplies the authored conductance. Zero closes
transport exactly and one reproduces the ordinary law. The original donor rule
and, for `for_gas`, its original throttle-entropy calculation are retained.
There is no mass, energy or species storage and no alternate stepping scheme.

The analytic JVP and VJP include the command derivative as well as both stream
inputs. At zero opening the donor for the right derivative is selected by the
pressure difference, including reverse pressure; at equal pressures the existing
donor kink remains. Endpoint derivatives in opening are one-sided. Commands
outside [0, 1], wrong vector extents, nonfinite inputs/directions, invalid steps
and unrepresentable rates refuse. Calls poll cancellation before work and while
checking carried data. Evaluation allocates no storage. The empty state codec
is retained; a mode tag distinguishes a modulated restriction, while ordinary
restriction identities and valid-input arithmetic are retained.

`tests/modulated_restriction.rs` adds conformance and directional/transpose
checks, full-open parity, closed-valve one-sided derivatives, source-basis and
codec checks, refusals, and a scheduled two-vessel graph with conserved mass,
energy and species and checkpoint replay. These authored tests are numerical
evidence only when executed on the corresponding revision. The signal source
owns sampling and event resolution. This is an ideal linear valve; it supplies
no actuator work, choking, leakage, geometric valve law or experimental claim.

## Explicit signal limits

`elements::SignalLimiter::new(lower, upper)` is the stateless, coordinatewise
control law `output[i] = clamp(input[i], lower[i], upper[i])`. Its vector
`input` and `output` ports have the same nonzero width. Both bound vectors
are finite and equally sized, with `lower[i] <= upper[i]`; equal bounds hold
a coordinate constant. The bounds are authored in the signal coordinates'
units. This is an explicit part of the controller, never silent clamping of
physical state or a waiver of an actuator's command admission rules.

The component is `NonSmooth`, has an empty native state codec, and binds every
bound bit and its width into the durable problem identity. It stores no energy
or mass and exchanges no power. Add it with `add_with_tangent` in a feedback
loop: the ordinary kernel solves the same graph through its declared diagonal
tangent. Entries are one strictly inside the bounds and zero outside. At a
bound the tangent selects zero from the generalized Jacobian, identically in
JVP and VJP; no classical derivative is claimed at a nondegenerate kink.
Validation refuses malformed extents, nonfinite inputs or directions, invalid
steps and cancellation. Finite input magnitudes outside the interval are
accepted. Evaluation allocates nothing and polls at bounded vector intervals.

`tests/signal_limiter.rs` checks component conformance and finite-difference
tangents away from kinks, the selected tangent at bounds, transpose duality,
fixed coordinates, native codec identity, refusal and cancellation. An
algebraic feedback graph changes between both saturation branches and its
interior while driving a real mass through a controlled effort source; its
motion and work follow independently solved forces, and checkpoint replay
crosses the active-set change bitwise. This law supplies no anti-windup,
smoothing, time-of-crossing search or feedback-stability guarantee.

## Field in a gap of variable width

`elements::GapField` stores `H = x^2 g / (2 s)` in a gap of width `g` [m]:
- the field of a parallel-plate capacitor (`GapField::capacitor`), with `x`
  the charge [C] and `s = epsilon A` [F m];
- the air gap of a winding whose return path stores no field
  (`GapField::inductor`), with `x` the flux linkage [Wb] and `s = mu N^2 A`
  [H m].

Its state is `GapFieldState { coordinate, gap_m }`. Port `field` receives the
rate of `x` and returns its conjugate. A capacitor's port is `FlowIn` (current
into the positive plate); an inductor's is `EffortIn` (voltage across the
winding). Port `gap` (`FlowIn`) receives the rate at which the gap opens and
returns the force the field holds against that opening, positive as the
plates or poles attract.

Both returned efforts are the exact mean of the gradient along the segment
between the step's end points: the average-vector-field discrete gradient,
the formula of fs-couple's moving field.

```text
e_x = x_m g_m / s + (dx / 12) (dg / s)
e_g = (x_m^2 + dx^2 / 12) / (2 s)
```

The receipt is lossless. Its stored change `e_x dx + e_g dg` is
`H(end) - H(start)` within rounding, for every step. The class is `Smooth`,
and the port tangent is the exact, symmetric derivative of the two efforts.

The gap is admitted strictly inside `gap_limits_m`, which must be finite with
`0 <= lower < upper`. A state or a step end outside that range is refused
(`ElementError::Law`), never clamped. The coupling scale and its reciprocal
must be positive and finite. The scale is constant: a coefficient that changed
with temperature would need the field's entropy, which the element does not
hold. The element has no fringing, saturation or leakage, and no motion of its
own. A `Mass` carries the moving boundary, joined to `gap` by a common-flow
junction with ratio -1, so that the mass's velocity opens the gap.

`tests/gap_field.rs` (4 tests):
- the conformance kit and the tangent check, over 64 states and inputs of
  both kinds. The efforts are held to Simpson's rule on the gradient (which
  is quadratic along the segment) within 16 eps of its size, and the stored
  change to `H(end) - H(start)` within 32 eps;
- a held field (no current into the plates, no voltage across the winding)
  joined to a mass and a held force, over ten steps:
  - the coordinate is held bit for bit;
  - the gap moves bit for bit with the mass;
  - both follow the uniform acceleration `(F - x^2 / (2 s)) / m` within
    1e-14 m and 1e-15 m/s;
  - the force's work is the field's and the mass's energy change;
- a gap discharging through a resistance (an inductor's `Damper`, a
  capacitor's `Admittance`) while a mass and a force move it, against the
  monolithic discrete-gradient step solved by hand with Newton:
  - twenty steps agree within 4.5e-15 (bounds 1e-13);
  - the resistance's heat reaches the heat capacity within 1.3e-14 J;
  - the work balance closes;
- refusals:
  - the scale and the limits at construction;
  - states outside the range (energy, codec);
  - a step end outside it, and a non-finite input, with the outputs and the
    candidate untouched;
  - the kind, the scale and each limit change the identity.

Nine defects planted one at a time in the element each failed at least one
of these tests, as did two in the graph below:
- in the element:
  - either effort without its increment term;
  - either tangent term wrong;
  - a closed lower limit;
  - the identity without a limit;
  - the inductor's port receiving flow;
  - a decode admitting any gap;
  - the gap advanced at the field's rate;
- in the graph:
  - the mass closing the gap;
  - the force reversed.

fs-couple's `tests/circuit_enthalpy_motion/graph.rs` declares the moving
field of `circuit_enthalpy::motion::step_candidate` as a graph:
- the module's netlist in a `KirchhoffCircuit`, with the gap's branch
  replaced by a source port (a voltage source for a capacitor, a current
  source for an inductor), and each resistor following its card at the
  temperature its heat port receives;
- a `GapField` bonded to that port;
- a `Mass` and a `Force` in one common-flow junction with the gap;
- the backward-Euler timeline and one heat route per resistor.

Against the module, at every interval:
- a gap inductor shorted through two heated resistors (8 x 0.1 s): within
  1.2e-9 K, 1.9e-12 Wb, 3.4e-15 m, 1.0e-14 m/s and 1.9e-12 J;
- a gap capacitor pushed open while it discharges (8 x 0.1 s): within
  2.7e-11 K, 1.4e-14 C, 1.3e-14 m, 1.8e-14 m/s and 1.1e-13 J;
- a linkage held by a 0 V port beside a driven branch (4 x 0.25 s): the
  linkage, the gap and the velocity are equal, and the sides agree within
  1.3e-8 K and 7.0e-11 J. Both meet the uniform acceleration within 1e-10.

The graph's work balance closes within 1e-8 J: the force's and the ports'
supply against the gap's, the mass's and the cells' energy change.

Not claimed:
- the module's layered axial support;
- its coupled multi-winding gap;
- its nonlinear dielectric;
- its conductivity and contact variants;
- its gap-limit and budget refusals as graph behaviour;
- a field coefficient that depends on temperature;
- fringing.

## Commutated battery-motor drive Component with heat ports

The optional `phs` feature also exposes
`commutated_drive::CommutatedBatteryMotorDrive` (bead
`frankensim-rc-1109-root-3gc12.6.6`, under `.9.2`). It is the drive that
`fs_couple::circuit_enthalpy::converter::commutation::battery_motor` coupled
by hand in its `capacitance` and `charge` lanes, which the graphs below
replaced. Its owner is the drive of
`fs_phs::circuit::converter::commutation::capacitance::battery_motor`: the
equivalent-circuit battery and its lead, a two-leg converter of four switches
(channel, body diode and output storage each) and a DC motor whose rotor and
stator both turn.

**Lanes.** The drive's type parameter, a sealed `DriveLane`, fixes the
switches' storage, the state and the owner the drive prepares:
- `DeviceStorage`, the default: the drive above, eleven coordinates. Code
  that names no lane is in this one;
- `NoDeviceStorage`: switches without output storage, the owner
  `PreparedCommutatedBatteryMotorCycle` of
  `fs_phs::circuit::converter::commutation::battery_motor`. Its nine
  coordinates are the battery's three, the converter's flux linkage and two
  bus charges, and the motor's three. It is the drive that
  `fs_couple::circuit_enthalpy::converter::commutation::battery_motor`
  itself couples by hand;
- `GateStorage`: switches with driven gates, the owner
  `PreparedGateBatteryMotorCycle`. Each switch keeps an intrinsic
  drain-source charge curve, a gate-source and a gate-drain capacitance and
  a `ChannelLaw`; each gate a driver with on and off resistances and a
  current limit. Either resistance may follow a `ResistanceLaw`
  (`driver_following(switch, rail, law)`). Its fifteen coordinates are the eleven's, with each node
  charge grouped with its high-side gate electrode, then four gate charges.
  The rail voltages are not in the lane: each step takes them from the rail
  ports, so a source holding a rail and the drive cannot disagree.

Ports, laws, programs, heat, balance and refusals are the same in every
lane; what follows holds for all unless it names one.

A step is the owner's whole gate program, `step_interval_commands_controlled`
over every interval and substep in one native workspace; nothing is
integrated here. The owner is built with the resistances and switches of the
step. A step is one period of the program; any other length is refused.

Ports, in order:
- `rotor`, `stator`: rotational, `EffortIn`. Each receives the torque on its
  shaft, held through the program, and returns the owner's time-weighted mean
  shaft velocity. Torque times mean velocity times the step is the owner's
  shaft work.
- In the gate lane, the eight rails of `RAIL_PORTS`, electrical, `EffortIn`:
  each switch's on rail, then its off rail. Each receives its voltage, held
  through the program, and returns the mean current its gate drew from it,
  from fs-phs's per-rail charge (`GatePowerAccount::rail_charge_c`). Voltage
  times mean current times the step is the rail's work exactly; one port per
  gate at its mean voltage and current would not be, since a gate switches
  rails within the program.
- The ten heat ports of `HEAT_PORTS`: the converter winding; the four
  switches, each with its channel and diode heat together; the battery's two
  resistances; the lead; the motor winding; the bearing. In the gate lane
  the four of `DRIVER_HEAT_PORTS` follow, each gate's on and off driver
  resistances together.
- `heat_auxiliary_load` when the program has an auxiliary load: the work the
  load took.

Each heat port receives the temperature where its heat leaves and returns
minus that heat over the step divided by the step.

State is the owner's coordinates, eleven, nine or fifteen by lane. Stored
energy is the owner's Hamiltonian: the battery's, the converter's (switch storage
included) and the motor's. Without switch storage the converter's is its
averaged bridge's: the inductor's and the two buses'. The switch storage is
one `DeviceStorage` parameter: four constant capacitances or four
charge-voltage curves. It is fixed; storage that changed with temperature
would need its entropy. fs-phs builds the gate storage from a whole
parameter set, though neither the channel nor the drivers enter its energy.
The gate lane builds it with a stand-in channel that refuses every
evaluation and the drivers at zero rail voltage, and admits each program on
it when built; no step runs on it.

**Laws.** `following(resistor, law)` gives one of the five resistances a
`circuit::ResistanceLaw`. `switch_following(position, law)` gives a switch a
`SwitchLaw`: its channel, diode and limits at a temperature. Each is asked at
the temperature its heat port receives. A law with no answer, or a resistance
that is not positive and finite, refuses the step with the port's name
(`CommutatedDriveError::Law`); nothing is clamped. Inductances, capacitances,
the OCV, the motor's coupling and inertias and the bearing drag are fixed.
In the gate lane each switch's `ChannelLaw` is asked at its heat port's
temperature too; the `SwitchLaw` there gives only the body diode and the
limits, since fs-phs takes the channel from the response. A driver
resistance following a law is given its answer at its gate's driver heat
port; a law that follows the operating point is refused there too.
The error type is the drive's own, as for the battery: `ElementError` has no
variant for an owner's refusal or a balance defect.

**Programs.** The drive runs the program it was built with every step unless
`scheduled(programs)` gives it further ones. It then has a signal input
`program` after its heat ports, and the load port if any program has a load.
The signal names the step's program: `0` the first, `k` the `k`-th further
one. It must be an exact whole number naming a program, as `IdealSwitch`
takes `closed`; anything else refuses the step. A `SignalSource` of
`Waveform::held_table` schedules it. Every further program has an auxiliary
load entry per interval, is admitted by the owner and has the drive's period;
the drive takes its further programs once.

**Heat and balance.** The heats, the load's work and the shaft works are the
owner's cycle accounts. Each of the following is refused (`Balance`) beyond
the rounding floor plus the caller's `energy_tolerance_j`:
- the owner's balance residual, booked as numerical loss;
- the differences of stored change, heat and supplied work from the owner's
  native record;
- the difference of each shaft port's work from the owner's shaft work;
- in the gate lane, the difference of each gate's two rail works from the
  owner's source work for that gate. The supplied work there adds the rails'
  works, and the heat the drivers'.

A refused step writes neither outputs nor state. The owner keeps each
substep's Newton residual within its budget and gives the cycle no
certificate, so the receipt reports the largest substep residual as both
residual and tolerance.

**Class and identity.** The drive is `Class::Smooth` with no tangent: the
kernel differences the port map. The map is smooth except where a body diode
starts or stops conducting, a switch reaches a limit or the load changes
branch as the inputs move. Identity words bind the curve, every parameter,
the storage law (a tag, 0 for capacitances, 1 for charge curves, 2 for
none and 3 for driven gates, then its words; for gates the chord budget,
the curves, both capacitances, the drivers, and each channel and driver
law's words, counted), the program, both controls, each law's words (counted) and the
tolerance, then any further programs. The codec is one
word of `f64` bits per coordinate; encode and decode refuse a state whose
energy the laws do not admit. Each lane has its own state charter.

`tests/commutated_drive.rs` (5 tests with both storage laws each, and one of
the lane without switch storage), on a program
of a dead time then input low and output high on, with an auxiliary load;
both the input low and output high diodes conduct:
- a step with all five resistances and four switches following laws against
  the owner built by hand at the laws' answers: state, mean velocities and
  every heat port bit for bit; stored change, entropy; the conformance kit;
- six periods on the kernel with two torque sources and eleven heat
  capacities. Proposed again at the temperatures the bodies returned, the
  drive takes the kernel's step within 1.5e-16 relative (bound 1.5e-15) and
  gives each body its heat within 2.6e-14 J (bound 2.5e-13). The ledger
  closes within 6.2e-18 J (bound 1e-16), and the input high body warms
  past 330 K, so its diode threshold feeds back;
- refusals:
  - construction;
  - a second law on the same resistance or switch;
  - a step that is not one period;
  - a non-finite torque or temperature, and a zero temperature;
  - a switch law out of its range (`Law` naming `heat_input_low`);
  - the owner's channel limit (`Step`);
  - cancellation.

  Each leaves outputs and candidate untouched. The codec refuses a short
  word list and a state outside the curve, and the identity changes with
  the storage law and with every law followed;
- a drive given a further program, with only the output high switch on after
  the dead time. Signalled `0` or `1`, it takes bit for bit the step of the
  drive built with that program alone, and the two steps differ. `0.5`,
  `2` and `-1` are refused, as are no programs, a program of another
  period, a program short of a load entry and a second list. The
  conformance kit passes with the signal. The identity changes with the
  further programs;
- the lane without switch storage, from the fixture's cell, motor and bus
  voltages, against `PreparedCommutatedBatteryMotorCycle` built by hand at
  the laws' answers: state, mean velocities and every heat port bit for
  bit. Its stored energy is the owner's whole Hamiltonian (measured equal,
  bound 1e-14 relative). The receipt's stored change is within 8.4e-14 J
  of the energies' (bound: the receipt's tolerance, 2.0e-10 J). The
  conformance kit passes. A step that is not one period, and a switch law
  out of its range (`Law` naming `heat_input_low`), are refused, leaving
  outputs and candidate untouched. The codec keeps nine words and refuses
  the eleven of a stored drive's state. The identity parts from a
  capacitance drive's at the storage tag, 2 against 0;
- in `tests/commutated_drive/gate.rs`, the gate lane equals
  `PreparedGateBatteryMotorCycle` built by hand at the laws' answers and
  the rails' voltages, bit for bit, with Miller coupling, channel laws that
  vary with temperature and six of the eight driver resistances following
  laws; each channel law is asked once, at its switch's heat port. It
  refuses a non-finite rail, a zero driver temperature, a rail past its
  driver's current limit (`Step`), and a channel or driver law with no
  answer or a negative one (`Law` naming the heat port), leaving outputs
  and candidate untouched. `driver_following` refuses a second law and one
  that follows the operating point. Its identity changes with the Miller
  capacitances, a driver resistance, the channel laws and the driver laws;
- without Miller coupling and with channels that conduct nothing, each
  gate's voltage, rail charges and driver heat meet the closed form of the
  owner's implicit midpoint at its driver laws' answers, (2RC - h)/(2RC + h) per substep, within 1e-12
  relative. Every other coordinate and port equals, bit for bit, that of a
  program whose gates never switch;
- on the kernel, with a source on every rail and shaft, the ledger closes
  within 1e-15 J.

Planted in the gate lane with the comparison to the owner by hand switched
off, one at a time, each failed at least one gate test: each gate's rails
swapped; each channel asked at the next switch's port; the driver heat
ports given nothing, or in reverse order; each rail returning its charge
for its mean current; the rails' work left out of the supplied work; the
driver heats left out of the drive's heat; each rail's work at the other
rail's voltage; the identity without the Miller capacitances, the drivers'
resistances or the channel laws; the Miller capacitances taken from the
gate-source ones, which fails the closed form and, with the closed form
switched off, the run whose gates never switch; each driver law asked at
the next gate's driver port; the on and off driver laws swapped; a driver
law's answer ignored; the identity without the driver laws; a second
driver law replacing the first; a driver law following the operating point
accepted; a driver law's refusal naming its switch's heat port; a negative
driver resistance passed to the owner. Not carried from fs-phs's
gate tests: the channel's refusals of a gate voltage outside its support
or a current past its limit, and final cancellation.

Planted defects, one at a time, each failed at least one test:
- a switch's heat taken with another switch's diode;
- the load's work not routed;
- the stator given the rotor's torque;
- two resistance laws asked at each other's ports;
- a switch law asked at its neighbour's temperature;
- a step of any length;
- a decode that admits any state;
- a switch law with no answer replaced by the authored switch;
- an identity without the resistance laws;
- the program's end velocity in place of its mean;
- entropy without the load;
- the program signal ignored;
- a fractional program signal rounded down;
- the identity without the further programs;
- a further program of any period;
- in the lane without switch storage:
  - the converter's energy with its two buses swapped;
  - a decode that takes a word list of any length;
  - the energy split where the stored lane's motor starts;
  - the owner prepared with one more substep;
  - the storage tagged as the capacitances';
  - the stored lane's charter.

fs-couple's `tests/commutated_converter_enthalpy/battery_motor/graph.rs`
declares the drive of the retired `capacitance` and `charge` lanes ("the
module") as a graph:
- this drive with the fixture's storage, its five resistances and four
  switches following the fixture's card at the temperature each heat port
  receives. A switch's law also asks its package's storage claim, so a
  temperature outside that claim is refused as the module refused it;
- two `Force` sources holding the command's shaft torques;
- the seven cells as the backward-Euler timeline returning end-of-step
  temperatures, and one heat route per heat port;
- where the module links two cells, the same link in the timeline. A
  test-local `ImplicitLinkLaw`, the stand-in for `FixedCooling` and not a
  material law, gives the module's fixed conductance and the flux's
  derivatives in both cells' enthalpies.

Until the module retired, three steps of each fixture were held to it, for
both storage files:
- steady motoring and regeneration (1/64 s, four substeps, the load's work
  into a package): bit for bit, held to 1e-14 relative and 1e-11 K;
- the hot dynamic fixture (cells of 0.01 J/K, a diode threshold rising
  20 mV/K, 1/1024 s):
  - constant capacitances: within 4.2e-15 relative and 2.2e-11 K;
  - charge curves: within 9.3e-15 relative and 4.9e-11 K;

  held to 1e-13 (2e-13 for charge curves) and 1e-9 K. The output high
  package passes 300.75 K, so its diode threshold feeds back;
- the same fixture as the module continues it, with only the output high
  switch on in the second step and the first gates again in the third. The
  graph schedules the drive's program with a `SignalSource`:
  - constant capacitances: within 2.8e-15 relative and 1.4e-11 K, held to
    1e-13 and 5e-10 K;
  - charge curves: within 5.5e-15 relative and 2.7e-11 K, held to 1e-13 and
    1e-9 K.

  Keeping the first gates ends the second step elsewhere;
- the module's heat-assembly fixture, a link of 0.25 W/K between the
  inductor's cell and the motor's (`FixedCooling` on the module's heat
  network), in steady motoring, then the hot fixture with the same link:
  - steady: the drive bit for bit and the cells within 1.7e-13 K, held to
    1e-14 relative and 1e-11 K;
  - hot, constant capacitances: bit for bit; charge curves: within
    7.8e-16 relative and 1.6e-11 K; held to 1e-13 and 1e-9 K.

  The graph's first steady step also meets the module's closed form for
  backward Euler through the link, X = hG(Q0 - Q6)/(1 + 2hG): within
  5.7e-14 K, held to 1e-12 K. The drive's laws at the two linked cells do
  not vary with temperature here, so the link reaches the cells only.

Each run is also checked without the module; these checks remain:
- the steady runs against `steady_drive`'s closed form over the three
  periods: each cell k times its steady heat, the drive back at its start
  but for the battery's charge, k dt `cell_current_a` lower. Within
  1.1e-15 relative and 2.7e-14 K, held to 2e-14 and 1e-12 K;
- the other runs against the owner's step built by hand from the graph's
  old state at its end temperatures, and each cell's backward-Euler
  balance (with the link's flux where linked). The owner is the card's
  authored laws in closed form: constant resistances, a diode threshold of
  0.5 V + 0.02 V/K (T - 300 K). The owner agrees bit for bit and the
  balances within 5.5e-16 J (4.4e-14 J with the steady link), held to
  1e-14 relative and 1e-14 J (1e-12 J with the link).

With the module comparison switched off, the lead's heat in the
polarization's cell and the timeline returning mean temperatures still
failed these checks. The references themselves failed with the owner at
the start temperatures, with half the threshold slope, and with the
battery's charge kept in the steady closed form.

The graph refuses two steps, as the module did:
- a storage claim qualified only up to 300.000001 K (the graph with `Law`
  naming `heat_input_high`);
- a step that is not one program.

After each refusal the graph's drive and cells keep their states. The
kernel's ledger closes within 1e-9 J.

Nine defects planted in the graph, one at a time, each failed at least one of
these tests:
- the lead's heat in the polarization's cell;
- the stator given the rotor's torque;
- a switch law that does not ask its storage claim;
- the timeline returning mean temperatures;
- the storage claims in reverse order;
- the load's work left outside;
- the schedule's third step on the second program;
- the link's conductance doubled;
- the timeline without the link's law.

`battery_motor/nine_state.rs` declares the lane `battery_motor` still
couples by hand (from here on "the module"), switches without output
storage, the same way: the drive in the `NoDeviceStorage` lane, its switches
asking no storage claim. The harness takes the lane as a type parameter
(`Lane`); each lane's owner is built by hand (the gate lane's in
gate/graph.rs, at its cells' temperatures). Against the module, three
steps of each fixture:
- steady motoring and regeneration (1/64 s, four substeps, the load's work
  into a package): bit for bit, held to 1e-14 relative and 1e-11 K;
- the feedback fixture (cells of 0.001 J/K, a diode threshold rising
  20 mV/K, 1/1024 s, one substep): within 4.4e-14 relative and 1.2e-9 K
  (1.2e-12 J in these cells), held to 5e-13 and 1e-8 K. The output high
  package passes 300.1 K in the first step.

Without the module, the steady runs meet the closed form within 1.1e-15
relative and 2.7e-14 K, held to 2e-14 and 1e-12 K. The feedback run meets
the owner, `PreparedCommutatedBatteryMotorCycle` built by hand, bit for bit,
and each cell's balance within 5.3e-17 J, held to 1e-14 relative and
1e-14 J.

Both sides refuse:
- a card qualified only up to 300.0001 K, cells and devices alike, which
  the step passes. The module refuses the cells' state
  (`MaterialDomain`). In the graph the cells' solve rejects the step, the
  timeline asks for half of it, and the drive refuses a step that is not
  one period; drive and cells keep their states. Qualified to 350 K, the
  graph takes the step;
- a channel limit of 0.01 A: the owner's refusal on both sides;
- a step that is not one program.

Planted with the module comparison switched off, one at a time, each failed:
the lead's heat in the polarization's cell; the owner at half the threshold
slope; the steady closed form keeping the battery's charge; the owner by
hand with one more substep.

The module's finite liquid coolant on its heat assembly (the fixture's
claims as a liquid, 2 kg at 299.5 K) is an eighth timeline cell with no
port and its own caloric law, joined to the inductor's cell by the
stand-in at the heat assembly's conductance. Over three steps of 1/64 s:
- against the module, the drive bit for bit and the cells within 1.1e-13 K,
  held to 1e-14 relative and 1e-11 K;
- against the owner by hand, the drive bit for bit and each cell's balance
  within 8.1e-14 J, the coolant's by its own caloric curve; held to
  1e-14 relative and 1e-12 J;
- the first step meets the module's backward-Euler closed form exactly,
  held to 1e-12 K, and the coolant warms.

Planted with the module comparison switched off, each failed: the link's
conductance doubled in the stand-in; the coolant's temperature read from
the inductor's cell; the coolant linked to the motor's cell; the owner's
coolant starting at 300 K. One plant cannot fire: swapping the coolant's
caloric law for the fixture's swaps a law identical by construction, since
the module's liquid card copies the fixture's claims. What holds the
coolant instead is the closed form from 299.5 K with a cell of 2 J/K.

Not carried from the module's nine-state tests: its clock budget and final
cancellation (a step out of circuit trials, or cancelled at its last poll,
publishes nothing), its single-claim binding, and its drive bound at one
state refusing a step from another.

Not claimed: torques and rail voltages that change within the program;
reverse recovery and avalanche; storage, OCV or bearing drag that change
with temperature. The
retired lanes' refusals of a source that changed or was bound another way
retired with them; the graphs do not cover them.

## Multirate stepping (bead `frankensim-rc-1109-root-3gc12.5.9`)

`System::plan_multirate` assigns components to rate groups from each component's
`limit_dt_s` and optional per-node `node_max_dt_s` / `node_substeps` bounds
(refusing overrides coarser than a component's own `limit_dt_s` or exceeding
`max_substeps`). Any load-line loop stays intact in one rate group (promoted to
its fastest member when `promote_load_line_blocks` is `true`, otherwise refused
with `MultirateError::SplitLoadLineBlock`).

`System::step_multirate` advances each rate group on its own local clock
`dt_g = H / substeps_g` over one macro window `[t_0, t_0 + H]`:

- **Bitwise single-rate identity**: when every rate group has `substeps == 1`,
  `step_multirate` and `step_multirate_to_events` delegate directly to
  `System::step` and `System::step_to_events`, reproducing state, clock,
  outputs, and ledger accounts bit for bit.
- **Conserved cross-rate accumulators**: fast-to-slow port outputs, port
  energies, entropies, stream masses, and carried inventories are summed across
  substeps in double-double (`fs_math::dd::Dd`). Cyclic cross-rate stages solve
  the macro window-average interface equations with dense safeguarded Newton
  (generalizing `fs_phs::multirate::MultiratePair` to arbitrary `Component`
  graphs and connection types).
- **Declared slow-to-fast interpolation**: `SlowInterpolation::ZeroOrderHold`
  holds slower outputs constant across fast substeps so cross-rate work, heat,
  impulse, mass, and carried quantities cancel to double-double roundoff
  (`interpolation_numerical_loss_j == 0.0`, certifying `Color::Verified`).
  `SlowInterpolation::Linear` ramps slower outputs linearly across fast
  substeps (after step 0) and computes the exact quadrature energy discrepancy
  in double-double, booking it explicitly into `interpolation_numerical_loss_j`
  and the receiving node's `numerical_loss_j` account (leaving
  `Ledger::energy_residual_j` closed while coloring the window
  `Color::Estimated`).
- **Synchronized clocks and atomic rollback**: component `next_event_s`,
  mid-substep shortening requests, scheduled boundary events, and
  `System::step_multirate_to_events` guards synchronize all rate groups to one
  shared `Sum` clock. Every refusal or cancellation restores the exact entry
  checkpoint, clock, and ledger without partial publication.
- **Limits**: held/interpolated cross-rate interfaces do not eliminate macro-
  window phase error on strongly coupled oscillators; refine the macro window
  alongside local clocks.

