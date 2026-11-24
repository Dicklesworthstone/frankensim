# fs-ports

## Purpose and layer

Layer: L3. The domain-neutral owner of the port vocabulary and of the
executable Component contract that coupled physics is composed through.

The L3 domain-neutral owner of the existing PortSchema v2, StreamPort v1,
closed-window accounting, relation admission, and scalar coupling primitives.
The implementation is extracted unchanged from fs-couple; that crate reexports
these exact types and remains the owner of assembled physical systems.

fs-iface, fs-qty and fs-exec are the direct runtime dependencies; fs-exec (L0)
is there only for the `Cx` that `component` passes through. Geometry and
physics operators may use the common types without importing acoustic, fluid or
thermal solvers. The schema, stream and accounting-window types introduce no
new quantities, wire versions, source authority, material laws, cancellation
behavior or whole-system passivity claims.

## Public types and semantics

- Port description: `PortSchema` (v2), `PortKind`, `PortValueShape`,
  `PowerPairing`, `CoordinateBinding`, `PortTimestamp`, `ConservationRole`,
  `StableId`. A schema is admitted at construction and is descriptive.
- Streams: `StreamPort` (v1) with its energy charts and
  `surface::SurfaceTransport`.
- Accounting windows: `AccountingWindowInterval`, `WindowBalance`,
  `WindowAuditReport`.
- Lossless relations: `ConservativeJunction`, `JunctionSchema`, `JunctionLaw`,
  and the adoption trait `PortedModule`.
- Iteration weight: `BoundedAitkenRelaxation`.
- Executable contract: `component::{Component, PortTangent, LoadLine,
  StageResidual, PortSpec, DualityMetric, PortLayout, Class, Receipt, Proposal,
  MassStepAudit, EntropyStepAudit, audit_step, audit_mass_step,
  audit_entropy_step}` and `component::conformance`.

The sections after the ten required ones give the semantics of each group.

## Invariants

- A `PortSchema` exists only if effort times flow, with the declared measure,
  has the dimensions of power and matches its kind.
- A lossless relation cancels interface power structurally: signs are exactly
  plus or minus one and a Dirac matrix is exactly skew. No tolerance is applied.
- Closed-window residuals are accumulated in one compensated, fixed-order
  reduction.
- In the Component contract, port energy is computed from the two values on the
  wire, and `audit_step` is the only first-law check of a step.

## Error model

Refusals are structured values. `CoupleError` covers schema, stream, window and
relation refusals. `component::ContractViolation` names the law a step broke. A
component's own `Error` type passes through the contract unchanged. Non-test
code contains no `unwrap`, `expect` or `panic!`.

## Determinism class

Deterministic. Non-test code is sequential, uses no threads, clocks, random
numbers or hash-ordered containers, and reduces in a fixed order.

## Cancellation behavior

Schema, stream, window and relation calls are bounded and do not poll. The
Component contract passes a `Cx` to `propose` and to the tangent and residual
calls; the component polls it and returns its own cancellation error.

## Unsafe boundary

None. The crate denies `unsafe_code`.

## Feature flags

None.

## Conformance tests

- `adoption::tests`, `surface/tests.rs` and `relaxation::tests`.
- `tests/component_reference.rs` for the Component contract and
  `tests/load_line_reference.rs` for its load-line tier, with the three
  `compile_fail` doctests in `component`.
- `tests/field_ports.rs` for field duality and independently computed traction
  work on an inertial surface.
- `tests/mass_ports.rs` for optional inventory closure, stream signs, an
  energy-balanced mass leak, cancellation and unrepresentable input refusals.
- `tests/entropy_ports.rs` for the optional entropy account: port entropy
  and booked generation against measured inventories, the strict tolerance
  and how an offset reference loosens it, refusals, and the conformance
  kit's audit and unit rescaling of a declaring body.
- `tests/heat_ports.rs` for the sign and size of a heat port's entropy and
  its refusal at a temperature that is not above zero; and for the field heat
  port: with one point it equals the scalar port, with three its energy and
  entropy are the sums over its points to the bit (400 drawn cases in each
  role), and one bad temperature or a vector of another width gives NaN.
- The fs-couple and fs-feec batteries named under "Purpose and layer".

## No-claim boundaries

No whole-system passivity, stability or convergence claim. No claim about the
validity of any constitutive law, equation of state or measured value that
flows through a port. Schema admission grants no energy closure. Each section
below states its own limits.

## Component contract v0

`component` is the executable composition contract (bead
`frankensim-rc-1109-root-3gc12.5.1`). It is version 0 and may change until the
five kernel slices of that epic pass on it. Its first implementer outside this
crate is `fs_mbd::component::RigidBody`.

**Tiers.** `Component` is mandatory: ports, class, `stored_energy_j`,
`limit_dt_s`, `next_event_s`, `propose`, `admit`. `PortTangent` adds the JVP and
VJP of the port map. `StageResidual` adds a residual in the component's own
unknowns. `LoadLine` is for a relation at power ports that is not a function of
the flow (dry friction holds any force up to its limit at zero slip). The
optional tiers require `Component`; three `compile_fail` doctests pin that.

**The load line.** A component of class `NonSmooth` cannot be iterated on the
flow it receives, because its effort is not a function of it. What its
neighbours do with an effort is a function, and frozen at an iterate it is
affine: `flow = released - compliance * effort`. `LoadLine::line_ports` names
the power ports of role `FlowIn` the line spans (one or more, in increasing
order: a pad on one body, a link between two). `propose_on_line(old, dt,
inputs, compliance, new, outputs, cx)` is `propose` with the released flows in
`inputs` at those ports and the line's square matrix `compliance`, row-major
over the line's coefficients. The component returns the efforts its relation
allows on that line; its candidate, its other outputs and its receipt are for
the flows on the line, not the released ones. The matrix need not be symmetric;
its symmetric part is not negative; a coefficient whose flow is prescribed has
a zero row and column. With a zero matrix the call is `propose`, to the bit.
The solver that supplies the line is the system kernel's; this crate fixes the
call and checks an implementation of it.

**The end of a step.** A port that returns a flow returns its MEAN over the
step: that is what a held effort does its work on. What the flow is at the END
of the step the mean does not say, and a relation that depends on it has to be
told (dry friction holds a pair that ends its step with no slip).
`Component::end_of_step(port)` is what a component promises about it, per power
port that returns a flow:
- `EndOfStep::Midpoint`, the default: the flow ends the step at twice the mean
  less the flow it began the step with. Exact for a component stepped by the
  midpoint rule. A component whose step is not that rule's must not keep the
  default.
- `EndOfStep::Told(q)`: the flow at the step's end is returned on signal
  output `q`, coefficient for coefficient, by the same proposal.
- `EndOfStep::Untold`: neither.

`LoadLine::line_ends` names, per port of the line, the signal input that
receives that port's flow at the step's end, or none (the default: none). A
port with none takes its flow to end the step by the midpoint relation. A port
with an end input is handed, after the line's square matrix, one more row of
`compliance` per coefficient of the input: `end[k] = released_end[k] -
sum(compliance[n + k][j] * effort[j])`, with the released end flows in
`inputs` there. A component that names none is handed the square matrix alone,
as before. Who makes the wire, and which pairings are refused, is the
kernel's (`crates/fs-system/CONTRACT.md`). `check_load_line` takes the matrix
with those rows and refuses another length.

**State.** The system owns every state. `propose(old, dt, inputs, new, outputs,
cx)` reads the accepted state and writes a candidate slot whose content on entry
is unspecified. It is pure: equal arguments give bitwise-equal results and
nothing else changes. `admit` is the last refusal, after every partner proposed.
A component may answer `Proposal::Shorten { dt_s }` instead of stepping.

**Snapshot codec.** A component may declare a codec of its state (bead
`frankensim-rc-1109-root-3gc12.5.4`):
- `Component::STATE_CHARTER`, an fs-exec `StateIdentityCharterV2`, `None` by default.
- `encode_state(state, words)`, which appends the state as plain 64-bit words (an `f64` as its bits). It returns `false` by default, for no codec.
- `decode_state(words)`, which returns the state those exact words hold, or `None` to refuse them.
- `identity_words(words)`, which appends the bits of the parameters that fix the component's physics. It appends nothing by default.

All four are defaulted, so a component that declares none is unchanged. A system is written to a durable snapshot only when every member has a codec. It refuses a snapshot whose members have other charter identities or other identity words (`crates/fs-system/CONTRACT.md`). `cargo run -p xtask -- check-state-charters` scans every charter literal for collisions.

**Port values.** One representative pair per port and step. The side the
component receives is held over the step; the side it returns is the step
average of the conjugate. Flow and heat rate are positive into the owner.
`Carrier::Power(kind)` books `dt * sum(effort * flow)`; `Carrier::Heat` carries
temperature [K] and heat rate [W] and books `dt * flow`, with entropy
`dt * flow / temperature`; `Carrier::Signal` books nothing. The kernel computes
port energy from the wire values (`PortSpec::energy_in_j`), so both ends of a
lossless connection book the same number.

`PortSpec::heat_field(id, role, width)` declares a heat port over a trace of
`width` points: one temperature [K] and one heat load [W] into the owner per
point. The loads are already integrated against the trace's basis, so the port
books `dt * sum Q_i` with no metric, and the entropy `dt * sum Q_i / T_i`, each
load at its own point's temperature (the nodal rule). With one point it is
`PortSpec::heat`, the same value, and the scalar entropy expression is its own
code path. The port carries no identity of its trace: two field heat ports
joined point by point must list their points in the same order, which only
their owners know.

`PortSpec::field_power` declares a field power port whose effort and flow
vectors have the same width and whose energy is `dt * e^T M f`.
`DualityMetric::new(width, entries)` takes row-major, finite, exactly symmetric
coefficients (signed zeros compare equal), requiring positive Cholesky pivots
and finite factors. It retains the factor; `multiply` and `solve` use caller
buffers without allocation. Shape, coefficient, asymmetric pair, pivot and
factorization refusals are structured; vector actions refuse wrong lengths or
nonfinite arithmetic. Ordinary power-port arithmetic is unchanged.
The metric does not infer a physical basis, spatial accuracy, units or
conditioning, and Cholesky admission is not an interval certificate of positive
definiteness. Rectangular and indefinite pairings are outside this API.
Connections must preserve the declared metric; matching carrier and width
alone is insufficient. Six G0/G3 tests cover off-diagonal work in both roles,
identity pairing, independent multiplication and solve, admission and arithmetic
refusals, and a field component whose Euclidean energy account is rejected.

`Carrier::Stream` (`PortSpec::stream(id, role, carried)`, width `2 + carried`)
carries material between a storage and a transport element. The storage
(`FlowIn`) returns what it holds is like: pressure [Pa], specific enthalpy
[J/kg], then one value per unit mass for each carried quantity. The transport
element (`EffortIn`) returns what crosses the port, positive into itself: mass
rate [kg/s], enthalpy rate [W], then the rate of each carried quantity. The
port books `dt * enthalpy rate` as energy and `PortSpec::mass_in_kg` gives
`dt * mass rate`; neither is a product of the two sides. A transport element
multiplies its mass rate by the values of the side the mass leaves (the donor
rule), which is its law and not the kernel's. Entropy carried by a stream is
not on the wire. `tests/component_reference.rs` checks widths, both roles'
energy and mass, the refusal of a signal role and of slices of the wrong
length.

`PortSpec::identified_stream` additionally binds a `StreamBasis`: an immutable
medium/source-law identity, an enthalpy-reference identity, and an ordered
`StreamQuantity` list. Species masses use the existing `SpeciesId`; other
conserved quantities have a stable semantic identity and SI inventory dimensions
(the dimensions of the amount, not its rate). Duplicate identities refuse even
when declared with different units. Empty lists still retain the medium and
energy datum. `PortLayout` rejects identified ports of one owner that disagree
on a shared inventory ordinal. It does not reorder or convert coordinates.

Connections must compare the entire basis. An unidentified `stream` is not a
wildcard for an identified port; two unidentified ports retain their original
width-only contract and its caller obligations. Width overflow is a structured
refusal. Basis admission does not authenticate source declarations, implement
reactions, or change wire arithmetic. `tests/carried_ports.rs` checks ordering,
dimensions, empty bases, unchanged signed transfers and declaration refusals.

**Mass.** `Component::stored_mass_kg` optionally measures mass from a state.
Its default is `None` (unavailable); zero storage must be declared as
`Some(0.0)`. `audit_mass_step` checks finite, nonnegative available inventories
against signed stream transfers. Step, vector and individual `dt * mass rate`
validation applies even when both inventories are unavailable. Availability
at only one endpoint is refused. The audit reduces scaled terms in a fixed,
compensated order with tolerance `128 * EPSILON * max(old, new,
sum(abs(port transfers)))`, avoiding overflow in balanced large transfers.
`MassStepAudit` reports net port input, inventory change, balance defect and
tolerance. Unavailable inventories return `None`, never a conservation claim.
The existing conformance check measures both endpoints and runs this audit;
energy receipts do not supply its tolerance. This adds no species, element,
momentum or charge audit to the Component contract.

**Entropy.** `Component::stored_entropy_j_per_k` optionally measures entropy
[J/K] from a state. Its default is `None` (unavailable), so an owner that
does not declare it is not audited and changes nothing; zero storage must be
declared as `Some(0.0)`. Entropy has a free reference: an owner that declares
it states its reference, and available values must be finite but may be
negative. `audit_entropy_step` checks that the entropy the ports carry in
(`dt q / T` on a heat port) plus the receipt's `entropy_generated_j_per_k`
equals the measured change. Step, vector, generation and individual transfer
validation applies even when both inventories are unavailable. Availability
at only one endpoint is refused, and so is an inventory change that is not
representable. The tolerance is the mass account's, strict and with no claim
from the owner: `128 * EPSILON * max(|old|, |new|, sum(abs(port transfers)) +
generated)`, reduced in the same fixed, compensated order. It is relative to
the inventories the owner reports, so a reference far from the states
loosens it in proportion. A declaring owner books what its step generates,
so its entropy closes to rounding; a time-discretization defect belongs in
`numerical_loss_j`, not in this account. `EntropyStepAudit` reports net port
input, generation, inventory change, balance defect and tolerance.
`check_component` measures both endpoints and runs this audit after the mass
audit; `check_unit_rescaling` compares the start's entropy in the new units.

**First law.** `audit_step` is the one place a step's energy is tested:
`sum(port energy in) = stored_change_j + numerical_loss_j` within the receipt's
`energy_tolerance_j`; the receipt's stored change against re-measured endpoint
energies; dissipation and entropy production not negative; residual within the
owner's tolerance. Both energy checks allow `f64::MIN_POSITIVE` J on top of
the receipt's tolerance: the f64 underflow allowance. Below the normal range
an f64 rounds by an absolute step, one subnormal (5e-324 J), which no
tolerance relative to the step's energies covers; a receipt claims only that
relative roundoff. At or above 2^-968 J (about 4.0e-292 J) the sum
absorbs the allowance, so no decision there changes. G0 tests in
`tests/component_reference.rs`: a balance and an endpoint one subnormal off
are admitted at a zero tolerance, and a defect of 1e-300 J at a zero
tolerance is still refused, as `EnergyBalance` and as `EndpointEnergy`. Physical dissipation stays in the component as stored
thermal energy or leaves through a heat port; it is never booked as numerical
loss. The audit trusts `energy_tolerance_j`. The owner must make it cover the
roundoff of `stored_energy_j` at both ends of the step, which scales with the
terms of the stored energy when they cancel (kinetic against potential), and
a tolerance claimed too wide lets a leak through the audit. A port list that
does not fit the layout is refused, not indexed.

**Errors.** `Component::Error` is the owner's type and reaches the caller
unchanged. Contract breaches are `ContractViolation`.

**Determinism and cancellation.** The contract code is sequential. `audit_step`,
`audit_mass_step` and `audit_entropy_step` do not allocate; `PortLayout::new` allocates its two offset vectors.
Components poll `Cx::checkpoint` and return their own cancellation error.

**Conformance.** `component::conformance::check_component` and
`check_port_tangent` are the checks every implementation runs.
`check_component` also asks `propose` once under a context built from the
caller's on a gate whose cancellation is requested before the call, and
requires an error: the law "propose polls cx: a step asked under a cancelled
context is refused". G0 tests in
`tests/component_reference.rs`: a spring-mass-damper with a mechanical and a
heat port passes; 20,000 driven steps close energy to 1e-11 relative; a heat
leak, a read of the candidate slot, a count kept in a `Cell` through `&self`
that moves an output by units in its last place, a wrong VJP sign, and a
`propose` that never polls `cx` each fail the law they break; a cache kept in
a `Cell` under the bits of the accepted state, the step and the force, which
`propose` never shows, passes and observes what the reference observes;
cancellation returns the owner's error; an audit given a port list longer or
shorter than its layout is refused. `check_exact_linear` holds a component
that declares `Class::ExactLinear` to that claim between two operating points
of one step. A solver may keep, at one step, what it derived from the slopes
at an earlier step of the same length, and a run that kept it must take the
steps of a run that derived it again; so along every unit input direction the
JVP at the two points must be the same bits, a zero's sign included. It
refuses a component of another class and two points at different steps.
`check_port_tangent` asks it of every component of the class, at its case and
at the case's stale state with each input `x` moved to `1.5 x + 0.25`. G0
tests in `tests/component_reference.rs`: the reference declares `ExactLinear`
only without its damper, whose heat is quadratic in the mean velocity, and
there has the same slopes at two states and two forces; a damped variant that
declares the class is refused by the slope of its heat with respect to the
force when only its state differs, when only the force differs, and by the
tangent tier; an undamped variant whose zero heat slope carries the sign of
the mean velocity is refused by a zero and a minus zero. The law compares two
points and is no proof of the claim at others, and it asks nothing of a
component with no declared tangent. `check_load_line` checks the
load-line tier at one point: the line's ports, the size of the matrix, purity
(including the stale candidate slot), the step audit at the flows on the line,
`admit`, and that a zero matrix gives `propose` to the bit. G0 tests in
`tests/load_line_reference.rs`: a dry pad on one port conforms over 81 released
flows, held (zero flow, zero loss, the force the line asks) and sliding (at its
limit, the way it was released); a dry link conforms on a line over two ports
with a matrix that is not symmetric; a loss accounted at the released flow
fails the energy balance, a mandatory tier that differs from the line tier at
a zero matrix fails that law, a line on a heat port and a matrix of the wrong
size are refused. The kit cannot check that the efforts satisfy the
component's own relation; only the component can state it.
`check_time_reversal` checks, at one point, a component that a test claims
is reversible. The test supplies three things: the reversal `R`, a view of
the state (`observe`) and the bounds. The laws are:
- `R` is an involution;
- every port is a power or heat port without a field metric;
- the step from `x0` passes `check_component` and dissipates nothing, and so
  does the step from `R(x1)` with the odd inputs negated;
- the backward outputs are the forward ones times their parity;
- `R(x2)` returns to `x0` as `observe` sees it, and in stored energy.

Parity comes from the port roles, not from the test: effort and temperature
are even, flow and heat rate odd. The bounds are derived in the test from the
rounding of the component's own arithmetic; the kit holds no tolerance of its
own for this check. G0 tests in `tests/component_reference.rs`:
- the reference without its damper reverses;
- with its damper it is refused at the step that dissipates;
- a reversal that keeps the velocity fails the outputs law;
- a reversal that is not an involution is refused first.

`check_unit_rescaling` checks, at one point, that a component keeps its law
in other units. The test supplies `build`, which makes the component in the
units that base-2 exponents name (`UnitExponents`, one per SI base dimension
`[m, kg, s, K, A, mol]`), each parameter times the factor of its dimensions
(`unit_factor`), and `scale_state` for the state. The kit scales the step,
and every input and output by the factor of its port's dimensions, taken
from the port alone: a power port's from its kind, a heat port's kelvin and
watts. Every factor is a power of two, so the laws hold to the bit, with no
tolerance:
- every exponent is at most `UNIT_EXPONENT_LIMIT` (8) in size, and every
  number scaled or compared stays zero or normal;
- every port is a power port without a field metric, or a heat port;
- the rescaled component has the base ports, and both cases pass
  `check_component`;
- the rescaled step is taken or shortened as the base step is, and its
  outputs, its end state and its receipt (stored change, dissipation,
  numerical loss, energy tolerance, entropy generated) are the base ones in
  the new units;
- so are the start's stored energy, mass and entropy, `limit_dt_s` and
  `next_event_s`.

G0 tests in `tests/component_reference.rs`:
- the reference passes in units changed in each base dimension alone and in
  all six together;
- with its damping left in the old units it is refused at the outputs, and
  it passes a change of the metre alone, which damping does not see;
- with a signal port, or with an exponent of 9, it is refused before a step.

Power pairing, effort times flow in watts, belongs to the seven port kinds,
not to a component: `PortSchema::try_new` refuses a pair whose product is
not W. fs-couple/tests/couple.rs pins every kind's canonical pair through it
(`port_schema_v2_scalar_seed_migration_goldens`,
`extended_port_kinds_have_watt_pairings_and_kind_identity`).

**No claim.** v0 has no entropy on stream ports and no multistage port
values. The kernel cannot check that a component's `identity_words` name every parameter of its physics; only the component knows them. The kit
does not yet check `StageResidual` or temporal order.
The audit's underflow allowance is absolute, in joules, so below 2^-968 J
(about 4.0e-292 J) the audit is not covariant under a change of units.
`check_unit_rescaling` is exact only while every non-polynomial function is
applied to a dimensionless group and nothing overflows or goes subnormal. A
function of a dimensional value (`exp(T / 300)`, `ln(p)`) breaks it just as
a hidden constant does, which is a finding too. A step that the audit
admits only within its underflow allowance can be refused in other units,
and `check_unit_rescaling` then fails at the rescaled `check_component`.
Which dimensions each parameter and each part of the state has, in `build`
and `scale_state`, stays the test's claim. Stream, signal and field-metric
ports are refused; the receipt's `residual` and `tolerance` are in the
owner's units and are not compared.
`check_time_reversal` sees a state only through the outputs, the stored
energy and what the test's `observe` returns, so what none of them shows is
not checked. Whether `R` is the physical reversal, and whether `observe` and
the bounds are complete, stays the test's claim. No `Component` declares
itself reversible. The kit proposes twice with equal arguments, so it sees a change made
through interior mutability only where that change alters what a call
returns; it cannot prove that a cache is keyed by the accepted state, which
stays the component's claim. The cancellation law sees only a context
cancelled before the call: that `propose` polls "at bounded intervals"
inside a long step is not checked, nor that the error it returns is the
owner's cancellation error rather than another refusal. A passing
conformance check says a component keeps its own accounts; it says nothing
about the accuracy of its physics.

The existing fs-couple/tests/couple.rs battery continues to exercise this owner
through its public reexports. fs-couple/tests/port_identity.rs pins type and
error identity across both entry points. FEEC's terminal-relative and integral
batteries consume the direct fs-ports dependency.

Closed-window first-law, mass, element, charge and entropy residuals accumulate
the two stored inventories and every signed boundary term in one compensated,
fixed-order reduction. Entropy mechanisms remain separate through that final
reduction and through exchange grouping. A weak heat or diffusive contribution
must not disappear merely because two separately rounded large totals agree.
The existing coupling battery covers one-unit deficits beside 1e16-unit
inventories/streams and cancellation within and across entropy rows. This is a
floating-point conservation check, not an interval bound or a recovery of
precision already lost by the producer of an individual inventory.

## Port-schema adoption contract

`PortedModule` binds a stable port order to the pre-integrated transfers of the
last accepted step, positive into the module. Rejection preserves the preceding
record; before the first accepted step it is `None`. Logical ticks never replace
physical-time integration. Schema admission alone grants no energy closure.

`JunctionSchema` describes 2..=4097 conjugate endpoints using common effort,
common flow, or an exactly skew finite Dirac matrix. Sign pullbacks apply to
both power factors. Pairwise metadata admission is the existing
`ConservativeJunction::new` contract, including orientation, clocks and roles.
The ordinary two-port common-effort instance retains its exact acceptance and
errors. Mixed domains require explicit power-preserving coordinate charts;
dimension relabelling is not a chart. No new solver or wire version is added.

G0 unit tests in `adoption::tests` cover N=2/3/64, cap and cap+1, exact two-port
admission parity, and non-skew/sign/identity refusal. They establish structural
power cancellation, not whole-system passivity or component audit closure.

Coordinate charts accompany canonical schemas. A mass-normalized chart retains
the raw N/sqrt(kg), m sqrt(kg)/s units and an explicit reference-mass conversion;
it does not claim a physical point force or measured modal mass. The signed
volume projection, when present, is retained separately. `ModuleWindowAudit`
distinguishes a real report, no accepted step, and a named unavailable physical
account. Observation errors are visible through `last_window_transfer_error`.

## Moving control-surface streams

`StreamPort::try_from_uniform_surface` evaluates the original source state on
an oriented uniform face or quadrature point. With material velocity u and cut
velocity w, mass uses `rho*A*(u-w).n`, momentum uses
`mdot*u+A*(p*n-tau*n)`, and the existing enthalpy chart includes
`mdot*(h+|u|²/2+gz)+p*A*w.n-A*(tau*n).u`. The swept pressure term is
separate from deviatoric work and is present at zero-mass moving walls.
Constituents and entropy use the same signed relative mass flow. Reusing the
retained chart refuses replacement of its geometry-derived mass or momentum
rates.

`SurfaceTransport` exposes relative/swept volume flow, advective and stress
momentum, material-side stress power, work on the moving cut, and angular flux
about an explicit origin. Open-stream convection is not silently interpreted
as a body force; slipping-wall work is not inferred heat. Source EOS/stress
validity, heat and species diffusion, spatial quadrature, mesh/GCL advancement
and finite-body evolution remain with their original owners.

G0/G3 native tests cover impermeable-wall impulse/work, Galilean boosts,
uniform-flow geometric conservation, off-origin angular flux and rigid-motion
power, orientation reversal, stationary-chart parity and invalid/rebound
inputs. These are local transport identities, not a whole-system stability
or physical-validation claim.

## Bounded vector relaxation

`BoundedAitkenRelaxation` chooses a convex iteration weight from a fixed-order
vector of finite applied/solved coordinates. Finite bounds satisfy
`0 < minimum <= initial <= maximum <= 1`. Each coordinate uses its running
maximum magnitude; both previous and current residuals are represented in that
same metric. Mixed physical units therefore do not enter raw dot products.
A common residual scale prevents squared tiny defects from underflowing.

Invalid controls, nonfinite coordinates and a changed coordinate count refuse;
a failed update preserves the full relaxation history. Zero coordinates and an
empty vector are admitted. Stagnant or unrepresentable secant ratios retain the
previous bounded weight. Candidate-local work and storage are O(n); the physical
caller owns its work/cancellation budget, source domains, physical residual
checks and publication. This is acceleration, not a stability, convergence,
uniqueness or passivity certificate. Original scalar relaxation is unchanged.

Five G0/G1/G3 unit regressions cover the exact oscillating affine map, a passive
resistor's independent operating point, independent unit/polarity changes,
zero/stagnant/extreme coordinates, and atomic invalid-input refusal.

## Independently measured carried inventories

`Component::stored_carried_quantity(state, carried)` optionally measures the
inventory of a zero-based stream coordinate after mass and enthalpy. `None`
means unknown, and zero storage is explicit `Some(0.0)`. Finite signed values
are admitted because a coordinate can represent momentum as well as species
mass. Declaring an inventory claims it is conserved under the declared stream
transfers; reactions and internal production are not silently inferred.

`PortSpec::carried_quantities` names the count and `carried_quantity_in` reads
the integrated signed transfer. `audit_carried_step` independently closes each
coordinate using scaled compensated sums and its own 128-epsilon rounding
allowance. Both endpoint measurements must be available or both unavailable.
Nonfinite transfers still refuse when measurements are unavailable. The
conformance kit performs the same audit. Calls allocate nothing.

The owner declares corresponding coordinate order and units across its stream
ports; an absent coordinate contributes zero. This does not match chemical
identities, balance reactions, or certify a mechanical momentum law. Three
G0/G3 regressions in `tests/carried_ports.rs` cover both port roles, signed and
small independent coordinates, missing/nonfinite inventories, invalid wires,
unrepresentable transfers, and compensated cancellation at extreme scales.
