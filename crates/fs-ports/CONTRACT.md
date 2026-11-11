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
- Executable contract: `component::{Component, PortTangent, StageResidual,
  PortSpec, PortLayout, Class, Receipt, Proposal, audit_step}` and
  `component::conformance`.

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
- `tests/component_reference.rs` for the Component contract, with the two
  `compile_fail` doctests in `component`.
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
unknowns. Both optional tiers require `Component`; two `compile_fail` doctests
pin that.

**State.** The system owns every state. `propose(old, dt, inputs, new, outputs,
cx)` reads the accepted state and writes a candidate slot whose content on entry
is unspecified. It is pure: equal arguments give bitwise-equal results and
nothing else changes. `admit` is the last refusal, after every partner proposed.
A component may answer `Proposal::Shorten { dt_s }` instead of stepping.

**Port values.** One representative pair per port and step. The side the
component receives is held over the step; the side it returns is the step
average of the conjugate. Flow and heat rate are positive into the owner.
`Carrier::Power(kind)` books `dt * sum(effort * flow)`; `Carrier::Heat` carries
temperature [K] and heat rate [W] and books `dt * flow`, with entropy
`dt * flow / temperature`; `Carrier::Signal` books nothing. The kernel computes
port energy from the wire values (`PortSpec::energy_in_j`), so both ends of a
lossless connection book the same number.

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

**First law.** `audit_step` is the one place a step's energy is tested:
`sum(port energy in) = stored_change_j + numerical_loss_j` within the receipt's
`energy_tolerance_j`; the receipt's stored change against re-measured endpoint
energies; dissipation and entropy production not negative; residual within the
owner's tolerance. Physical dissipation stays in the component as stored
thermal energy or leaves through a heat port; it is never booked as numerical
loss. The audit trusts `energy_tolerance_j`. The owner must make it cover the
roundoff of `stored_energy_j` at both ends of the step, which scales with the
terms of the stored energy when they cancel (kinetic against potential), and
a tolerance claimed too wide lets a leak through the audit. A port list that
does not fit the layout is refused, not indexed.

**Errors.** `Component::Error` is the owner's type and reaches the caller
unchanged. Contract breaches are `ContractViolation`.

**Determinism and cancellation.** The contract code is sequential. `audit_step`
does not allocate; `PortLayout::new` allocates its two offset vectors.
Components poll `Cx::checkpoint` and return their own cancellation error.

**Conformance.** `component::conformance::check_component` and
`check_port_tangent` are the checks every implementation runs. G0 tests in
`tests/component_reference.rs`: a spring-mass-damper with a mechanical and a
heat port passes; 20,000 driven steps close energy to 1e-11 relative; a heat
leak, a read of the candidate slot and a wrong VJP sign each fail the law they
break; cancellation returns the owner's error; an audit given a port list
longer or shorter than its layout is refused.

**No claim.** v0 has no entropy on stream ports, no field ports with a
duality matrix, no multistage port values and no snapshot identity. The kit
does not yet check `StageResidual`, unit rescaling, time reversal or temporal
order. A passing conformance check says a component keeps its own accounts; it
says nothing about the accuracy of its physics.

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
