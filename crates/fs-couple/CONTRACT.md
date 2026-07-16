# CONTRACT: fs-couple

Multiphysics composition through a versioned port-thermodynamic vocabulary: a
lossless interface relation, explicit storage/dissipation/source primitives,
and caller-supplied scalar balance instrumentation without false passivity.

## Purpose and layer

Layer L3 (multiphysics coupling). Dependency-light: `fs-qty` supplies the
canonical six-base `Dims`/power dimension and `fs-iface` supplies neutral field
function-space roles. No domain solver is a dependency.

## Public types and semantics

- `PORT_SCHEMA_VERSION = 2`; `StableId` admits a canonical transport-safe ID.
- `PortValueShape` distinguishes scalar, non-empty vector/tensor, and field
  values (`fs-iface::SpaceType`); `PowerPairing` distinguishes scalar product,
  Euclidean dot product, and field duality with explicit integration-measure
  dimensions.
- `CoordinateBinding` makes basis, frame, and positive orientation explicit;
  `PortTimestamp` carries a named logical clock and tick;
  `ConservationRole` is canonicalized into a sorted, duplicate-free list.
- `PortSchema::try_new` binds stable ID, legacy/current `PortKind`, six-base
  effort/flow dimensions, value/field shape, coordinates, power pairing,
  timestamp, and conservation roles. PR-1 admission checks that the contraction
  matches the shape, effort × flow has watt dimensions without exponent
  overflow, and energy is named as a conservation role.
- `PortKind::scalar_seed_schema` migrates the existing mechanical, fluid, and
  thermal scalar seeds without inventing identity, frame/basis, or clock data.
- `ConservativeJunction::iterate_added_mass_fixed` /
  `iterate_added_mass_aitken` are schema-bound migration bridges for the
  retained mechanical scalar fixture; they must remain bitwise equal to the
  legacy results and are not general FSI operators.
- `ConservativeJunction` admits two conjugate schemas and evaluates the finite
  scalar seed as `SchemaInterconnection`/`SchemaPort` (shared effort, opposite
  flow). It localizes metadata mismatches and refuses non-scalar/non-finite
  evaluation.
- Four closed primitive variants are distinct in `PortPrimitive`:
  `ConservativeJunction`; `StorageElement` with Hamiltonian/free-energy state
  schema and constitutive-gradient operator; `DissipativeRelation` with a
  constitutive family plus mandatory monotonicity/nonnegative-production
  evidence reference; and `SourceOrReservoir` with an explicit signed
  `AccountingBoundary` carrying the same basis/frame/orientation as its port
  plus included-source/external-reservoir treatment.
- The legacy scalar `PortKind`/`Port`/`interconnect` surface remains as the
  migration oracle: `Port { effort, flow, kind }` has `power` = effort × flow
  and `conjugate_to` checks the same scalar kind.
- `interconnect(kind_a, kind_b, effort, flow) -> Result<Interconnection,
  CoupleError>` — a Dirac structure (shared effort, opposite flow) whose
  `interface_power` is `0` exactly (power-conserving by construction); refuses
  incompatible ports. `interface_power(&[Port])` = `Σ effort·flow`.
- `EnergyAudit` — `record`, `max_generation`, `is_passive(tol)`: the legacy
  `is_passive` name checks only caller-supplied scalar interface imbalance at
  each recorded exchange. A nonzero balance is a bug alarm, not a proof of
  whole-system or closed-window passivity.
- `AitkenRelaxation::new(omega_init, omega_max)` + `next_omega(residual)` — the
  scalar Δ² dynamic relaxation factor, magnitude-capped.
- `iterate_fixed_relaxation` / `iterate_aitken` — the added-mass interface
  fixed point under fixed vs Aitken relaxation → `FsiResult { converged, steps,
  solution, final_residual }`.

## Invariants

- The Dirac interconnection conserves interface power EXACTLY (to roundoff) —
  the G0 law; incompatible ports are refused.
- Every admitted v2 schema has a shape-compatible power contraction, a checked
  watt-dimensional effort/flow product, and an explicit energy role. Field
  duality includes its integration-measure dimensions in that checked product.
- Stable relation IDs cannot alias their owned port/boundary IDs. A
  conservative pair requires distinct port IDs, matching physical metadata and
  clock/timestamp, plus outward-from-owner conventions on both sides. PR-1
  refuses common-frame orientations until a public pullback can prove them.
- Storage, dissipation, and source ownership remain separate typed claims. A
  conservative junction cannot stand in for any of them; a dissipative
  relation cannot exist without an evidence reference; a source cannot exist
  without a named accounting boundary whose coordinate/sign binding exactly
  matches the port.
- The energy audit reports an interface-balance failure exactly when some
  caller-supplied exchange has absolute imbalance above `tol` or is non-finite.
- On the added-mass fixture (`μ ≥ 1`): naive staggering (`ω = 1`) diverges while
  Aitken-relaxed coupling converges to `x* = c/(1+μ)`; Aitken never takes more
  steps than a stable fixed under-relaxation.

## Error model

Structured `CoupleError`; no panics. In addition to legacy incompatible kinds,
errors cover invalid IDs/empty shapes, pointwise or field-measure dimension
overflow, non-power products, shape/pairing mismatch, missing energy role,
localized schema conjugacy mismatch, identity aliasing, accounting-boundary
coordinate mismatch, scalar/non-scalar misuse, wrong added-mass port kind, and
non-finite schema-bound values.

## Determinism class

Fully deterministic: schema admission canonicalizes roles by enum order;
interconnection, relation construction, audit, and iterations are pure
functions of their inputs.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/couple.rs` (16 cases): v2 scalar-seed migration goldens for all three
legacy kinds and bitwise whole-result migration of the added-mass fixture;
schema fail-closed metadata; localized junction mismatch; non-scalar refusal
by the scalar evaluator; field-duality measure dimensions; all four primitive
descriptors and identity-alias refusals;
legacy power-conjugate ports; exact interface power conservation and
incompatible-port refusal; energy-audit imbalance and non-finite alarms; the
Aitken Δ² factor; naive staggering diverges where Aitken stays stable; Aitken
accelerates over stable fixed relaxation; light added mass converges naively;
determinism.

## No-claim boundaries

- The FSI fixture is the classic LINEARIZED added-mass interface map
  (`H(x) = −μx + c`) — enough to reproduce the instability and its fix; a full
  nonlinear FSI solve over real fluid/structure subsystems is the consumer.
- `AitkenRelaxation` is the scalar Δ² relaxer; the vector INTERFACE
  QUASI-NEWTON (IQN-ILS) accelerator and MULTIRATE co-simulation are staged.
- PR-1 proves generic effort × flow power dimensions and preserves only the
  three scalar seed kinds. PR-2 owns the rotational/electrical/magnetic/
  chemical vocabulary and kind-specific dimensional admission.
- The scalar evaluator does not execute vector/tensor/field or general
  multi-port Stokes–Dirac operators. `StorageElement` and
  `DissipativeRelation` carry durable public operator/evidence references; this
  crate does not execute or validate the referenced domain law in PR-1.
- Field-duality admission checks dimension arithmetic including the declared
  measure; it does not prove quadrature, trace pullback, orientation, or the
  numerical duality operator.
- `AlongFrame`/`AgainstFrame` can describe schemas and audit boundaries, but
  PR-1's conservative junction refuses them rather than assuming an unproved
  effort/flow pullback. A later neutral transfer API must make that transform
  explicit before admission.
- PR-3 owns `StreamPort`, energy-accounting charts, and exact crosswalk refusal.
  PR-4 owns closed-window first-law, mass/element/charge, entropy-generation,
  and optional exergy audits. An `AccountingBoundary` descriptor is not one of
  those audits.
- `fs-iface::SpaceType` records a neutral field role; it does not supply a
  mortar/Nitsche/harmonic transfer operator or prove inf-sup compatibility for
  a particular domain adapter.
- The energy audit's balances are supplied by the caller each exchange; wiring
  them onto the ledger is the coupling driver's integration.
- Dirac interface losslessness does not establish passivity of component
  storage/dissipation/source laws, spatial or temporal discretization,
  interface transfer, nonlinear iteration, multirate windows, or the coupled
  system. Those obligations require a signed, closed-window energy audit.
