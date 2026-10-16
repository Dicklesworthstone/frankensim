# CONTRACT: fs-contract

Assume-guarantee component contracts (plan addendum, Proposal E): certified
design motifs whose `(envelope ⇒ guarantee, certificate)` compose into system
claims via envelope containment.

## Purpose and layer

Layer L3. Depends on `fs-iface` (L3, `SpaceType` for typed interface
quantities) and `fs-evidence` (UTIL, the `Color` lattice). Pure, deterministic
composition logic — it does not run solves or produce certificates, it
composes existing ones. The `deployment` module also owns the typed contract
between a floating design twin and one exact deployed target; it freezes proof
assumptions and adjudicates proof evidence but does not implement a controller,
plant solver, compiler, or theorem prover.

## Public types and semantics

- `Interval { lo, hi }` — a closed, finite, ordered interval;
  `Interval::new` rejects non-finite/inverted; `contains(other)` is inclusive
  on the boundary.
- `Envelope` — an interval box over named, `SpaceType`-typed quantities
  (`with(quantity, space, interval)`); `OperatingConditions` — the system
  model's computed interval per quantity.
- `Contract { name, interface: SpaceType, linear, envelope, guarantee,
  certificate: Color, requires }` — an assume-guarantee motif. `color_ok()`:
  a NONLINEAR contract may not carry a verified-color certificate.
- `ContractLibrary` — the certified-motif catalog (`insert`, `get`).
- `compose(&lib, root, &ops) -> Result<SystemClaim, ContractError>` — resolves
  `root` + its transitive `requires`, checks each member's operating conditions
  land inside its envelope and its color discipline holds, and returns a
  `SystemClaim` whose certificate is the WEAKEST member's color.
- `ContractError` — `BadInterval` / `UnknownContract` / `MissingCondition` /
  `OutsideEnvelope` / `ColorDiscipline` / `CircularDependency`.
- `deployment::DeploymentRefinementSpec` names both transition systems, the
  exact device and toolchain, state/observation maps, units/frames, clocks,
  quantization/saturation, plant, environment/disturbances, fault and safe
  state, horizon/invariant/objective, permitted error, relation strength,
  assumptions, capabilities, and bounded offline proof resources.
- `DeploymentRefinementProblem::admit` validates that complete seam and seals
  it. v1 refuses missing maps, schema drift, implicit unit/frame conversion,
  non-commensurate clocks, contradictory numeric/timing envelopes, omitted
  faults, different safe states, unavailable capabilities, unbounded proof
  resources, and malformed identities or bounded sets.
- `RefinementRelation` keeps trace inclusion, approximate simulation,
  approximate bisimulation, robust invariant, and performance bound distinct.
  `ProofAxis` independently names numeric, temporal, functional, and safety
  obligations; evidence on one axis cannot discharge another.
- `DeploymentRefinementManifest` retains schema-v1 canonical bytes and their
  FNV-1a replay root. A live target retains exact identities, maps,
  numeric/timing/fault/safe-state semantics, relation, budgets, and assumptions;
  environment and disturbance sets may narrow but never enlarge.
- `discharge_universal_claim` is the only constructor for a
  `DeploymentRefinementReceipt`. It requires exactly one established static or
  exhaustive proof for every axis against the frozen root. Measurement,
  Unknown, Refuted, missing/duplicate axes, zero hashes, and stale roots are
  typed refusals.

## Invariants

- SOUNDNESS: the composed certificate is never tighter than the weakest
  member's (its `ColorRank` is the minimum over members) — the Gauntlet
  contract-composition property.
- Envelope containment is inclusive on the boundary; a quantity with no
  operating condition, or one outside its envelope, blocks composition (fail
  closed — a guarantee is only asserted where provably inside).
- A nonlinear contract cannot be verified-color.
- The requires-graph is acyclic; a shared sub-contract (diamond) is resolved
  once, not flagged as a cycle.
- Deployment evidence is target- and assumption-specific. Changing a target or
  toolchain version, saturation, relation strength, environment domain, fault
  model, or safe state invalidates reuse. Measured traces cannot manufacture a
  universal-refinement receipt.
- Numeric, temporal, functional, and safety proof axes remain distinct. A
  receipt carries all four records in stable order and the exact relation it
  established.

## Error model

Structured `ContractError` and `DeploymentRefinementError` values (refusals
that teach), never panics for admitted input.

## Determinism class

Fully deterministic on one toolchain/ISA: composition is a pure function of
`(library, root, ops)`; members are returned sorted. Deployment manifests use
ordered maps/sets and exact IEEE-754 bits, so identical admitted inputs produce
identical canonical bytes and roots. Cross-toolchain proof equivalence is not
claimed; the toolchain is a semantic identity field.

## Cancellation behavior

Composition is synchronous pure work. Deployment admission/identity work is
bounded to 64 entries per collection and 65,536 aggregate UTF-8 bytes. The
offline proof budget carries nonzero work, memory, wall-time, and cancellation
poll-stride limits; external proof engines must honor them. This crate does not
run those engines or poll cancellation itself.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/contract.rs` (Proposal E, 10 cases): interval validation + inclusive
boundary containment; composition when conditions land inside; the
weakest-member soundness invariant (verified + estimated → estimated);
outside-envelope and missing-condition rejection; nonlinear-cannot-be-verified
color discipline; circular-dependency rejection; the diamond shared
sub-contract; unknown-contract rejection; determinism.

`tests/deployment.rs` (I05.6a; G0/G3): four-axis obligation separation;
missing maps; observation-schema, clock, unit, and frame mismatch; fault
omission; safe-state and capability refusal; saturation, target-version,
environment, fault, and safe-state drift; sound environment narrowing;
schema/root/byte replay; stable serialization; missing/duplicate proof axes;
measured/Unknown/Refuted refusal; and non-replay across relation strengths.

## No-claim boundaries

- Composition rule v1 is deliberately primitive ENVELOPE CONTAINMENT over
  interval boxes — it does not model general assume-guarantee refinement,
  nonlinear superposition, or probabilistic envelopes.
- Contracts CARRY certificates (`Color`); this crate does not produce or verify
  the certificates themselves (the solvers + fs-evidence do). The composed
  color is the weakest member's, by rank.
- Operating conditions are supplied by the caller (the system model); this
  crate does not compute them.
- The `SpaceType` interface tag is carried on quantities; cross-contract
  coupling-type compatibility checking (via fs-iface's checker) is a later
  integration.
- Deployment v1 supports identity-preserving units and frames only. An explicit
  certified conversion/transform contract is future work; names that merely
  look convertible are refused.
- Integer-commensurate clocks are an admission prerequisite, not proof of WCET,
  schedulability, latency, or closed-loop timing. Those remain temporal-axis
  obligations for an independent checker.
- A `DeploymentRefinementReceipt` says only that supplied evidence discharged
  the frozen obligations. It does not validate the physical plant abstraction,
  certify hardware, qualify a compiler, prove regulatory compliance, or extend
  beyond the exact target, environment, fault set, horizon, objective, and
  relation in its manifest.
