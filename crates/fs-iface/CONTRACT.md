# CONTRACT: fs-iface

Interface types + the coupling-graph static checker (plan addendum,
Proposal 13). Turns ill-posed couplings into compile-time errors before any
solve, using the FEEC periodic table as the interface type lattice.

## Purpose and layer

Layer L3 (FLUX-adjacent). Pure static analysis — NO numerical runtime
dependencies (std only). The FEEC periodic table (`H(grad)`/`H(curl)`/
`H(div)`/`L²`) is the type lattice; a static [`check`] validates a coupling
graph against an inf-sup-certified pairing registry.

## Public types and semantics

- `SpaceType` — the de Rham complex spaces in exact-sequence order
  (`HGrad`→0, `HCurl`→1, `HDiv`→2, `L2`→3 by `form_degree`); `H(grad)`
  full-continuity, `H(curl)` tangential, `H(div)` normal, `L²` none.
- `CouplingRole` — `Continuity` (matching the SAME field across an interface;
  legal iff both sides share a trace space) or `Saddle` (a mixed/saddle-point
  block; legality is inf-sup stability from the registry).
- `InterfaceField { id, space }`, `Coupling { id, trial, test, role }`,
  `CouplingGraph` (builder: `.field(id, space)`, `.couple(id, trial, test,
  role)`).
- `PairingRegistry` — a DECLARATIVE literature table of inf-sup (LBB) results.
  `standard()` seeds the well-established ones: certified `(H(div), L²)`
  (RT/BDM mixed Poisson/Darcy) and `(H(grad), L²)` (Taylor–Hood Stokes);
  known-unstable `(L², L²)` and `(H(grad), H(grad))` (equal-order LBB
  violations). `classify_saddle` returns `Certified{cite}` / `Unstable{reason}`
  / `Unknown`. `certified_partners` teaches a legal alternative.
- `check(&CouplingGraph, &PairingRegistry) -> CheckReport` — pure,
  deterministic. `CheckReport { admitted, findings }`; each `CheckFinding`
  carries the offending coupling id (localization), a `check` slug, `Severity`,
  a diagnosis, and a teaching `fix`.

## Invariants

- CONSERVATIVE BY DEFAULT: an `Unknown` saddle pairing is REJECTED
  (illegal-until-certified) — the checker never silently admits a pairing it
  does not recognize.
- Every finding is localized to a specific coupling id.
- An empty graph is vacuously legal (admits with no findings).

## Error model

Findings, not panics: malformed graphs (a coupling referencing an undeclared
field → `graph.field`; a field coupled to itself → `graph.self-loop`) and
illegal couplings (`coupling.continuity`, `coupling.infsup`) are `Reject`
findings with teaching fixes. `admitted` is false iff any `Reject` exists.

## Determinism class

Fully deterministic: couplings are processed in declaration order and findings
emitted in order; `check` is a pure function of `(graph, registry)`.

## Cancellation behavior

None — `check` is a bounded, synchronous pure function (no `Cx`, no I/O).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/checker.rs` (Proposal 13): the FEEC exact-sequence type lattice;
empty graph vacuously legal; same-space continuity legal; cross-space
continuity illegal (pressure-to-displacement); certified saddle pairings
(Darcy H(div)–L², Stokes H(grad)–L²) legal; known-unstable (P1–P1) rejected
with a certified-partner fix; UNKNOWN pairing rejected conservatively; missing
field and self-loop localized errors; multi-fault per-coupling localization;
determinism; empty-registry certifies nothing.

## No-claim boundaries

- The registry is a LITERATURE TABLE of inf-sup results. The checker
  guarantees coupling-graph LEGALITY against it; it does NOT verify that the
  corresponding saddle-point element families or solvers are implemented
  (those are fs-feec high-order / the solver-stack beads). Registry entries
  may name discretizations not yet built.
- Space typing is at the function-space ROLE level (de Rham periodic table);
  polynomial-degree/regularity refinement (e.g., Taylor–Hood degree pairing
  vs. equal-order) is captured coarsely — the registry marks the stable
  FAMILY, not per-degree admissibility.
- Symmetry harvesting (the second half of Proposal 13) is a separate bead.
- Wiring this checker into fs-ir admission / the FrankenScript coupling
  surface is a later integration.
