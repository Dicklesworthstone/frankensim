# CONTRACT: <crate-name>

> Status: SKELETON. This contract gains force as the crate's beads land; a crate
> must have a complete contract before it becomes a dependency target of other
> crates' real code (AGENTS.md).

## Purpose and layer
<one-line purpose>. Layer: <UTIL|L0..L6>.

## Public types and semantics
Skeleton only: `VERSION`. Populated by this crate's implementation beads.

## Invariants
None claimed yet. No-claim until stated.

## Error model
No fallible APIs yet. Errors will be structured values with machine-readable
diagnoses (Decalogue P10), never panics across the crate boundary.

## Determinism class
Target: Deterministic (bit-stable across runs and thread counts per ISA) unless
a section here documents a narrower class. Not yet verified — see G5.

## Cancellation behavior
All future hot paths poll cancellation at tile boundaries (Decalogue P7).
No compute paths exist yet.

## Unsafe boundary
None. `unsafe_code` is denied workspace-wide; any future capsule must be
registered per docs/CONVENTIONS.md and ship a SAFETY.md.

## Feature flags
None. Frontier features use `frontier-*`, moonshots `moonshot-*`, default off.

## Conformance tests
tests/conformance.rs (placeholder). Any reimplementation must pass this suite.

## No-claim boundaries
Everything: this is a skeleton. No numerical, performance, or safety claims.

## Golden couplings (bead y4pt)

Every golden/replay hash this crate pins must have a row in
`golden-couplings.json` declaring its upstream semantic surfaces; bumps
follow docs/GOLDEN_POLICY.md (committed tree, both build modes,
plausible root cause, registry updated in the same commit).
