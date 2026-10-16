# fs-selfknow-e2e — CONTRACT

The STRUCTURE & SELF-KNOWLEDGE end-to-end battery (bead knh1.7): Layer
4 exercised as one script.

## Purpose and layer
Layer L4 (dev-facing e2e harness). Six stages over the six selfknow
crates: fs-iface, fs-symmetry, fs-spectral, fs-surrogate (ladder),
fs-adjoint (explain), fs-plan (voi).

## Public types and semantics
`VERSION` and `STAGES` (the ledger keys, in battery order). The battery
itself lives in `tests/suite.rs` behind `selfknow-e2e`.

## Invariants
Every stage's FAIL-SAFE assertion:
- unknown pairings rejected (illegal-until-certified);
- no false speedups (asymmetric inputs flagged; the perturbation bound
  contains the true asymmetry);
- no crying wolf (healthy gaps keep their color; degraded gaps demote
  merges);
- no infinite descent (leaking queries terminate at full order with a
  complete leak trail);
- no confabulated explanations (gate AND narrative refuse together);
- no false-zero VoI claims (a zero sampled flip fraction stays explicitly
  grid-qualified and non-authoritative).
Each stage appends REAL fs-ledger events (the audit trail).

## Error model
Test panics; the harness `expect`s ledger IO.

## Determinism class
Fully deterministic fixtures; ledger event times are logical.

## Cancellation behavior
Not applicable (a test battery; each stage is independent).

## Unsafe boundary
No `unsafe`.

## Feature flags
`selfknow-e2e` (default OFF) gates the battery for nightly lanes.

## Conformance tests
tests/suite.rs — stage-1 interface types (G0), stage-2 symmetry
falsifier (G3), stage-3 degraded-gap refusal (G4), stage-4 ladder
descent, stage-5 explanation reconciliation (G3), stage-6 VoI laws
(G0/G2).

## No-claim boundaries
- Both-ISA nightly scheduling is the CI lane's job; this crate ships
  the battery, not the cron.
- Stage fixtures are the member crates' canonical smalls — flagship-
  scale batteries live with the flagships.
