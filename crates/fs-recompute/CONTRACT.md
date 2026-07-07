# fs-recompute — CONTRACT

## Purpose and layer

L6 (HELM). Proposal 2's STORE: a content-addressed Merkle DAG whose
nodes record `(op_id, input_hashes, params, code_version_hash,
rng_seed, achieved_error, required_tolerance)`, with the gap
`required_tolerance − achieved_error` as first-class SLACK — the
resource incremental recompute spends. The Error Ledger becomes a
build graph with a soundness certificate for every skip, and
DETERMINISM is promoted from implementation detail to CERTIFIED
CONTRACT (risk R2 owned here).

## Public types and semantics

- `NodeRecord` (the seven-field schema) with `slack()` (negative
  representable — over-budget nodes are first-class and never satisfy
  skips), `content_hash()` (canonical serialization: params sorted by
  key, floats by BITS, inputs in order, fs-ledger's Blake3-class tree
  hash), and `to_row()` (all seven fields + slack).
- `Store::put(record, artifact_bytes)`: content-addressed insert;
  identical record + identical artifact is a write-time memo hit
  (`Deduped`); identical record + DIFFERENT artifact bytes is
  `StoreError::DeterminismViolation` — the trip-wire that makes the
  determinism contract self-policing. STOP-THE-LINE, not a warning:
  tolerance-level memoization is unsound until the op is fixed.
- `Store::can_skip(record, new_tolerance)`: the skip-soundness oracle.
  Identity for skips excludes the recorded tolerances (a node cached
  under a looser requirement still hits if it ACHIEVED enough);
  `Hit{slack}` is the certificate, `ToleranceTightened{deficit}` names
  the recompute reason, `Miss` is honest absence.
- `Store::pin(node, PinReason::{EvidencePackage, Contract})`: pinned
  nodes are NEVER evicted; `evict_unpinned(keep)` removes oldest
  unpinned first (deterministic) and cannot touch pins by
  construction.
- `snapshot()`: canonical text form (fork/round-trip stability).

## Invariants

1. Node hashes are repeat-stable, param-order canonical, and sensitive
   to EVERY one of the seven fields (floats by bits); negative slack
   is first-class; 1000-deep chains are hash-stable; empty/single-node
   stores behave (rcs-001).
2. The determinism trip-wire: identical (record, artifact) dedupes;
   identical record with different bytes errors with both artifact
   hashes named (rcs-002).
3. Skip decisions carry slack certificates, the exact boundary is a
   zero-slack hit, deficits are named, and skip identity ignores
   recorded tolerances (rcs-003).
4. THE CERTIFICATION (G5-at-scale primitive): a fixture study —
   deterministic tile reduction (fs-exec `det_sum` per tile +
   order-fixed `pairwise_fold`) — produces BIT-IDENTICAL artifacts
   across {1,2,4,8} REAL worker threads and adversarial permuted
   completion orders; every re-put is accepted as a dedup by the
   contract (rcs-004).
5. Pins survive eviction; eviction is deterministic oldest-unpinned-
   first; pinning unknown nodes teaches (rcs-005).
6. Ledger rows carry all seven fields + slack; rows and snapshots are
   bitwise-deterministic across builds (rcs-006).

## Error model

`StoreError::DeterminismViolation` (stop-the-line, with likely-cause
teaching text: unordered reduction, unstable sort, uninitialized
padding) and `UnknownNode`. Nothing panics across the boundary.

## Determinism class

The crate's whole point. Store operations are BTree-ordered and
sequence-numbered; hashing is canonical; the conformance battery
certifies worker-count and completion-order independence of the
fixture study through the store's own trip-wire.

## Cancellation behavior

Store operations are O(log n) point operations (no long loops); the
fixture study's cancellation discipline belongs to fs-exec.

## Unsafe boundary

None. `#![forbid(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs`, cases rcs-001..rcs-006 — JSON-line verdicts,
seeded LCG randomness, the fs-obs slack-table event. Any
reimplementation must pass the suite unchanged.

## No-claim boundaries

- Cross-ISA certification: rcs-004 certifies across worker counts and
  completion orders on the host; the both-reference-ISA gate rides
  the perf/CI lane's remote runners (the fs-la golden-hash pattern).
- Invalidation traversal (dirty propagation through the DAG) and the
  cache-policy surface are the recompute-invalidate / recompute-api
  beads; this store supplies their pinning hooks.
- The SQLite-backed persistent form (fs-ledger schema v3 tables) is
  deferred; `snapshot()` is the interim durable form.
- Slack SPENDING policies (which skips to take under a budget) are
  the recompute-api bead's.
