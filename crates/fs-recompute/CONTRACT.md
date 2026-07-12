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
  skips), `content_hash()` (versioned length-prefixed binary canonical
  serialization: params sorted by key, floats by BITS, inputs in order,
  fs-ledger's Blake3-class tree hash), and `to_row()` (valid JSON carrying
  all seven fields + slack, including bit forms). Caller strings cannot
  inject field delimiters or collide with structured fields.
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
  the recompute reason, malformed tolerances return `InvalidTolerance`,
  and `Miss` is honest absence.
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
   zero-slack hit, deficits are named, malformed requested tolerances
   fail closed, and skip identity ignores recorded tolerances (rcs-003).
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
7. Error budgets and slack burns are finite, non-negative magnitudes;
   malformed values refuse without mutation (rcs-007).

## Tolerance-aware invalidation (bead lmp4.7, feature-gated)

`invalidate::plan` computes the recompute frontier for a perturbation:
deltas flow `Σ L_e · δ(parent)` through EVERY frontier node (skipped
nodes are STALE by their bound — staleness reaches consumers scaled by
their sensitivities), each node absorbs against its OWN effective
slack, and the recompute set is closed UPWARD along delta-carrying
edges (`PulledByDescendant`: fresh bytes need fresh inputs). Skip
verdicts carry VERIFIED-color interval claims in their rows;
`apply_plan` BURNS absorbed bounds into runtime state (`burned`,
SEPARATE from the immutable record identity — the suite caught an
early design where burning mutated the hashed record and broke
identity), so repeat perturbations see the spent slack. Fail-closed
hardening: exact ties recompute; non-finite sensitivities force
recompute; finite negative sensitivities also force recompute because
they are not magnitude bounds; duplicate perturbations at one source
sum by the triangle inequality; negative slack never skips; δ = 0 is an
empty frontier.
Skip YIELD is the R4 health metric; loose bounds degrade gracefully
to hash-memoization behavior, still correct.

Invariants: flow-through absorption + upward closure (inv-001);
fail-closed zoo (inv-002); the G3 SOUNDNESS battery — over seeded
traces on an executable DAG, EVERY node's final value (cached or
fresh) lies within its tolerance of full-recompute truth, and the
falsifier's forced recomputes agree within their certified bounds;
any violation is Sev-0 (inv-003); graceful degradation with yield
measured (inv-004); verified-color claims + slack burning (inv-005).

## perturb() API + cache policy (bead lmp4.8, same feature gate)

`api::RecomputeApi` is the operator-facing surface: `perturb(node, δ)`
returns a FIRST-CLASS `PerturbPlan` — the minimal frontier, its
estimated cost from MEASURED per-node costs (Proposal 8's planner
input), the hash-memoization baseline cost, and the verified-color
certificates for everything skipped — pure until `commit` (which burns
slack and updates telemetry). Cache policy: `ensure_capacity` evicts
by COST-WEIGHTED score (recompute-cost × measured hit-probability,
lowest first, deterministic seq tie-break), pins untouchable, and a
pinned population exceeding the capacity is the STRUCTURED
`CacheFullOfPins` refusal — never an OOM. `SkipYield` is the per-op R4
dashboard with worst-first ordering (where bound-tightening effort
goes).

Invariants: diamond plans recompute exactly the un-absorbable
{source, tight} set with certificates for the rest, leaf/root
boundaries behave, plans are pure until commit (api-001); slack is
spendable through the API — repeated absorptions exhaust it (api-002);
cost-weighted eviction preserves hot expensive nodes that insertion-
order LRU would destroy, pins survive, saturation is structured
(api-003); per-op yields separate never-absorbing ops from absorbers,
dashboard live via fs-obs (api-004); the kill-criterion replay
machinery measures certified-vs-memo cost on a 100-variant trace
(fixture-scale; the production decision runs on recorded agent
traces) (api-005). Every plan keeps its verdicts and evidence rows
read-only and captures both the store mutation revision and a deterministic
fingerprint of certificate-relevant state. Committing another plan,
substituting a plan from different state, or otherwise changing the store
makes the plan a structured `StalePlan` refusal, validated before any burn
or telemetry mutation. Duplicate burns are aggregated before the atomic
fit check (api-006).

## Error model

`StoreError::DeterminismViolation` (stop-the-line, with likely-cause
teaching text: unordered reduction, unstable sort, uninitialized
padding), `UnknownNode`, malformed error/burn refusals, and `StalePlan`.
Nothing panics across the boundary.

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

`tolerance-invalidation` — the lmp4.7 invalidation algorithm AND the
lmp4.8 perturb()/cache-policy API, OFF by default per the Ambition-Tag
gating rule until the Gauntlet tier and kill-metric (≥2× median
wall-clock speedup vs plain memoization on recorded agent traces)
stay green. Adds fs-evidence (the verified-color skip claims).

## Conformance tests

`tests/conformance.rs`, cases rcs-001..rcs-007 — JSON-line verdicts,
seeded LCG randomness, the fs-obs slack-table event. Any
reimplementation must pass the suite unchanged.

## No-claim boundaries

- Cross-ISA certification: rcs-004 certifies across worker counts and
  completion orders on the host; the both-reference-ISA gate rides
  the perf/CI lane's remote runners (the fs-la golden-hash pattern).
- Invalidation traversal (dirty propagation through the DAG) and the
  cache-policy surface are the recompute-invalidate / recompute-api
  beads; this store supplies their pinning hooks.
- The SQLite-backed persistent form (fs-ledger schema v5 tables) is
  deferred; `snapshot()` v2 is the interim durable form.
- Slack SPENDING policies (which skips to take under a budget) are
  the recompute-api bead's.
- Sensitivity bounds are SUPPLIED (interval-derived by callers);
  adjoint-sharpened bounds (Proposal 1) tighten the loose ones.
- Path-sum accumulation is conservative for shared subpaths (no
  common-subexpression tightening yet).
