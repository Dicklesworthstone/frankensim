# CONTRACT: fs-ledger

> Status: ACTIVE (Design Ledger, schema v2). Owns the core schema + Rev S
> extension tables, BLAKE3 content addressing, the WAL/snapshot concurrency
> contract, and — since schema v2 — forkable worlds, `at(t)` views,
> `explain()`, the replay audit, and unreferenced-artifact GC (`travel`
> module).

## Purpose and layer

The Design Ledger (plan §11.2, Bet 10): FrankenSQLite-backed system of record
for content-addressed artifacts, event-sourced ops with the frozen Five
Explicits, lineage edges, metric time series, the autotuner cache, and the
fine-grained event stream. Layer: L6 (HELM). Runtime deps: `std` + `fsqlite`.

## Public types and semantics

- `Ledger` — one connection + the pragma contract (WAL, `synchronous=FULL`,
  `busy_timeout`, enforced foreign keys) + versioned migrations
  (`PRAGMA user_version`; idempotent DDL batches in `schema::MIGRATIONS`).
- `ContentHash`, `Blake3`, `hash_bytes` — in-house BLAKE3 (plain hash mode,
  32-byte output), pure safe Rust; artifact identity everywhere. The
  implementation is OWNED by the UTIL crate `fs-blake3` (bead 7uq9) and
  re-exported here unchanged — same paths, same bits.
- Artifacts: `put_artifact` (≤ `STORAGE_CHUNK_LEN` inline; larger stored as
  `artifact_chunks` rows because fsqlite has no incremental-blob API),
  `ArtifactWriter` (streaming; hashes incrementally, stages chunks under a
  provisional key inside a writer-owned transaction, promotes on `finish`),
  `get_artifact` / `read_artifact_chunks` / `artifact_info`,
  `verify_artifact_integrity` (full re-hash), `corrupt_artifact_for_test`.
- Ops/lineage: `begin_op` (validates the Five Explicits field-by-field;
  units travel inside the typed IR, the other four are mandatory columns),
  `finish_op` (exactly once; `ok|error|cancelled`), `op`, `link` (FK-checked
  `in|out` edges).
- Streams: `record_metric` (finite REAL only), `append_event` /
  `append_events` (batched, atomic), `tune_put`/`tune_get` (upsert keyed
  kernel × shape-class × machine fingerprint).
- Rev S extension tables (sparse v0, uniform `(name UNIQUE, body JSON)`
  shape): `put_extension`/`get_extension` over `requirements`, `model_cards`,
  `evidence`, `scenarios`, `constraints`, `capability_probes`, `imports`,
  `unsafe_capsules`.
- Hygiene: `lint()` (orphan edges/metrics/chunks, storage-shape and length
  invariants, half-finished ops, dangling branch references) — all-zero on
  any healthy or crash-recovered ledger.
- Time travel (`travel` module, schema v2): `fork`/`branches`/`branch_diff`
  (a fork is a new op-log branch sharing every artifact by hash; visibility
  = own ops + ancestors' up to each fork point), `begin_op_on` (branch +
  recorded `ExecMode`), `at_time` (consistent views at arbitrary instants:
  outcomes not yet written are masked, unfinished ops' outputs invisible),
  `explain` (full causal trees, depth-limited, DAG-deduped, loud on orphan
  inputs), `replay_verdict` (deterministic ops must reproduce output hashes
  exactly; fast divergences reported without failing),
  `gc_unreferenced_artifacts` (edge-less artifacts only; referenced
  artifacts are immortal).

Schema divergences from plan Appendix D, both deliberate: `JSON` columns are
STRICT-legal `TEXT` with `json_valid()` CHECKs (Appendix D as written is not
valid STRICT SQL), and `artifacts` gains `len`/`chunk_count` +
`artifact_chunks` for bounded-memory large-field storage.

- `tombstone` module (addendum Proposal E, bead lmp4.13): the TOMBSTONE
  LEDGER — swarm memory's cheap half. `Descriptor` (name + dimensioned
  params) computes a π-space signature via fs-regime's exact Buckingham
  machinery (the PRIMARY, domain-native index: dimensionally-equivalent
  deaths collide across raw parameters) and a deterministic hashed
  feature-vector embedding (tokens + magnitude decades; Franken-only, no
  external model). `TombstoneIndex`: automatic appends on falsification
  kills (carrying the Proposal-6 falsifier JSON) and on abandoned
  branches ABOVE a cost threshold; `pre_exploration_check` (the
  orchestrator gate — π-space first, embedding second);
  `fund_with_distinguisher` VALIDATES the cited feature (must name a
  real parameter differing by ≥ 0.05 decades — free text refused) and
  logs accepted distinguishers on the tombstone so they accumulate;
  `re_exploration_rate` is the proposal's kill-criterion metric;
  `flush_to_ledger` persists rows as `tombstone` events.

- `vcs` module (addendum Proposal 10 base verbs, bead lmp4.9): VERSION
  CONTROL FOR PHYSICS — commits/branches/checkout over Merkle roots,
  free-riding on `travel`'s forkable worlds. A COMMIT is the binary
  Merkle root of a branch's visible frozen ops (leaf = BLAKE3 of the
  Five Explicits + outcome + diag + sorted linked-artifact hashes;
  wall times and rowids EXCLUDED, so logically identical histories
  produce identical roots across ledgers and runs); commits chain to
  their branch predecessor and persist as `vcs-commit` events.
  CHECKOUT materializes a committed view; `checkout_delta` returns the
  symmetric-difference op frontier (the `perturb()`-style delta a
  recompute solver consumes — nearby checkouts cost |delta|, not
  |history|). `merge_views` splits base/only-A/only-B for the
  diff/bisect/merge consumers; `storage_audit` measures the
  "N branches ≈ 1× + deltas" sharing claim; `op_artifact_hashes` and
  `commit_leaf` are the public leaf surface. GC safety is inherited:
  `gc_unreferenced_artifacts` walks lineage reachability, and the VCS
  suite proves no live-branch artifact is ever collected.

## Waiver trust boundary (bead qmao.1.1)

Annotation and authorization are DISTINCT: `Waiver` (id/signer/reason
strings) is a human memo that travels in provenance but authorizes
nothing — `derive` refuses any upgrade claim regardless of annotation.
The only path past a laundering refusal is `derive_waived` with a
`WaiverGrant`: a versioned, length-prefixed canonical payload bound to
the node name, exact parent provenance hashes (replay to another node
or lineage fails), the claimed color, the color-upgrade scope, a
signer key id, and an expiry day — verified through the caller-
supplied `WaiverVerifier` capability before any write. A grant carries and
signs the full `Color::canonical_bytes` payload (signing encoding v2), not only
the color rank name, so authorization for one interval, validity regime, or
estimator payload cannot authorize another. The in-tree
default is `NoWaiverVerifier` (refuses everything): no cryptographic
capability ships in this crate, so promotion is impossible until a
Franken-compliant signature verifier is wired in (the no-crypto
no-claim). Node provenance hashes use a versioned length-prefixed encoding
(v3), including the bit-exact canonical color bytes — the former v2 encoding
used display-rounded color JSON and could alias distinct floating-point
payloads. Refusals are structured
(`WaiverRejection`: scope/node/color/lineage mismatch, expiry, bad
signature) and grants re-verify from the stored ledger node.

## Invariants

1. Artifact identity = BLAKE3 of content; identical bytes dedupe to one row
   (concurrent duplicate insert resolves to dedupe, never an error).
2. Storage shape: inline XOR chunked; `len` always equals stored byte count;
   chunk `seq` is dense from 0. Enforced by CHECKs and re-checked by `lint`.
3. Ops are event-sourced facts: `(t_end IS NULL) = (outcome IS NULL)` is a
   table CHECK; an op finishes at most once (`DoubleFinish` otherwise).
4. Edges only reference existing ops and artifacts (enforced FKs).
5. A crash-recovered ledger lints clean: transactions make op+edges+metric
   groups all-or-nothing (kill -9 battery, `ledger_007`).
6. Wall-clock timestamps are provenance envelope, never content identity.

## Error model

All fallible APIs return `LedgerError` — structured variants with stable
`code()` strings and actionable Display text: `Open`, `FutureSchema` (newer
file refused, never clobbered), `Sql`, `Busy` (retryable contention —
busy/locked/write-conflict; retry with backoff), `MissingExplicit` (names the
offending Five Explicits field), `Invalid` (names the field),
`Corrupt`, `NotFound`, `DoubleFinish`, `WriterInTransaction`. Never panics
across the crate boundary.

## Determinism class

Content hashing is bit-stable across runs, thread counts, and ISAs (pure
function). Row ids, timestamps, and physical file bytes are NOT deterministic
and are excluded from identity. Deterministic replays should pass logical
times to `begin_op`/`append_event` (caller-controlled `t`).

## Cancellation behavior

No compute kernels; all calls are short transactions. A dropped
`ArtifactWriter` rolls back its transaction leaving zero residue (tested).
Once fs-exec lands, ledger writes stay on the latency lane per plan §5.2;
scope-tree integration is the fs-obs sink bead's scope.

## Unsafe boundary

None. Safe Rust only (workspace `deny(unsafe_code)`); the BLAKE3
implementation is pure safe Rust.

## Feature flags

None. All v0 behavior is `[S]` default-path.

## Conformance tests

`tests/conformance.rs`: official-vector BLAKE3 battery (0 B → 2 MiB+1,
covering multi-level trees), seeded streaming-split property, versioned
migration + future-version refusal, dual-path chunked dedupe + round trip,
corruption-fails-loudly (inline + chunked), concurrent snapshot readers
during a write sweep (monotone + internally consistent), kill -9 crash
battery (6 seeded rounds → lint-clean + integrity-clean), and an events/sec
throughput smoke ledgered as a metric. `tests/travel.rs`: genuine-v1 →
v2 migration with history intact, fork storage audit (N forks = 1× artifacts
+ deltas) + branch independence, replay audit battery (clean /
deterministic-failure / fast-divergence), explain() full-lineage
reconstruction with loud orphan-input failure, at(t) monotone mid-sweep
consistency, and a kill -9 battery during fork traffic. Unit tests in
`src/lib.rs`, `src/hash.rs`, and `src/travel.rs` cover the API surface and
edge cases.

## Speculation telemetry (bead lmp4.3, schema v3)

Schema v3 adds the `speculation` extension table (uniform Rev S
shape) carrying the four solve-node fields `(proposer_id, accepted,
bound, iterations_saved)` keyed by solve-op identity — ADDITIVE: the
migration regression test proves every pre-existing table still
answers queries. The economics control loop lives in fs-verify
(HELM-side); this ledger stores telemetry, it does not drive solves.

## Three-color write gate (bead qmao.1)

`colors::ColorGraph` is the WRITE-TIME gatekeeper over fs-evidence's
color schema: derived nodes' colors are COMPUTED from their parents
(with regime re-checks against the current execution state,
auto-demoting validated parents whose regime the state has exited,
demotion events logged); a claimed color that outranks the derivation
REFUSES with the capping parents named (the laundering refusal, G3
gauntlet-tested); the only override is a signed `Waiver` that appears
in the ledger row AND participates in the node's provenance hash —
it cannot be dropped without changing history. Rows are canonical
JSON lines for the event stream. Note: this module adds fs-evidence
as a runtime dependency (the colors are its types).

## No-claim boundaries

- Multi-process multi-writer access: unclaimed (FrankenSQLite documents this
  as partial; use one process, one connection per thread).
- BLAKE3 keyed hashing, key derivation, XOF output beyond 32 bytes: not
  implemented.
- Cryptographic security claims: the implementation matches official vectors
  but has no side-channel or performance hardening (scalar, unoptimized).
- Throughput numbers are smoke floors, not roofline claims (§14 discipline:
  real claims need machine fingerprints and acceptance bands).
- Branch DELETION and cross-branch merge (as opposed to merge-view
  queries): not provided; branches are append-only history.
- `at_time` trusts caller-supplied op timestamps; callers that write
  non-monotone `t_start`/`t_end` get views consistent with what they
  recorded, not with wall-clock truth.
- Multi-GiB single artifacts: chunk storage bounds row sizes, but the
  streaming path is verified at the tens-of-MiB scale only so far; fsqlite
  transaction memory behavior at multi-GiB scale is unmeasured.

## No-claim boundaries (tombstones)

- Retrieval is exact-scan over in-memory indexes (linear); ANN/sublinear
  retrieval and the FrankenTorch encoder upgrade land when volume
  demands them — the deterministic feature vector is the documented
  degradation path (polish note honored in reverse).
- π-signature comparison requires the SAME group structure; explorations
  with different physics never collide (and are never suppressed).
- The orchestrator PROTOCOL (querying before funding) is enforced by the
  agent-orchestration layer; this module provides the gate, the
  validation, and the metric.
- Descriptor parameters must be positive (π-space is multiplicative);
  signed features belong in the embedding text, not the signature.

## No-claim boundaries (vcs)

- The `Vcs` registry is in-session (commits also persist as events);
  cross-session registry reconstruction from event rows lands with the
  diff/bisect beads that need it.
- `checkout_delta` names the ops to reconcile; executing the delta-solve
  is fs-recompute's contract (lmp4.7/8), not this module's.
- Merge ADJUDICATION (the sheaf Hodge split) is lmp4.12's crown jewel;
  this module supplies its base/only-A/only-B views.
- Commit leaves fold linked-artifact HASHES (content-addressed);
  artifact bytes are shared by the store, not re-hashed per commit.
