# CONTRACT: fs-ledger

> Status: ACTIVE (Design Ledger, schema v3). Owns the core schema + Rev S
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
  `busy_timeout`, enforced foreign keys) + versioned migrations. Each
  `PRAGMA user_version` marker commits atomically with its DDL batch, and a
  fresh (v0) file initializes the ENTIRE ladder plus its final marker in one
  transaction — a crash mid-init leaves an empty v0 file, never a partial
  schema. SCHEMA ATTESTATION BEFORE ADVANCEMENT (bead gp3.18): a v0 file
  must contain NO user objects (pre-existing tables refuse initialization);
  a file claiming v>0 is attested object-for-object against a reference
  built from the shipped DDL — sqlite_master SQL text (covers STRICT,
  CHECKs, and foreign keys) plus per-table `PRAGMA table_info` (name,
  declared type, not-null, default, primary key) and index presence.
  Divergence refuses with structured `LedgerSchemaMismatch { claimed_version,
  violations }` BEFORE any migration runs, so `CREATE TABLE IF NOT EXISTS`
  can never launder an alien or mangled schema into a labeled one. RECOVERY
  TOLERANCE: objects or columns from a LATER version are accepted iff they
  match the current shipped definition exactly — this generalizes the
  historical v2 crash window (committed DDL, stale marker), which still
  heals; incompatible same-name early objects still fail closed.
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
  `in|out` edges), `edge_exists` (exact role-qualified verifier query).
- Streams: `record_metric` (finite REAL only), `append_event` /
  `append_events` (batched, atomic), `tune_put`/`tune_get` (upsert keyed
  kernel × shape-class × machine fingerprint), and `tune_put_if_absent`
  (insert-only conflict preservation for evidence-ledger adoption).
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
  `explain` (strict-JSON causal trees even for hostile artifact-kind text,
  retaining producer outcome and diagnostic, depth-limited, DAG-deduped, loud
  on orphan or malformed input identities), `replay_verdict` (IR,
  all frozen explicits, execution mode, input lineage, outcome, and diagnostic
  must agree; both studies must be drained and finalized before a clean
  verdict is possible, and empty branches refuse because no executed study was
  compared; deterministic ops must then
  reproduce output hashes exactly; fast hash divergences are reported without
  failing; row/branch/session/time envelopes are excluded),
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
  free-riding on `travel`'s forkable worlds. A COMMIT drains first and is the
  v2 domain-separated, length-framed Merkle root of a branch's visible frozen
  ops (leaf = IR + Five Explicits + outcome + diagnostic + execution mode +
  sorted role-qualified linked-artifact hashes; node/root domains and leaf
  count are distinct; wall times, rowids, branch ids, and session envelopes
  are EXCLUDED, so logically identical histories produce identical roots
  across ledgers and runs); unchanged recommits are idempotent (never a
  self-parent cycle), while changed commits chain to their branch predecessor
  and persist as `vcs-commit` events. ENVELOPE VS SEMANTIC IDENTITY (bead
  gp3.17): the Merkle root is the SEMANTIC-STATE identity; the commit
  ENVELOPE is `CommitId { ledger, branch, root }` where the ledger identity
  is persisted on first use as a `vcs-identity` event (copies of a file
  share it — same lineage; independent databases differ; the first event by
  rowid is the authority under concurrent minting). The registry, heads,
  lookup, and checkout are all envelope-keyed, so equal roots reached by
  different branches or ledgers coexist without clobbering
  (`lookup_semantic` lists every envelope sharing a root). CHECKOUT returns
  the exact in-session frozen op/artifact view captured by that commit
  (envelope-scoped, since snapshots carry ledger-local op ids), so later
  ops and later links to an old op cannot leak future artifacts;
  `checkout_delta` compares SEMANTIC LEAF MULTISETS — portable across
  ledger instances and import orders, never local row ids — and reports
  each differing op as `DeltaOp { leaf, local_op }` with the id from that
  side's own ledger (the `perturb()`-style delta a recompute solver
  consumes — nearby checkouts cost |delta|, not |history|). The root binds
  history SEQUENCE; the delta binds the SET: a reordered import has a
  different root but an empty semantic frontier. `merge_views` splits base/only-A/only-B for the
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
or lineage fails), the exact `IntervalOp` (an Add grant cannot authorize
Mul), the claimed color, the color-upgrade scope, a
signer key id, and an expiry day — verified through the caller-
supplied `WaiverVerifier` capability before any write. A grant carries and
signs the full `Color::canonical_bytes` payload (domain-separated signing
encoding v3), not only
the color rank name, so authorization for one interval, validity regime, or
estimator payload cannot authorize another. The in-tree
default is `NoWaiverVerifier` (refuses everything): no cryptographic
capability ships in this crate, so promotion is impossible until a
Franken-compliant signature verifier is wired in (the no-crypto
no-claim). Authentication runs for EVERY `derive_waived` call, including a
claim that does not upgrade rank; choosing the waived path can never turn an
invalid signature into a provenance-bearing grant. Node provenance hashes use
a domain-separated, versioned length-prefixed encoding (v4), including source
versus derived status, the exact operation, and bit-exact canonical color bytes.
The former v3 omitted the operation, and v2 used display-rounded color JSON.
Color-write row schema v2 serializes the canonical v3 signing payload,
signature, key id, node name, parent hashes, exact claimed-color bytes,
operation, scope, and expiry so persisted authorization can be independently
reverified rather than trusting an `authorized:true` assertion. Refusals are structured
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
migration + future-version refusal, the schema-attestation gauntlet
(valid-empty atomic init, conflicting-object-at-v0, partial-schema,
wrong-column, wrong-affinity, and missing-index all refused fail-closed
with the file untouched), dual-path chunked dedupe + round trip,
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
edge cases. `tests/vcs.rs` locks framed, role-qualified, mode-bound commit
identity and proves checkout cannot expose later ops or later artifact links;
in-flight commits are refused. The travel migration battery also reconstructs
the old v2
post-DDL/pre-version-marker crash state and proves reopen heals it without
duplicating columns or losing v1 history.

## Speculation telemetry (bead lmp4.3, schema v3)

Schema v3 adds the `speculation` extension table (uniform Rev S
shape) carrying the four solve-node fields `(proposer_id, accepted,
bound, iterations_saved)` keyed by solve-op identity — ADDITIVE: the
migration regression test proves every pre-existing table still
answers queries. The economics control loop lives in fs-verify
(HELM-side); this ledger stores telemetry, it does not drive solves.

## Three-color write gate (bead qmao.1)

`colors::ColorGraph` is the WRITE-TIME gatekeeper over fs-evidence's
color schema. `source()` accepts Estimated leaves only. A Verified leaf must
carry a `SourceOrigin::Certificate`; the gate reruns `verified_from` and writes
the rederived interval. A Validated leaf must carry an
`SourceOrigin::Anchoring` with dataset content hash and exact regime; the gate
reconstructs the complete color and refuses blank datasets, empty or malformed
regimes, and any claimed dataset/regime drift.
The exceptional source path is an authenticated `WaiverGrant` under the
distinct `source-color` scope and v4 source signing payload. A derive grant
cannot be replayed as source authority.

Derived nodes' colors are COMPUTED from their parents with regime re-checks
against the current execution state. Every parent that exits its regime emits
a retained `ColorDemotion`, keyed by parent position and id in canonical
parent-list order. This remains unambiguous when a parent id occurs twice.
Demotions and typed source origins participate in the domain-separated v6 node
hash. An ordinary explicit claim must equal the exact canonical derived color:
equal rank alone is insufficient because it could narrow an interval, widen a
validity regime, or shrink dispersion; unsupported rank weakening is likewise
refused until a formal weakening relation exists. Claims that outrank the
derivation REFUSE with the capping parents named (the laundering refusal, G3
gauntlet-tested). The derived override is an authenticated `WaiverGrant` under
the `color-upgrade` scope and operation-bound v3 signing payload.

`ColorNode` fields are private and exposed read-only. `ColorGraph::node` is a
checked `Option` lookup, never an indexing panic. `verify_replay()` checks
append ids, backward parent references, canonical demotions, source-origin
rederivation, grant-to-node/lineage binding, ordinary derived colors, and every
node hash. Canonical schema-v3 JSON rows include the typed origin and the exact
v3/v4 signed payload needed for an independent verifier. G3/G5 tests cover
forged positive sources, source/derive grant separation, invalid ids,
multi-parent demotion preservation, deterministic replay, origin substitution,
and signed-payload tampering. Note: this module adds fs-evidence as a runtime
dependency (the colors are its types).

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
- `ColorGraph::verify_replay()` structurally re-earns colors and hashes but does
  not re-authenticate stored waiver signatures: the graph does not retain the
  verifier capability or admission day. The schema-v3 row retains the exact
  signing payload, key id, signature, and expiry for an independent verifier.
  Regime demotion records retain the offending value and are hash-bound, but
  the complete execution-state map is not persisted by this in-memory gate.

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
