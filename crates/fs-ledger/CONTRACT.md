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
  Every byte-returning retrieval re-hashes stored content against its key.
  Retrieval treats signed database integers as untrusted: it never
  preallocates from declared `len`, uses fallible materialization, and performs
  metadata-only SQL preflights for bounded artifact-kind/metadata envelopes,
  BLOB presence, per-row size bounds, dense chunk sequence, count, and checked
  byte totals before selecting variable-size values. Guarded envelope and BLOB
  queries repeat the size predicates so a mutation between preflight and
  materialization still cannot deliver an oversized row. Same-length byte
  tampering is necessarily a late hash failure: streaming callbacks retain
  effects for the delivered prefix, while `get_artifact` returns no bytes.
- Ops/lineage: `begin_op` (validates the Five Explicits field-by-field;
  units travel inside the typed IR, the other four are mandatory columns),
  `finish_op` (exactly once; `ok|error|cancelled`), `op`, `link` (FK-checked
  `in|out` edges), `edge_exists` (exact role-qualified verifier query).
- Streams: `record_metric` (finite REAL only), `append_event` /
  `append_events` (batched, atomic), `tune_put`/`tune_get` (single-statement
  atomic upsert keyed kernel × shape-class × exact machine fingerprint),
  `tune_put_if_absent` (insert-only conflict preservation for evidence-ledger
  adoption), and `tune_rows` (deterministic `(shape_class, machine)` order).
  Tune kernel/shape identities are 1..=64 KiB visible ASCII bytes; machine
  identities are exact opaque 1..=256-byte BLOBs; parameter and measurement
  JSON are each at most 1 MiB. Reads metadata-preflight stored values and
  repeat bounds in guarded payload queries. A kernel scan refuses before
  payload materialization above 1,024 rows or 16 MiB aggregate output.
- Rev S extension tables (sparse v0, uniform `(name UNIQUE, body JSON)`
  shape): `put_extension`/`get_extension` over `requirements`, `model_cards`,
  `evidence`, `scenarios`, `constraints`, `capability_probes`, `imports`,
  `unsafe_capsules`.
- Hygiene: `lint()` (orphan edges/metrics/chunks, artifact and tune storage
  bounds, storage-shape and length invariants, half-finished ops, dangling
  branch references) — all-zero on any healthy or crash-recovered ledger.
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
Non-authorizing annotations still cross an audit boundary: `derive` and replay
require canonical bounded id/signer identities and a non-placeholder reason of
at most 4,096 bytes with control and bidirectional-format characters refused.
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
invalid signature into a provenance-bearing grant. Before verifier dispatch,
machine identities, human reason text, signature presence, claimed-color bytes,
and lineage size are structurally validated and bounded; an accepting callback
cannot authorize malformed epistemic metadata. Verifiers return a sealed,
atomic `PolicyDecision`; the accepting policy fingerprint and historical
admission day travel with the direct grant and every transitive
`WaiverDependency`. Callback panic is a structured fail-closed refusal. Node
provenance hashes use a domain-separated, versioned length-prefixed encoding
(v9): v9 binds color-algebra v2 in both the node domain and canonical color
bytes; v8 added certificate artifact identity, source and waiver policy
fingerprints, admission days, and those fields in every transitive dependency;
v7 added canonical transitive waiver dependencies, v6 added all demotions and operation-correct
grant payloads, v5 added typed source origins, v4 added source/derived status
and the exact operation, v3 added bit-exact canonical color bytes, and v2 used
display-rounded color JSON. Color-write row schema v7 adds the exact canonical
color bytes, exact canonical typed-origin bytes, and node-hash encoding version;
v6 added the explicit color algebra version, and v5 serialized typed origins,
the transitive waiver closure, the canonical v3 derived or v4 source signing
payload, signature, key id, node name, parent hashes, exact claimed-color bytes,
operation, scope, expiry, admission day, and policy fingerprint. Demotion-row
schema v1 stores the exact offending IEEE-754 bits in addition to display JSON.
The in-module row reconstruction test strictly decodes the lowercase exact hex
fields and substitutes them into the same hash reconstruction alongside the
in-memory typed metadata for certificate, anchoring, Estimated, and derived
rows. This proves exact color/origin field fidelity, not that a general
persisted-row importer exists: this in-memory graph still does not parse event
rows, resolve policy fingerprints, or dispatch external authority. Refusals are
structured (`WaiverRejection`:
malformed/bounded field, scope/node/color/lineage mismatch, expiry, policy
refusal, or verifier panic).

## Invariants

1. Artifact identity = BLAKE3 of content; identical bytes dedupe to one row.
   ENVELOPE AGREEMENT (bead gp3.19): the dedupe applies only under an
   AGREEING envelope — the offered `kind` must match exactly, and offered
   metadata (when a claim is made) must canonically equal the stored
   metadata (engine `json()` comparison: whitespace-insensitive, key order
   significant). Offering `meta: None` makes no claim and accepts the
   stored envelope; offering metadata against a row stored without any is
   a conflict. Disagreement refuses with structured
   `LedgerArtifactEnvelopeConflict { hash, field, stored, offered }` at
   EVERY dedupe site (pre-check, concurrent duplicate-key race, streaming
   writer finish) — provenance never depends on insertion order, and
   content identity stays bytes-only (no schema change; byte dedup
   retained). An existing row is fully shape-checked and re-hashed before
   dedupe can succeed, so corruption is never preserved under a successful
   receipt. Concurrent duplicate insert with an agreeing, intact row still
   resolves to dedupe, never an error.
2. Storage shape: inline XOR chunked; `len` always equals stored byte count;
   every inline/chunk BLOB is at most `STORAGE_CHUNK_LEN`; artifact kinds are
   1..=`MAX_ARTIFACT_KIND_BYTES` UTF-8 bytes, metadata JSON is at most
   `MAX_ARTIFACT_META_BYTES` UTF-8 bytes, and chunk `seq` is dense from 0. The
   schema CHECK enforces the storage-form XOR, JSON metadata, and non-negative
   storage counters. Canonical writes enforce the remaining bounds; ordinary
   reads fail closed through metadata-only preflight plus guarded selection,
   and `lint` detects envelope, length, row-bound, count, and sequence violations.
3. Tune rows have canonical, bounded envelopes: kernel and shape-class are
   non-empty visible ASCII within `MAX_TUNE_KERNEL_BYTES` and
   `MAX_TUNE_SHAPE_CLASS_BYTES`; the opaque machine BLOB is non-empty and at
   most `MAX_TUNE_MACHINE_BYTES`; params/measured are valid JSON within their
   1 MiB bounds. Both write APIs share this admission gate. `tune_get` and
   `tune_rows` inspect type, exact BLOB byte length, and JSON validity before
   payload selection; guarded queries repeat those predicates. `tune_rows`
   additionally caps rows and checked aggregate output bytes before selecting
   variable-size values. Raw-SQL rows outside this contract fail closed as
   `TuneCorrupt`; valid but excessive histories fail as `TuneReadLimit`.
4. Ops are event-sourced facts: `(t_end IS NULL) = (outcome IS NULL)` is a
   table CHECK; an op finishes at most once (`DoubleFinish` otherwise).
5. Edges only reference existing ops and artifacts (enforced FKs).
6. A crash-recovered ledger lints clean: transactions make op+edges+metric
   groups all-or-nothing (kill -9 battery, `ledger_007`).
7. Wall-clock timestamps are provenance envelope, never content identity.

## Error model

All fallible APIs return `LedgerError` — structured variants with stable
`code()` strings and actionable Display text: `Open`, `FutureSchema` (newer
file refused, never clobbered), `Sql`, `Busy` (retryable contention —
busy/locked/write-conflict; retry with backoff), `MissingExplicit` (names the
offending Five Explicits field), `Invalid` (names the field),
`Corrupt`, `TuneCorrupt` (stored tune envelope violated), `TuneReadLimit`
(bounded scan refused), `NotFound`, `DoubleFinish`, `WriterInTransaction`.
Never panics across the crate boundary. Signed database metadata that
represents a length or count is converted with an explicit non-negative check;
physical corruption cannot reinterpret `-1` as `u64::MAX`.

## Determinism class

Content hashing is bit-stable across runs, thread counts, and ISAs (pure
function). Row ids, timestamps, and physical file bytes are NOT deterministic
and are excluded from identity. Deterministic replays should pass logical
times to `begin_op`/`append_event` (caller-controlled `t`). Tune scans use the
total order `(shape_class, machine)` and refuse rather than truncate.

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
throughput smoke ledgered as a metric. `ledger_013` races two file-backed
connections through the same atomic tune upsert after forcing tune-table leaf
splits, then proves a single untorn row, bounded scan, clean lint, and identical
reopen. Tune unit regressions cover every exact field limit and limit+1,
empty/NUL/non-ASCII identities, hostile oversized raw-SQL rows, lint detection,
and deterministic row/aggregate scan caps. `tests/travel.rs`: genuine-v1 →
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
color schema. Every public write rejects a blank, padded, placeholder,
control-bearing, non-canonical, or oversized node name before color work,
authority callbacks, hashing, row formatting, or cloning. The public
`MAX_COLOR_NODE_NAME_BYTES` limit equals the shared
`fs_evidence::MAX_COLOR_IDENTITY_BYTES` bound, and replay rechecks the same
grammar. `source()` accepts Estimated leaves only, rejecting blank or
placeholder estimator identities, surrounding identity whitespace, identities
longer than `fs_evidence::MAX_COLOR_IDENTITY_BYTES`, and NaN/negative
dispersion (positive infinity remains the explicit
no-spread-claim sentinel). The reserved `derived:v2:` identity namespace cannot
enter through `source()` because such diagnostics and compositions require
retained parent lineage. fs-evidence composition uses the same owner constant
and domain-separated, length-framed bounded identities; replay calls the
fs-evidence demotion-identity helper instead of duplicating that grammar. Thus
legitimate long pipelines do not fail only because provenance names grew by
concatenation. A Verified leaf must carry a
`SourceOrigin::Certificate` with the retained certificate artifact's content
hash; the gate reruns `verified_from` and writes the rederived interval. A
Validated leaf must carry an
`SourceOrigin::Anchoring` with dataset content hash and exact regime; the gate
reconstructs the complete color and refuses blank, placeholder, or padded
producer/dataset/axis identities, empty or malformed regimes, and any claimed
dataset/regime drift. Validity boxes are bounded at 1,024 axes. A multi-parent
fold first merges the distinct axes of effective Validated parents into a
bounded preflight map; the 1,025th axis refuses before parent colors are cloned
or an oversized intersection is constructed. A regime exit or Estimated parent
short-circuits that axis work because Estimated already absorbs the fold.
Claimed colors and origin regime counts are validated before origin cloning or
canonical-byte comparison, and structural validation completes before the
source authority is invoked. Shape-valid public
fields are
not authority: `source_with_origin` also requires an injected
`SourceOriginVerifier`. Its read-only request and canonical payload cover the
node name, exact claimed color, certificate artifact hash, and every other
certificate/anchor field. Its sealed `PolicyDecision` returns acceptance and
the exact policy fingerprint atomically; that fingerprint is hash-bound into
the node and row. Verifier panic fails closed before append;
`NoSourceOriginVerifier` is the fail-closed default.
The exceptional source path is an authenticated `WaiverGrant` under the
distinct `source-color` scope and v4 source signing payload. A derive grant
cannot be replayed as source authority. Authentication does not bypass payload
structure: the shared color validator rejects NaN or inverted Verified
intervals, invalid Validated identities/regimes, and invalid Estimated
identities/dispersion before either `source_waived` or `derive_waived` can
append a node. Ordinary composition is checked by the same validator before
append. Ordered infinite endpoints remain a sound but vacuous enclosure and
must not be mistaken for decision-grade tightness. A waiver authorizes claim
policy, never malformed epistemic data.

Derived nodes' colors are COMPUTED from their parents with regime re-checks
against the current execution state. Every parent that exits its regime emits
a retained `ColorDemotion`, keyed by parent position and id in canonical
parent-list order. This remains unambiguous when a parent id occurs twice.
Demotions, typed source origins and their admitting policies, and canonical
transitive waiver dependencies with their policies/admission days participate
in the domain-separated v9 node hash. An ordinary explicit claim
must equal the exact canonical derived color:
equal rank alone is insufficient because it could narrow an interval, widen a
validity regime, or shrink dispersion; unsupported rank weakening is likewise
refused until a formal weakening relation exists. Claims that outrank the
derivation REFUSE with the capping parents named (the laundering refusal, G3
gauntlet-tested). The derived override is an authenticated `WaiverGrant` under
the `color-upgrade` scope and operation-bound v3 signing payload.

`ColorNode` fields are private and exposed read-only. `ColorGraph::node` is a
checked `Option` lookup, never an indexing panic. `verify_replay()` checks
append ids, backward parent references, canonical demotions, source-origin
rederivation, structural color validity (including waived nodes),
grant-to-node/lineage binding, ordinary derived colors, and every node hash.
Ordinary and waived derivations retain the sorted, duplicate-free union of
their parents' historical waiver dependencies plus each parent's own grant;
fan-in and retained authority closure have both count and aggregate-byte limits
before cloning, hex serialization, or append. The current closure cap is
`MAX_WAIVER_CLOSURE_BYTES` (8 MiB), in addition to 1,024 distinct authorities.
`scientific_color()` returns `None` for every directly or transitively
waived node, while raw declaration inspection is deliberately named
`declared_color_unverified()`. `waiver_dependencies()` and
`depends_on_waiver()` expose the exact reason for that refusal. Replay resolves each
dependency to its earlier authorizing node and recomputes the exact
parent-derived closure, including the original policy fingerprint and
admission day. Canonical schema-v7 color rows include the color-algebra and
node-hash versions, exact canonical color/origin bytes, typed origin and
certificate artifact identity, direct policy/admission context, transitive
dependencies, and the exact v3/v4 signed payload needed for an independent
verifier. G3/G5 tests cover
forged positive sources, source/derive grant separation, invalid ids,
multi-parent demotion preservation, deterministic replay, origin substitution,
policy and certificate-artifact substitution, callback panic, composed-bound
overflow, signed-payload tampering, padded source identities, and authenticated
attempts to admit malformed colors. Invalid node-name tests cover empty,
control-bearing, placeholder, and oversized inputs and prove rejection occurs
before an injected callback is invoked. The aggregate-axis regression drives
all 1,024 parent slots with 1,025 distinct axes and also proves that a prior
regime demotion correctly avoids a false-positive refusal. Note: this module
adds fs-evidence as a runtime dependency (the colors are its types).

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
- The v1 tables do not encode the per-BLOB `STORAGE_CHUNK_LEN` or artifact
  envelope byte bounds as schema CHECKs. Existing databases are protected by
  bounded canonical writes, metadata-only read preflights, guarded variable-size
  queries, and lint; resistance to a hostile client executing arbitrary SQL is
  not claimed as a DDL property.
- The v1 `tune` DDL checks JSON syntax but does not encode the canonical
  identity, machine, JSON byte, row-count, or scan-byte bounds. The public API
  enforces them on writes, metadata-preflights and guards reads, and reports
  bounded envelope violations through lint; arbitrary raw SQL is therefore
  detected and refused, not prevented as a DDL property.
- `ColorGraph::verify_replay()` structurally re-earns colors and hashes but does
  not itself re-run external source-origin or waiver capabilities. It retains
  the complete request/artifact fields, exact policy fingerprints, waiver
  admission day, signing payload, key id, signature, and expiry so an
  independent verifier can resolve the named policy and re-authenticate.
  Replay re-applies the exact Estimated-leaf identity and annotation-validation
  rules and refuses orphan human waiver annotations on source leaves. Bounded,
  audit-safe annotations on real derived operations retain their documented
  non-authorizing meaning. Regime demotion records retain the offending value
  and its exact IEEE-754 bits and are hash-bound, but
  the complete execution-state map is not persisted by this in-memory gate.
- Waiver expiry is checked at the authorizing node's admission day. Descendants
  preserve that historical grant and remain tainted indefinitely; they do not
  silently renew it, and `verify_replay()` has no caller-supplied current day
  with which to make a new policy decision.
- Transitive waiver visibility currently stores the complete unique grant
  closure on each descendant, bounded by `MAX_WAIVER_DEPENDENCIES` and
  `MAX_WAIVER_CLOSURE_BYTES`. This is deliberately inspectable and replayable;
  compact/sublinear waiver-lineage storage and a general persisted-row importer
  are not claimed.

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
