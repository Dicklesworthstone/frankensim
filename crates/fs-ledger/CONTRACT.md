# CONTRACT: fs-ledger

> Status: ACTIVE (Design Ledger, schema v12). Owns the core schema + Rev S
> extension tables, BLAKE3 content addressing, the WAL/snapshot concurrency
> contract, and — since schema v2 — forkable worlds, `at(t)` views,
> `explain()`, the replay audit, and unreferenced-artifact GC (`travel`
> module).

## Purpose and layer

The Design Ledger (plan §11.2, Bet 10): FrankenSQLite-backed system of record
for content-addressed artifacts, event-sourced ops with the frozen Five
Explicits, lineage edges, metric time series, the autotuner cache, and the
fine-grained event stream. Layer: L6 (HELM). Runtime deps: `std`, `fsqlite`,
and the lower-layer Franken crates declared in `Cargo.toml`, including
`fs-exec` for opaque drain/finalize reports and snapshot-envelope validation.

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
- `LedgerInstanceId` / `Ledger::instance_id()` — opaque, move-stable identity
  of one physical ledger. Schema v4 stores exactly one 16-byte UUID in
  `ledger_identity`; fresh initialization and v1-v3 migration seed it inside
  the same transaction as the v4 version marker. Schema v5 adds attested
  update/delete refusal triggers plus a reinsert guard that permits only the
  initial empty-table seed. Reopenings and path aliases of one file
  agree, while replacement files at the same path and independent in-memory
  handles differ. `instance_id()` is the cached open-time value;
  `checked_instance_id()` first performs one bounded `sqlite_master` lookup
  for the three reserved guard names and requires their normalized definitions
  to equal the shipped v5 DDL, then re-reads the current row and refuses a
  missing, malformed, changed, or extra identity row. Row validation reads at
  most two rows across the whole table before requiring exactly one
  `singleton = 1` row, so constraint-bypassed rows cannot hide outside a
  filtered query or force an unbounded corruption scan. `lint()` performs the
  same checked boundary; a missing or weakened guard refuses even while the
  identity row remains unchanged. A v4+ ledger with missing or malformed
  identity refuses open rather than silently rotating authority. First
  initialization refreshes `user_version` after beginning its transaction, so
  a peer that canonically initialized after the preflight v0 read is attested
  and accepted instead of being misclassified as alien unversioned state.
- `ContentHash`, `Blake3`, `hash_bytes` — in-house BLAKE3 (plain hash mode,
  32-byte output), pure safe Rust; artifact identity everywhere. The
  implementation is OWNED by the UTIL crate `fs-blake3` (bead 7uq9) and
  re-exported here unchanged — same paths, same bits.
- Artifacts: `put_artifact` (≤ `STORAGE_CHUNK_LEN` inline; larger stored as
  `artifact_chunks` rows because fsqlite has no incremental-blob API),
  `ArtifactWriter` (streaming; hashes incrementally, stages chunks under a
  provisional key inside a writer-owned transaction, promotes on `finish`),
  `get_artifact` / `read_artifact_chunks` / `artifact_info`, plus
  `get_artifact_bounded` / `read_artifact_chunks_bounded` for consumers that
  must refuse above a caller-supplied payload cap before any byte callback or
  result-buffer allocation,
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
- Ops/lineage: `begin_op` (validates the Five Explicits field-by-field and
  byte-bounds session, IR, seed, versions, budget, and capability fields at
  1 MiB each before JSON validation and stores execution mode only
  as the exact `deterministic|fast` enum;
  units travel inside the typed IR, the other four are mandatory columns),
  `finish_op` (exactly once; `ok|error|cancelled`; optional diagnostic JSON is
  bounded at 1 MiB before validation), `op` (metadata-only
  type/length/CASE-gated JSON preflight followed by the same guarded payload
  query), `op_execution_context` (fixed-size typed branch/mode read after the
  same op-envelope preflight), `link` (FK-checked `in|out` edges accepted only
  while the target op is bounded, branch-valid, and unfinished),
  `edge_exists` (exact role-qualified verifier query), plus
  `artifact_producer_ops_bounded` and `op_artifact_edges_bounded`. The bounded
  lineage reads accept caller caps through 1,024 rows, issue `LIMIT cap+1`
  through the schema-v9 covering indexes `(artifact, role, op)` and
  `(op, role, artifact)`, return only the capped deterministic prefix, and
  expose `truncated` so verifier paths can refuse extra producer/edge fan-out
  without scanning, sorting, or materializing an unbounded DAG. Selected edge
  values are SQL-CASE-sanitized to fixed role/hash envelopes before Rust can
  materialize them. A zero cap is a bounded existence probe. `link` uses one
  conditional `INSERT ... SELECT`, so it serializes with `finish_op`: an edge
  either lands before the terminal outcome or changes zero rows and returns
  `OpLineageSealed`. Missing, corrupt, and already-finished ops remain distinct
  structured refusals, and the zero-row path never inserts partial lineage.
  `seal_artifact_output` atomically binds an artifact to an already-existing
  sole output producer; the immutable `artifact_output_seals` row and attested
  edge triggers then reject every different producer. Exact same-op sealing is
  idempotent, input reuse remains legal, and `artifact_output_seal` is the
  fixed-size verifier read. `seal_op_artifact_edges` independently freezes one
  operation's complete bounded edge set after an exact-cardinality `cap+1`
  probe. Its immutable `op_artifact_edge_seals` row blocks every later edge
  insert, update, or delete for that op; `op_artifact_edge_seal` revalidates the
  stored count with at most `count+1` covering-index rows before returning it.
  Both seal accessors verify their parent edge/op state, so constraint-bypassed
  orphan or fan-out corruption fails closed rather than becoming idempotent.
  Migration to schema v9 deliberately leaves both seal tables empty: exclusive
  provenance is a consumer claim that the ledger cannot infer from historic
  edges alone. Consumers may atomically adopt missing seals only after
  revalidating their complete bounded historic lineage; conflicting seals are
  immutable and must fail closed.
- Streams: `record_metric` (finite REAL only), `append_event` /
  `append_events` (batched, atomic), `tune_put`/`tune_get` (single-statement
  atomic upsert keyed kernel × shape-class × exact machine fingerprint),
  `tune_put_if_absent` (insert-only conflict preservation for evidence-ledger
  adoption), and `tune_rows` (deterministic `(shape_class, machine)` order).
  Tune kernel/shape identities are 1..=64 KiB visible ASCII bytes; machine
  identities are exact opaque 1..=256-byte BLOBs; parameter and measurement
  JSON are each at most 1 MiB. Reads metadata-preflight stored values and
  repeat bounds in guarded payload queries; nested SQL `CASE` gates JSON
  validation behind type and byte-length checks, so an oversized raw row is
  refused without reparsing its payload. A kernel scan refuses before payload
  materialization above 1,024 rows or 16 MiB aggregate output.
- Durable session registry (`session_registry`, schema v6-v8):
  `SessionMutationClaim` binds one mutation authority to the checked physical
  `LedgerInstanceId`, durable governor, session-open authority, kind, session,
  exact scope, generation, optional causal ordinal, and exact payload bytes.
  `claim_session_mutation` commits that claim before caller work and returns
  exactly one of fresh `Claimed { permit }`, verified `Pending`, or verified
  `Terminal`; only the fresh caller receives the sealed positive permit needed
  to terminalize an existing Pending claim. New `submission` claims must carry
  a unique governor/kind admission ordinal in `1..=i64::MAX` and cannot bypass
  this preclaim boundary. `append_session_terminal_batch` atomically commits
  each typed receipt, its dense authority-owned global audit events, and a
  deterministic complete ordered batch witness. Exact batch replay appends
  nothing. A mixed retry may re-witness an existing terminal beside a new one,
  so one terminal may have up to 1,024 distinct witnesses; every witness is
  rehashed over the complete member preimage and totals. Reads verify bounded
  storage types, claim/payload/receipt hashes, dense event ownership, global
  event bytes, every batch marker, and every complete membership preimage.
  The reciprocal generation fence rejects a new old-generation submission
  after a terminal successor pause acknowledgement and rejects a pause terminal
  while an omitted draining-generation submission remains Pending; partial or
  corrupt terminal-looking storage fails closed. Generation recovery probes
  use indexed keyset pages, verify terminal witnesses rather than trusting raw
  row presence, and cap inspection at 8,192 claims generally and 4,096
  submission predecessors for a pause. One batch is capped at 1,024 terminals,
  1,024 events, and 4 MiB encoded bytes; claim payloads and terminal receipts
  are each capped at 1 MiB. Schema v8 mirrors every claim's bounded discovery
  envelope in an independently indexed immutable witness. Exact authority
  reads require one row in each table and compare every copied field after
  authenticating the primary claim hash; filtered recovery and generation
  fences take a deduplicated union of both indexes before authentication.
  Single-table deletion, key drift, or semantic-column corruption therefore
  cannot turn a dangerous claim into a trusted negative lookup.
  `session_governor_claim_snapshot` strengthens the restart boundary from a
  count to authenticated membership: under one self-owned stable read
  transaction it keyset-scans the unfiltered union of all authorities, applies
  SQL type/length/canonical-text guards before materializing each compact
  envelope, independently recomputes both copies' complete claim-hash preimage
  from their stored payload hash, requires exact cross-copy equality, and only
  then selects the requested governor. The private-field snapshot binds the
  checked physical ledger, governor, count, and a domain-separated root over
  sorted authority bytes. The compatibility count accessor delegates to this
  full scan. A bounded `OFFSET 8192` probe first caps each physical mirror
  before malformed-key aggregates or multiplicity checks; the deduplicated
  union independently admits at most 8,192 total authorities, pages in
  1,024-authority chunks, and refuses authority 8,193 before issuing either of
  its compact-envelope queries. An exact `N`-authority scan consumes the
  checked-identity query, two preflight queries, `max(1, ceil(N/1024))` page
  queries, and four compact-copy queries per authority. No partial snapshot
  escapes on any refusal. An already-open caller transaction is refused
  without being committed or rolled back. Migration from v7
  backfills inside the version transaction and hash-verifies every surviving
  source claim before publishing v8. V8 also splits the two OR-based immutable
  reinsert guards into one point lookup per unique key, avoiding dependence on
  multi-index-OR planning at exact read-cap fixture scale.
- Solver checkpoint receipts (`session_registry`, schema v10):
  `attest_solver_checkpoint` accepts an existing artifact only when its exact
  kind is `solver-state`, its declared size is at most 64 MiB, its complete
  bytes pass generic fs-exec snapshot-envelope length/checksum validation, and
  its provenance equals an opaque executor-minted `DrainFinalizeReport` run.
  The immutable private-field `SolverCheckpointReceipt` binds physical ledger,
  session, run, exact pause authority, gate generation, artifact hash, drain
  report hash, and registered/drained counts. One pause authority owns one
  receipt; exact retry after response loss returns it, while changed fields
  conflict atomically. Fixed-width transport decoding earns only a candidate:
  `verify_solver_checkpoint_receipt` must revalidate identity, physical-ledger
  membership, immutable row equality, artifact integrity, and run provenance.
- Quantity dimension crosswalks (`crosswalk`, schema v11):
  `record_qty_dimension_crosswalk` persists only fs-qty's typed
  legacy-five-to-six `AppendMoleZero` receipt. Both exact JSON artifacts must
  already exist, each is admitted under a 4 KiB materialization budget, and
  the ledger independently decodes the historical bytes, reproduces the
  canonical target, and compares every receipt field before one immutable row
  commits. An exact retry is idempotent; a different target for the same old
  hash refuses. `qty_dimension_crosswalk` repeats artifact integrity and
  fs-qty replay before returning evidence, so row presence alone is never a
  migration claim. Migration from v10 creates an empty table and infers no
  historical rows.
- Semantic state checkpoints (`state_checkpoint`, schema v12):
  `record_state_checkpoint` mints one portable private-field receipt over a
  typed state-slot identity, caller-asserted fs-matdb law/version/state-schema
  metadata and canonical parameter-block hash, one retained
  `constitutive-runtime-state` artifact, and an injected L3 contract/code
  identity. The identity preimage is the receipt version as little-endian
  `u32`, the 32-byte state slot, a little-endian `u16` law-id byte count plus
  exact UTF-8, both semantic versions as little-endian `u32`, and the three
  32-byte artifact/parameter/code hashes in canonical order; the portable
  transport appends the resulting 32-byte identity. Exact retry dedupes by
  receipt hash; one stable slot may accumulate
  successive immutable checkpoints. `load_state_checkpoint` accepts the
  receipt hash plus caller-known law/parameter/code semantics, compares the
  complete semantic tuple before materializing state bytes, and returns a
  structured `LedgerUnknownStateSemantics` refusal naming the stored and expected law,
  versions, parameter hash, and code hash on disagreement. Runtime-state reads
  are integrity checked under a 64 MiB cap. Migration infers no semantic
  receipts from generic solver snapshots.
- Rev S extension tables (sparse v0, uniform `(name UNIQUE, body JSON)`
  shape): `put_extension`/`get_extension` over `requirements`, `model_cards`,
  `evidence`, `scenarios`, `constraints`, `capability_probes`, `imports`,
  `unsafe_capsules`.
- Storage-envelope hygiene: `lint()` (orphan edges/metrics/chunks/crosswalk artifacts; artifact,
  op, tune, crosswalk, and semantic-checkpoint storage bounds; storage-shape
  and length invariants; half-finished ops; dangling branch references) —
  all-zero on any healthy or crash-recovered ledger). Lint deliberately does
  not replace semantic receipt replay or full artifact re-hashing; callers use
  the receipt verifiers and integrity scan for those stronger claims.
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
  `gc_unreferenced_artifacts` (artifacts with neither a lineage edge, a
  solver-checkpoint receipt, either side of a quantity dimension crosswalk,
  nor a semantic state-checkpoint receipt; every supported root makes an
  artifact immortal).

Schema divergences from plan Appendix D, all deliberate: `JSON` columns are
STRICT-legal `TEXT` with `json_valid()` CHECKs (Appendix D as written is not
valid STRICT SQL), and `artifacts` gains `len`/`chunk_count` +
`artifact_chunks` for bounded-memory large-field storage. Schema v4 adds the
singleton `ledger_identity` table so higher-layer sink authority is tied to
the database instance rather than a path string or Rust object address;
schema v5 makes ordinary SQL update/delete/replacement mutation of that row
fail closed. Schema v6 adds immutable session claims, terminal receipts,
owned-event links, and flush-batch witnesses. Schema v7 appends causal-ordinal
ownership/range indexes and insert guards without rewriting the shipped v6
tables. Schema v8 adds the immutable dual-copy claim-discovery witness and
independently indexed reinsert guards; genuine v7 backfill is authenticated
before its marker advances. The tracked v6 table shape had no wired public
registry writer; the v2
batch/event hash domains in `session_registry` are therefore the first
supported writer format. NULL submission ordinals in immutable v6-shaped rows
are read only as defensive compatibility: Pending remains indeterminate and a
terminal consumer must recover its authenticated ordinal from the receipt.
Schema v9 adds immutable sole-producer and exact-operation-edge-set seals.
Schema v10 adds the immutable `session_checkpoint_receipts` table, one unique
pause-authority index, artifact foreign key, and update/delete/reinsert guards;
migration infers no receipts from historic free-form acknowledgement payloads.
Schema v11 adds immutable `qty_dimension_crosswalks` rows with source and
target artifact foreign keys, fixed v1-to-v2 `AppendMoleZero` semantics, a
target lookup index, and update/delete/reinsert guards. Migration deliberately
does not infer semantic evidence from old artifact pairs.
Schema v12 adds portable immutable `semantic_state_checkpoint_receipts`,
bounded UTF-8 law-id bytes, exact u32 law/state-schema versions, a runtime-state
artifact foreign key, parameter and contract/code hashes, indexes by stable
slot/law/artifact, and update/delete/reinsert guards. State slots are indexed
but not unique so successive immutable snapshots do not overwrite history.

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
`verify_color_row_stream` is the independent persisted-row verifier: it parses
strict JSON without any `ColorGraph`, `ColorNode`, or in-memory `Color`, resolves
parent hashes only from earlier accepted rows, folds immediately preceding
demotion rows in canonical parent-position order, reconstructs signed waiver
payloads from their persisted fields, then rebuilds the exact v9 node preimage
and compares the BLAKE3 hash. It accepts only color-write schema v7, demotion
schema v1, node-hash v9, and color algebra v2; older/future formats are explicit
structured refusals, never display-JSON fallbacks. This re-earns persisted hash
identity; it does not rerun source/waiver authorities or import rows back into a
scientifically admitted graph. Refusals are
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
   payload selection; guarded queries repeat those predicates, and JSON checks
   are `CASE`-guarded so evaluation order cannot parse an over-limit value.
   `tune_rows` additionally caps rows and checked aggregate output bytes before
   selecting variable-size values. Raw-SQL rows outside this contract fail
   closed as `TuneCorrupt`; valid but excessive histories fail as
   `TuneReadLimit`.
4. Ops have canonical bounded envelopes: optional session, IR, seed, versions,
   budget, capability, and optional diagnostic are each at most 1 MiB; required
   BLOB/JSON fields are nonempty, JSON fields are valid, and execution mode is
   exactly `deterministic` or `fast`. Canonical writes
   apply byte limits before JSON parsing. `op` first reads only storage types,
   byte lengths, and nested-CASE JSON-validity bits, then repeats every bound in
   the payload query; hostile raw rows fail as `OpCorrupt` without materializing
   their variable-size values. `lint().malformed_ops` reports the same contract.
5. Ops are event-sourced facts: `(t_end IS NULL) = (outcome IS NULL)` is a
   table CHECK; an op finishes at most once (`DoubleFinish` otherwise).
6. `finish_op` is the public lineage-sealing transition. `link` and
   `finish_op` race through database-serialized write statements: a link wins
   wholly before finish, or observes zero changed rows after finish and returns
   `OpLineageSealed`. No public API removes or role-shifts an edge, so after a
   successful finish the role-qualified edge set is stable.
7. Edges only reference existing ops and artifacts (enforced FKs).
8. A crash-recovered ledger lints clean: transactions make op+edges+metric
   groups all-or-nothing (kill -9 battery, `ledger_007`).
9. Wall-clock timestamps are provenance envelope, never content identity.
10. Physical ledger identity is a persisted opaque 16-byte UUID. It survives
   handle moves, file aliases, and reopenings, but never transfers to a new
   database merely because that database occupies the same path. New
   identities use 122 bits from the operating system's `/dev/urandom` source
   on supported Unix targets, with RFC 4122 version/variant bits overlaid.
   Schema v5 refuses every UPDATE, DELETE, or non-initial INSERT through
   attested triggers. The checked accessor re-attests all three exact trigger
   definitions before trusting the row and detects drift against an already-open
   handle after a bypassed guard is restored.
11. A durable session terminal is valid only as the conjunction of its exact
    immutable claim, receipt hash, dense owned-event sequence, rejoined global
    event bytes, and at least one complete authenticated batch witness. Claim,
    terminal, event ownership, batch marker, and batch membership commit in one
    transaction. Missing, extra, reordered, foreign-ledger, future-schema, or
    hash-mismatched state is corruption; raw terminal-row presence never proves
    completion. Recovered Pending work is explicitly indeterminate and receives
    no terminalization permit.
12. A solver checkpoint receipt is valid only as the conjunction of its
    domain-separated fixed-field identity, one immutable row in the checked
    physical ledger, a present `solver-state` artifact whose complete bounded
    bytes re-hash correctly, a valid generic snapshot envelope, matching run
    provenance, and equal nonzero registered/drained worker counts. Row or
    artifact presence alone never proves pause completion.
13. A quantity dimension crosswalk is valid only when its immutable row names
    two retained, integrity-clean artifacts and replaying the exact historical
    JSON through fs-qty reproduces the same typed receipt and exact canonical
    six-base bytes. The two byte identities remain distinct even when multiple
    historical spellings converge on one canonical target.
14. Semantic runtime-state bytes are replayable only when their immutable
    receipt identity and row agree, the retained artifact re-hashes under the
    explicit read cap, and the caller supplies the exact stored law id, law
    version, state-schema version, canonical parameter-block hash, and
    contract/code hash. Any disagreement withholds the state bytes as unknown
    semantics rather than attempting a stale decode.

## Error model

All fallible APIs return `LedgerError` — structured variants with stable
`code()` strings and actionable Display text: `Open`, `FutureSchema` (newer
file refused, never clobbered), `Sql`, `Busy` (retryable contention —
busy/locked/write-conflict; retry with backoff), `MissingExplicit` (names the
offending Five Explicits field), `Invalid` (names the field),
`Corrupt`, `TuneCorrupt` (stored tune envelope violated), `TuneReadLimit`
(bounded scan refused), `NotFound`, `DoubleFinish`, `WriterInTransaction`.
`ArtifactReadLimit` refuses an artifact whose stored metadata declares a length
above the caller's explicit validation/materialization budget before payload
delivery; it makes no independent content-integrity claim.
`OpCorrupt` refuses a stored op envelope that violates its type, byte, JSON,
finish-state, branch, execution-context, or role-qualified edge contract before
materialization. `OpLineageSealed` is the typed, non-retryable refusal for any
edge offered after `finish_op`; unknown ops remain `NotFound`, so callers never
confuse missing work with immutable terminal provenance. `Invalid { field:
"cap", .. }` refuses a bounded lineage cap
above 1,024 before issuing SQL; malformed producer identities surface as
`Corrupt`. `Invalid { field: "artifact_output_seal", .. }` refuses a missing,
ambiguous, or conflicting producer at seal time, while `Invalid { field:
"op_artifact_edge_seal", .. }` refuses an excessive, missing, or mismatched
edge-set seal. `Invalid { field: "edge", .. }` refuses an output link that
conflicts with an immutable producer seal or any link to a sealed operation.
`InstanceIdentityCorrupt` refuses a v4+ database whose singleton identity is
missing, malformed, or differs from the handle's cached open-time authority.
`InstanceIdentityUnavailable` refuses to mint a new identity when the safe
std-only OS entropy source is unavailable; it never falls back to process ids,
timestamps, addresses, or counters.
Checkpoint attestation and verification report malformed artifacts, foreign or
missing receipts, run disagreement, drain-count disagreement, response-loss
conflicts, and transport identity mismatch as fail-closed `Invalid`,
`ArtifactReadLimit`, `Corrupt`, or underlying storage errors; none is silently
treated as a negative lookup or completion.
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
and deterministic row/aggregate scan caps (including an exact 16 MiB boundary
that counts the cloned kernel identity once per returned row).
Artifact unit regressions cover inline and chunked exact caller caps, cap+1
refusal with zero payload callbacks, and the explicit metadata-declaration
precedence for a tampered length. Op unit regressions cover exact and cap+1 canonical writes
for every variable-size field, raw-SQL oversized IR/version rows, guarded read
refusal, typed execution-context corruption/missing behavior, deterministic
role-qualified lineage ordering, exact-cap/cap+1/zero-cap truncation, explicit
covering-index query plans without temporary sorts, hostile edge-identity
sanitization, typed missing/corrupt/finished link refusal with zero row changes,
caller-transaction rollback, a real two-connection link/finish ordering race,
immutable sole-producer and exact-op-edge-set seals, raw trigger and orphan
detection, a real two-connection seal/link race, v8-to-v9 migration
including stale-marker healing, and `malformed_ops` lint detection.
The `ledger_003b`/`ledger_003c`/`ledger_003d`/`ledger_003e` identity battery covers handle
movement, independent memory ledgers, file reopen and aliasing, same-path file
replacement, genuine v3 and v4 migrations, UUID shape, v5 update/delete/insert
refusal for valid UUID-shaped replacements, fail-closed missing/malformed
identity without advancing the marker, and checked old-handle/lint refusal
after a deliberate DDL bypass plus restoration of the shipped trigger. It also
drops each guard while leaving the identity row unchanged and substitutes a
same-name weakened definition, proving that checked authority and lint require
the exact shipped guard set. The inline G4 first-open barrier fixes one opener's
v0 observation before a peer initializes and proves the stale observer accepts
the peer's attested schema and physical identity.
`tests/session_registry.rs` covers preclaim/Pending/terminal state, exact and
mixed-batch replay, submission admission ownership, reciprocal pause fences,
real-file reopen, foreign-ledger and altered-byte conflicts, exact cap/limit+1
claim/receipt/event/batch byte budgets, and transaction rollback. Nested
registry tests use deliberate in-memory trigger/table bypasses to prove future
schema, hash, event-link, batch-membership, batch-total, and partial-terminal
corruption fails closed, including both directions of the generation fence.
They also prove claim-side, discovery-side, and missing-witness corruption
cannot hide a row from filtered recovery. Inline restart-snapshot regressions
cover empty and exact membership, same-count/different-member roots,
concordant two-table governor rekeying, deterministic lowest-authority
corruption, unchanged refusal of a caller-owned transaction, exact-cap
pagination, and cap+1 refusal. The migration battery accepts an authenticated
genuine-v7 claim, heals exact v8 objects under a stale v7 marker, and rolls back
without advancing when a v7 claim's semantic bytes no longer match its hash.
Canonical bulk fixtures exercise the exact and limit+1 read boundaries for the
8,192-claim recovery probe, 4,096-submission pause fence, and 1,024-witness
terminal lookup without weakening the production constants under test.
Nested registry tests also cover checkpoint response-loss idempotency,
fixed-width transport forgery, foreign-ledger refusal, conflicting artifact
replay, wrong artifact kind, and snapshot/run mismatch.
`tests/color_battery.rs` `col-018` freezes exact canonical-byte sentinels for
positive/negative zero, infinite dispersion, infinite interval endpoints, a
maximum-length identity, and validated/derived rows; it independently rehashes
the full row stream and proves schema/hash/algebra version refusal plus exact-byte
tamper detection. Existing signed-grant, transitive-waiver, hostile-JSON, and
non-finite demotion cases also pass through the public row-only verifier.
`tests/travel.rs`: genuine-v1 →
v2 migration with history intact, fork storage audit (N forks = 1× artifacts
+ deltas) + branch independence, replay audit battery (clean /
deterministic-failure / fast-divergence), explain() full-lineage
reconstruction with loud orphan-input failure, at(t) monotone mid-sweep
consistency, and a kill -9 battery during fork traffic. Unit tests in
`src/lib.rs`, `src/hash.rs`, and `src/travel.rs` cover the API surface and
edge cases. `tests/vcs.rs` locks framed, role-qualified, mode-bound commit
identity and proves a terminal op refuses later artifact links while checkout
remains on its frozen commit; in-flight commits are refused. `tests/travel.rs`
also proves a refused late link leaves the `explain()` causal tree byte-stable.
The travel migration battery also reconstructs
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
verifier. The row-only verifier treats those exact lowercase byte fields as the
hash preimage and rejects duplicate JSON keys, sparse/out-of-order node ids,
forward parents, mismatched demotion positions, non-canonical hex, malformed
grant payload reconstruction, unsupported versions, and hash divergence.
G3/G5 tests cover
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

## Color admission authority (bead 6pf9, stage S1)

The graph is the minting authority for `fs_evidence::AdmittedColor`:

- `ColorGraph::admission_receipt(id)` mints an `AdmissionReceipt` (node
  provenance hash + row schema v7 + color-algebra v2 +
  `color_admission_policy_fingerprint()`) ONLY for nodes that are known,
  unwaived (`scientific_color()` present — direct and transitive waiver
  taint both refuse), positively ranked (Verified/Validated; Estimated
  refuses), and replay-clean (`verify_replay_node` re-earns the provenance
  hash and stored state at mint time — a tampered node refuses with the
  exact `ColorReplayError`).
- `admission_receipt_in_regime(id, state)` additionally refuses a Validated
  node whose regime excludes the CURRENT execution state, returning the
  exact demotion the regime check derived: regime exit demotes structurally
  and never converts at the stale rank.
- `LedgerColorAdmissionVerifier` is the injected `AdmissionVerifier`
  capability: acceptance re-derives everything — receipt versions and
  policy fingerprint must match this build, the node hash must name a live
  node, the candidate must equal the node's scientific color in CANONICAL
  BYTES (display JSON is never trusted), and the node must still replay. A
  receipt minted before a tamper dies with the tampered graph.
- No-claim: authority is capability injection, not cryptography — a lying
  verifier at the composition root can admit anything, exactly like a lying
  `WaiverVerifier`. Receipts do not bind a graph instance identity; two
  graphs replaying identical writes mint interchangeable receipts (the node
  hash chain IS the identity). Consumer-API migration to require
  `AdmittedColor` is staged in bead 6pf9 S2-S4.

## No-claim boundaries

- Multi-process multi-writer access: unclaimed (FrankenSQLite documents this
  as partial; use one process, one connection per thread).
- Unbounded restart snapshots are not claimed: a physical ledger with more
  than 8,192 rows in either session-claim mirror, or more than 8,192 distinct
  authorities across their union, is refused before governor construction
  rather than allocating or scanning without a fixed limit.
- BLAKE3 keyed hashing, key derivation, XOF output beyond 32 bytes: not
  implemented.
- Cryptographic security claims: the implementation matches official vectors
  but has no side-channel or performance hardening (scalar, unoptimized).
- `LedgerInstanceId` is a collision-resistant uniqueness token, not a secret,
  signature, or authentication credential. Byte-for-byte database copies
  intentionally retain one identity because they are copies of one lineage.
- A checkpoint receipt proves that the referenced artifact and executor report
  are mutually bound and durably present. It does not independently prove that
  every solver thread was registered with the executor `DrainTracker`, that the
  snapshot contains sufficient application-specific resume state, or that its
  numerical contents are scientifically correct.
- A quantity crosswalk proves only the exact supported schema transformation:
  five exponents were preserved and an exact zero mole exponent was appended.
  It does not claim broader physical equivalence, provenance, or authority for
  either quantity artifact; it neither mints nor preserves a `SemanticType`
  or quantity-kind identity and cannot authorize dimensionally equal
  substitutions.
- A semantic state-checkpoint receipt proves exact persistence and equality of
  an injected semantic tuple. The nominal state-slot wrapper prevents implicit
  role exchange inside fs-ledger, but raw-hash adaptation and
  `KnownStateSemantics` remain caller assertions: fs-ledger does not prove that
  fs-ir, fs-matdb, or an L3 owner minted them. It also does not interpret the
  runtime-state codec, prove that the bytes are a sufficient restart snapshot,
  or prove cancel/resume/bit-replay. Those authorities remain with the
  executable L3 law and the E0d machine lifecycle. A generic v10 solver
  checkpoint is not silently upgraded into this stronger claim.
- Safe std-only identity generation is implemented through `/dev/urandom` on
  Unix. Fresh identity creation on non-Unix targets is explicitly refused;
  existing v4+ ledgers remain readable when their persisted identity and
  schema attest. Checked boundaries detect a guard while it is missing or
  different, and already-open handles detect row drift after the exact guard is
  restored. A client with arbitrary DDL authority can still rewrite an identity
  before a new handle establishes its cache, so the identity is not
  cryptographically authenticated against a hostile database owner.
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
- The v1 `ops` DDL checks JSON syntax and seed non-emptiness but does not encode
  the API's per-field 1 MiB ceilings. Canonical writes enforce them, `op`
  metadata-preflights and guards reads, and lint reports violations; arbitrary
  raw SQL is detected and refused rather than prevented as a DDL property.
- The v1 `edges` DDL does not itself freeze an operation's edge set when
  `outcome` becomes non-NULL. The public API is the sealing authority: `link`
  atomically selects only a canonical unfinished op, and no public edge
  update/delete API exists. A client with arbitrary SQL/DDL authority can still
  bypass that API boundary. Schema-v9 edge-set triggers provide a stronger raw
  SQL freeze only when a consumer explicitly creates an
  `op_artifact_edge_seals` row; `finish_op` does not infer or create that
  consumer-specific seal.
- The v1 `tune` DDL checks JSON syntax but does not encode the canonical
  identity, machine, JSON byte, row-count, or scan-byte bounds. The public API
  enforces them on writes, metadata-preflights and guards reads, and reports
  bounded envelope violations through lint; arbitrary raw SQL is therefore
  detected and refused, not prevented as a DDL property.
- Registry rows produced by the earlier uncommitted, unwired session-registry
  scaffold are not a compatibility claim. In particular, no dual verifier
  auto-trusts its unpublished v1 batch/event hash domains; such rows fail
  closed. Compatibility covers the tracked v6 table shape and the supported
  v2 writer preimages described above.
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
- `verify_color_row_stream` proves row-shape completeness and exact v9 hash
  reconstruction only. It deliberately does not decode canonical color bytes
  into scientific `Color` values, rerun typed-origin or waiver capabilities,
  decide present-day grant expiry, or confer admission authority. A complete
  persisted-row-to-`ColorGraph` importer remains unclaimed.
- Waiver expiry is checked at the authorizing node's admission day. Descendants
  preserve that historical grant and remain tainted indefinitely; they do not
  silently renew it, and `verify_replay()` has no caller-supplied current day
  with which to make a new policy decision.
- Transitive waiver visibility currently stores the complete unique grant
  closure on each descendant, bounded by `MAX_WAIVER_DEPENDENCIES` and
  `MAX_WAIVER_CLOSURE_BYTES`. This is deliberately inspectable and replayable;
  compact/sublinear waiver-lineage storage is not claimed.

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
