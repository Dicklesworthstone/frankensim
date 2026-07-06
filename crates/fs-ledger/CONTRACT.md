# CONTRACT: fs-ledger

> Status: ACTIVE (Design Ledger v0). Owns schema v0 + Rev S extension tables,
> BLAKE3 content addressing, and the WAL/snapshot concurrency contract.
> Time travel, forkable worlds, and `explain()` belong to the follow-on
> fs-ledger time-travel bead and are NOT provided here.

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
  32-byte output), pure safe Rust; artifact identity everywhere.
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
  invariants, half-finished ops) — all-zero on any healthy or crash-recovered
  ledger.

Schema divergences from plan Appendix D, both deliberate: `JSON` columns are
STRICT-legal `TEXT` with `json_valid()` CHECKs (Appendix D as written is not
valid STRICT SQL), and `artifacts` gains `len`/`chunk_count` +
`artifact_chunks` for bounded-memory large-field storage.

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
throughput smoke ledgered as a metric. Unit tests in `src/lib.rs` and
`src/hash.rs` cover the API surface and edge cases.

## No-claim boundaries

- Multi-process multi-writer access: unclaimed (FrankenSQLite documents this
  as partial; use one process, one connection per thread).
- BLAKE3 keyed hashing, key derivation, XOF output beyond 32 bytes: not
  implemented.
- Cryptographic security claims: the implementation matches official vectors
  but has no side-channel or performance hardening (scalar, unoptimized).
- Throughput numbers are smoke floors, not roofline claims (§14 discipline:
  real claims need machine fingerprints and acceptance bands).
- Time travel, forks, `explain()`: follow-on bead.
- Multi-GiB single artifacts: chunk storage bounds row sizes, but the
  streaming path is verified at the tens-of-MiB scale only so far; fsqlite
  transaction memory behavior at multi-GiB scale is unmeasured.
