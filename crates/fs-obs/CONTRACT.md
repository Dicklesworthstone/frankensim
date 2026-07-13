# CONTRACT: fs-obs

## Purpose and layer
Structured observability: the ONE event schema for kernels, solvers, test
suites, and (once fs-ledger lands) the ledger `events` table. Layer: UTIL.

## Public types and semantics
- `Event { session, scope, seq, severity, kind, wall_ns }` — envelope + typed
  payload; `to_jsonl()` canonical single-line serialization; `content_hash()`
  (FNV-1a 64) over content EXCLUDING wall-clock.
- `EventKind` v1 registry: solver_residual, tile_complete, cancellation,
  budget_delta, gradient_check, conformance_case, benchmark_result,
  storm_assertion, custom (pre-serialized JSON escape hatch).
- `Severity` (trace/info/warn/error), `Emitter` (per-scope monotone seq),
  `validate_line` (strict structural validator), `lint_failure_record`
  (failure-records-must-reproduce lint v1), `fnv1a64`, `SCHEMA_VERSION`.
- `ident::ReplayIdentity` and `ident::IdentityBuilder` — schema-v1 canonical
  replay identity and the original infallible builder for already-bounded
  internal producers. `ident::BoundedIdentityBuilder` emits exactly the same
  bytes and root while requiring an explicit canonical-byte ceiling. Its
  consuming typed appends return `IdentityBuildError`, reserve a complete
  field before mutation, and cannot finish after a refusal.

## Invariants
- One event = one line; strings escaped so no literal newlines appear.
- Canonical field order: v, session, scope, seq, severity, kind, payload,
  [wall_ns last]. The golden-line test freezes this shape; changing it
  requires a SCHEMA_VERSION bump with semantic justification.
- `content_hash` is independent of `wall_ns` (deterministic-mode logs from
  two runs of the same seed hash identically).
- Additive schema evolution only: kinds may be added, fields never repurposed.
- Non-finite floats serialize as tagged strings ("non-finite:NaN"), never
  invalid JSON.
- Replay identity v1 frames the kind and each typed field with little-endian
  u64 lengths. The bounded builder checks native-to-u64 framing, checked total
  length, and the producer's byte cap before reserving and appending. Field
  tags, order, duplicate keys, float bits, and child identities are semantic;
  documented exclusions remain outside the canonical byte stream.

## Error model
`SchemaError { at, message }` with fix guidance. `IdentityBuildError`
distinguishes canonical-byte-cap refusal, unrepresentable framing, length
overflow, and allocation failure. No panics cross the bounded identity or
event-validation boundaries; they reject rather than repair.

## Determinism class
Deterministic: pure functions; no clocks (callers supply `wall_ns`), no I/O,
no RNG.

## Cancellation behavior
All operations O(event size). No Cx required.

## Unsafe boundary
None.

## Feature flags
None.

## Conformance tests
Unit suite (11 cases): every-kind serialize+validate round-trip, golden line,
envelope/content hash split, monotone sequences, corruption rejection,
failure-record lint, hostile-string escaping, non-finite tagging, FNV known
answers.

## No-claim boundaries
- FNV-1a is NOT cryptographic; ledger-grade content addressing (BLAKE3-class)
  arrives with fs-ledger.
- Tracing-on-scope-trees integration (needs fs-exec), ledger `events` sink
  (needs fs-ledger), and overhead budgeting (needs roofline harness) are this
  crate's REMAINING bead scope — not yet claimed.
- The validator is structural, not a full JSON parser (the writer is ours;
  external JSON is out of scope).
- The legacy `IdentityBuilder`, `ReplayIdentity::clone`, and `hex()` remain
  allocation-infallible conveniences. Public admission paths that can receive
  resource-driving input must use `BoundedIdentityBuilder`; later migrations
  must make their owner APIs fallible rather than wrapping it in `expect`.
