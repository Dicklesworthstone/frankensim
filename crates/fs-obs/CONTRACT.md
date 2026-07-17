# CONTRACT: fs-obs

## Purpose and layer
Structured observability: the ONE event schema for kernels, solvers, test
suites, and (once fs-ledger lands) the ledger `events` table. Layer: UTIL.

## Public types and semantics
- `Event { session, scope, seq, severity, kind, wall_ns }` — envelope + typed
  payload; `to_jsonl()` canonical single-line display transport;
  `content_identity()` is the version-2 typed, exact-bit replay encoding and
  `content_hash()` is its FNV-1a root, EXCLUDING wall-clock.
- `EventIdentityReceipt { declared identity version, canonical bytes, root }`
  retains the complete event-identity proof. `content_identity_receipt()`
  captures it, `from_retained_parts()` reconstructs untrusted stored parts,
  and `Event::admit_content_identity()` accepts only a supported version, a
  root derived from the retained bytes, and bytes exactly equal to the
  admitted event's current canonical identity.
- `EventKind` v1 registry: solver_residual, tile_complete, cancellation,
  budget_delta, gradient_check, conformance_case, benchmark_result,
  storm_assertion, custom (pre-serialized JSON escape hatch whose UTF-8 is
  identity-bearing opaque bytes, not alleged canonical JSON).
- `Severity` (trace/info/warn/error), `Emitter` (per-scope monotone seq),
  `validate_line` (strict structural validator), `lint_failure_record`
  (failure-records-must-reproduce lint v1), `fnv1a64`, `SCHEMA_VERSION`.
- `ident::ReplayIdentity` and `ident::IdentityBuilder` — schema-v1 canonical
  replay identity and the original infallible builder for already-bounded
  internal producers. `ident::BoundedIdentityBuilder` emits exactly the same
  bytes and root while requiring an explicit canonical-byte ceiling. Its
  consuming typed appends return `IdentityBuildError`, reserve a complete
  field before mutation, and cannot finish after a refusal.
- `check_event_content_identity_version` refuses retained event identities
  from any version other than `EVENT_CONTENT_IDENTITY_VERSION`; owner-local
  declarations for both event content and the replay frame feed the generated
  `identity-schemas.json` policy gate.

## Invariants
- One event = one line; strings escaped so no literal newlines appear.
- Canonical field order: v, session, scope, seq, severity, kind, payload,
  [wall_ns last]. The golden-line test freezes this shape; changing it
  requires a SCHEMA_VERSION bump with semantic justification.
- `content_hash` is independent of `wall_ns` (deterministic-mode logs from
  two runs of the same seed hash identically). It never hashes `to_jsonl`:
  exact float bits, including distinct NaN payloads, remain semantic even when
  the human-readable JSON tags both values as `non-finite:NaN`.
- Retained event identities are never accepted from a naked root. Admission
  binds the declared `EVENT_CONTENT_IDENTITY_VERSION`, the full canonical
  bytes (which independently frame the artifact domain, that identity version,
  and `SCHEMA_VERSION`), the FNV root of those bytes, and the exact event being
  admitted. `from_retained_parts()` is deliberately a raw transport
  constructor; its declared version and bytes remain untrusted semantic input
  until `Event::admit_content_identity()` validates them.
- `Custom::json` is exact opaque UTF-8 identity material. Whitespace and
  object-member order therefore move the root. fs-obs does not invoke or claim
  an unchecked JSON canonicalizer; the identity battery decodes the typed
  byte field and proves it round-trips exactly.
- Additive schema evolution only: kinds may be added, fields never repurposed.
- Non-finite floats serialize as tagged strings ("non-finite:NaN"), never
  invalid JSON.
- Replay identity v1 begins with the exact bytes of
  `ident::REPLAY_IDENTITY_DOMAIN`, then frames the schema version, kind, and
  each typed field with little-endian u64 lengths. There is no separate magic
  constant that can drift from the declared domain. The bounded builder checks
  native-to-u64 framing, checked total length, and the producer's byte cap
  before reserving and appending. Field tags, order, duplicate keys, float
  bits, and child identities are semantic; documented exclusions remain
  outside the canonical byte stream.

## Error model
`SchemaError { at, message }`, `EventIdentityVersionError`, and
`EventIdentityAdmissionError` (unsupported semantics, retained root/preimage
mismatch, or exact-event mismatch) provide structured refusal. `IdentityBuildError`
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
The unit and integration batteries cover all nine kinds' serialize+validate
round-trip; an ordered count-locked mutation of all 26 payload fields; golden
line; envelope/content hash split; every top-level semantic event field;
independent event-identity-version and wire-schema bytes/roots;
same-display/different-bit NaNs; exact opaque Custom bytes; retained receipt
admission and each refusal class; stale identity versions; monotone sequences;
corruption rejection; failure-record lint; hostile-string escaping;
non-finite tagging; FNV known answers; and exhaustive replay-frame mutations,
including direct domain-prefix binding.

## No-claim boundaries
- FNV-1a is NOT cryptographic; ledger-grade content addressing (BLAKE3-class)
  arrives with fs-ledger.
- Scope-tree mirroring is IN-CRATE (huq.16): `Emitter::enter_scope` /
  `exit_scope` mirror the asupersync scope tree explicitly (the scope tree IS
  the trace tree), refusing path-forging segments (`/`, control characters)
  and unbalanced exits, with one monotone `seq` stream so interleaved child
  scopes replay in exact emission order. Runtime layers walk their live scope
  trees through these calls; this crate stays deterministic and I/O-free.
- The ledger `events` sink surface is the documented pair
  (`EventKind::kind_name`, `Event::to_jsonl`) — sinks store both without
  re-encoding. The sink implementation itself lives with fs-ledger's owners;
  committed `*.events.jsonl` fixtures are schema-enforced in CI by
  `xtask check-obs-events`, whose validator authority is
  `fs_obs::validate_line` (never a second dialect). Overhead budgeting
  (roofline harness) remains unclaimed bead scope.
- The validator is structural, not a full JSON parser (the writer is ours;
  external JSON is out of scope). The JSON line is presentation/transport, not
  the exact-bit content-identity preimage. Public `Custom` construction still
  requires the caller to supply one valid JSON object and run `validate_line`;
  its identity attests only the exact opaque UTF-8, not validity or semantic
  JSON equivalence.
- The legacy `IdentityBuilder`, `ReplayIdentity::clone`, and `hex()` remain
  allocation-infallible conveniences. Public admission paths that can receive
  resource-driving input must use `BoundedIdentityBuilder`; later migrations
  must make their owner APIs fallible rather than wrapping it in `expect`.
