# fs-blake3 CONTRACT

## Purpose and layer

Layer: **UTIL**. This crate is the single in-tree owner of both the BLAKE3
hash function (plain and derive-key modes, fixed 32-byte output) and the
schema-typed canonical identity substrate. The hash implementation was
extracted verbatim from `fs-ledger` (bead 7uq9), which re-exports the
compatibility types unchanged. Zero dependencies and safe Rust are deliberate:
solver-free distribution cones such as `fs-checker` may depend on this crate
without gaining solver, geometry, FFI, or license surface.

## Public types and semantics

- `Blake3` — streaming hasher. `new()`, `update(&[u8])` (any split
  pattern; digest equals hashing the concatenation), `finalize() ->
  ContentHash`. `finalize` borrows immutably: it may be called
  repeatedly and interleaved with further `update`s.
- `hash_bytes(&[u8]) -> ContentHash` — one-shot convenience.
- `hash_domain(domain: &str, payload: &[u8]) -> ContentHash` — the
  canonical domain-separation scheme for every typed 32-byte root in the
  workspace. It uses BLAKE3's standard derive-key construction:
  `DERIVE_KEY_CONTEXT` for the domain and `DERIVE_KEY_MATERIAL` for the
  payload. Mode flags separate typed roots from plain `hash_bytes`
  artifacts as well as from other domains.
- `ContentHash(pub [u8; 32])` — Copy identity value. `as_bytes`,
  `to_hex` (64 lowercase hex chars), `from_hex` (either case, exactly
  64 chars, else `None`), `from_slice` (exactly 32 bytes, else
  `None`). `Display` renders `to_hex`; `Debug` wraps it.
- `ContentId` — plain BLAKE3 identity of exact retained bytes. It carries no
  semantic, origin, or authority claim.
- `CanonicalSchema`, `FieldSpec`, `Presence`, `WireType`, and `Field` — a
  static, nominal schema descriptor and the exact field handles admitted by an
  encoder. `SchemaId<D>::for_schema()` names the declared descriptor; it does
  not by itself admit that descriptor against structure or resource limits.
- `IdentityRole` and the sealed strong-identity families (`SemanticId`,
  `WireContentId`, `EvidenceNodeId`, `EntityId`, `SourceByteId`, `SourceId`,
  `ModelId`, `CheckerId`, `VerifierId`, `KeyPolicyId`, and
  `ProblemSemanticId`), together with `SchemaId` — role- and schema-specific
  32-byte values. Nominal roles and schema marker types are
  non-interchangeable. Strict parsing checks byte shape only and adds no trust.
  In v1, explicit optional encoding is exposed only for `Bytes` through
  `optional_bytes`.
- `CanonicalEncoder`, `CanonicalLimits`, and `CanonicalError` — streaming,
  fail-closed construction under explicit frame, field, field-count,
  collection/chunk, and cancellation-poll budgets. Defaults are a 1 MiB frame,
  256 KiB field, 256 fields, 16,384 items per collection or chunks per streamed
  byte field, and a 4,096-byte cancellation poll stride.
- `IdentityReceipt` and `IdentityAuditRecord` — the successfully published
  typed identity, exact canonical-frame root, bounded counters/limits, and
  fixed-size audit metadata. Payload bytes are not retained in the record.
- `AuthorityRef` and the `Presented`, `Verified`, and `Admitted` typestates —
  consuming authority transitions driven by separately injected verifier and
  admission capabilities.
- `ByteObservation`, `ObservedIdentity`, `SameIdDifferentBytes`, and
  `adjudicate` — retained-evidence comparison for refusing one typed ID that is
  observed with different byte roots or lengths.

## Invariants

- Output is bit-identical to the official BLAKE3 specification for
  plain-mode hashing (spec/oracle vectors in the unit tests and in
  fs-ledger's `ledger_001` conformance suite).
- Streaming and one-shot hashing of the same byte sequence produce the
  same digest for every split pattern.
- `hash_domain(d1, p1) == hash_domain(d2, p2)` implies (up to a BLAKE3
  collision) `d1 == d2 && p1 == p2`; equality between a tagged hash and a
  plain `hash_bytes` artifact likewise requires a cross-mode BLAKE3 collision
  because the compression flags differ.
- `from_hex(to_hex(h)) == Some(h)` for all `h`.
- Canonical frame v1 binds its magic and version, identity role, exact-finite
  float policy, schema domain/name/ID/version/context, complete declared field
  schema and order, canonical field stream, and final field count.
- Required field ordinal, name, wire type, and presence must match exactly.
  `None` differs from empty bytes, and variants bind both numeric tag and
  payload.
- Ordered collections preserve caller order. Canonical sets must be strictly
  byte-lexicographic and duplicate-free. Typed children bind their role,
  complete schema descriptor identity, and all 32 digest bytes.
- Finite `f64` values are encoded by exact IEEE-754 bits. Signed zero remains
  distinct. NaN and infinities refuse before publication.
- Successful construction publishes both a derive-key typed ID and a plain
  BLAKE3 root of the complete canonical frame. Resource budgets are retained
  receipt metadata. Neither budgets nor a successful cancellation schedule are
  hash inputs; the schedule itself is not retained.
- Fallible encoder operations consume the encoder. A refused or cancelled
  construction cannot be resumed or finished and publishes no receipt or root.
- Authority advances only `Presented -> Verified -> Admitted`. Anchor,
  verifier, or policy presence alone is untrusted; verification and admission
  are separate consuming trait calls over the exact stored `AuthorityRef`. The
  same concrete capability may implement both interfaces; correctness of the
  injected checks remains caller/policy responsibility.
- Same typed ID plus differing caller-supplied byte-root or length observations
  yields a refusal preserving both observations; it is not itself proof of a
  cryptographic collision.

## Error model

Compatibility parsers (`from_hex`, `from_slice`) return `Option` and refuse
malformed input with `None`. Canonical construction returns `CanonicalError`
for invalid limits or schemas, arithmetic overflow, resource limits, field
mismatch or incompleteness, declared-length mismatch, nonfinite floats,
invalid set order or duplicates, and cancellation. No partial identity is
published on any of these paths. Crate-controlled inputs do not panic;
caller-provided iterators and cancellation callbacks may themselves panic.
Plain hashing inputs longer than 2^54 chunks are outside the supported
envelope (vastly beyond any in-tree artifact size).

## Determinism class

The typed ID and canonical-frame root are pure functions of the canonical input
bytes — bit-stable across runs, thread counts, platforms, and ISAs. There is no
floating-point arithmetic: finite values are represented by exact bits.
Successful stream partitioning, larger admissible limits, and non-cancelling
probe schedules do not move identity. Receipt metadata may differ when
admissible limits differ. Time, host state, locale, and I/O are not hash inputs.

## Cancellation behavior

`Blake3`, `hash_bytes`, and `hash_domain` remain synchronous and
non-cancellable. Canonical construction takes an explicit `CancellationProbe`;
`NeverCancel` is the deliberate opt-out. Encoding polls during initialization,
schema validation and descriptor emission, header sizing and absorption,
streamed chunks, long comparisons, payload absorption, and immediately before
publication. `cancellation_poll_bytes` must be positive. Cancellation returns
`Cancelled { absorbed_bytes }`, consumes the encoder, and publishes neither
root. No latency claim applies while caller iterator or probe code is blocked.

## Unsafe boundary

None. 100% safe Rust; no capsule.

## Feature flags

None.

## Conformance tests

Unit tests in `src/lib.rs`: official empty-input spec vector, oracle
`abc` vector, hex round-trip and rejection, streaming-vs-one-shot
across block/chunk edges, single- and multi-chunk derive-key oracle
vectors, and raw/tagged namespace separation. The historical
multi-chunk / multi-level tree vectors continue to run in fs-ledger's
`ledger_001` conformance suite through the re-exported paths.

`tests/identity.rs` covers independent manual frame/schema parity and mutation
sensitivity; explicit tags and roles; field, schema, role, and child
non-confusability; exact-finite float policy; ordered/set/stream encodings;
hostile bounds and invalid schemas; collision refusal; authority transitions
and refusals; bounded audit records; typed parsing; legacy quarantine; and
compatibility behavior. Cancellation regression cases include
`cancellation_at_every_checkpoint_publishes_no_partial_identity` and
`cancellation_covers_schema_validation_and_long_set_comparisons`. These tests
alone do not establish the complete G4 tier.

## No-claim boundaries

Digest equality is conditional on BLAKE3 collision resistance; it is not proof
of mathematical identity. Content, schema, and semantic IDs do not prove
origin, authenticity, external authority, or scientific correctness. Parsing
validates shape only. Even `Admitted` retains
`ScientificCorrectnessNotProven`: injected authority capabilities are policy
decisions, not built-in cryptographic verification.

No Unicode, unit, JSON, locale, or domain-specific normalization is supplied.
Collision adjudication compares caller-supplied retained observations; it
establishes neither their independence nor authenticity and cannot detect
discarded distinctions or prove collision absence. `IdentityAuditRecord` is
bounded descriptive metadata, not signed or authenticated evidence. Legacy
FNV provenance remains quarantined with no conversion or equality bridge.

General-purpose keyed hashing, a public KDF API, and extended (XOF) output are
NOT implemented. `hash_domain` uses the standard derive-key modes only for
public, non-secret identity namespaces. No constant-time claim is made
(content addressing, not secret handling). No SIMD or multithreaded throughput
claim is made.
- Domain strings are caller-owned protocol identifiers. Callers must use
  hardcoded, globally unique, versioned contexts; this crate does not register
  or deduplicate names dynamically.
