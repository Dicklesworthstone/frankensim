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
- `DomainHasher` — the sealed streaming form of `hash_domain`. Callers choose
  a public protocol domain and feed payload chunks, but cannot select raw mode
  flags or secret-keyed hashing. Its final root is byte-identical to hashing
  the concatenated payload through `hash_domain`.
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
  fail-closed construction under explicit frame, field, per-schema field-count,
  recursive schema-expansion, collection/chunk, and cancellation-poll budgets.
  Defaults are a 1 MiB frame, 256 KiB field, 256 fields per descriptor, 16,384
  items per collection or chunks per streamed byte field, and a 4,096-byte
  cancellation poll stride. `CanonicalEncoder` admission of a recursive
  descriptor is additionally capped at
  `max_fields * (MAX_SCHEMA_CHILD_DEPTH + 1)` field occurrences.
- `CanonicalRowSink`, `OrderedBytesStreamError`,
  `OrderedBytesStreamDiagnostic`, `OrderedBytesStreamPhase`, and
  `OrderedBytesStreamDisposition` — allocation-free, fallible production of
  `OrderedBytes` rows. A separate fallible length source declares each row
  before its higher-ranked, encoder-owned sink absorbs borrowed chunks.
  Diagnostics retain schema/field/row progress, canonical and row byte
  counters, prior collection items, chunk-budget units admitted through the
  first refusal, refusal phase, and the consumed-without-publication
  disposition. Exact caller errors remain typed and distinct from
  `CanonicalError`.
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
- Plain and domain-separated streaming each match their one-shot form for the
  same byte sequence under every split pattern; empty chunks are nonsemantic.
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
  payload. Canonical-frame v1 exposes optional presence only through
  `FieldSpec::optional_bytes`, matching its one explicit optional encoder;
  unsupported optional schema declarations are unrepresentable instead of
  becoming admitted but unencodable schemas.
- Ordered collections preserve caller order. Canonical sets must be strictly
  byte-lexicographic and duplicate-free. Each set item's byte and aggregate
  field budgets are admitted before potentially long ordering comparisons.
  `ordered_bytes_stream` emits the same v1 grammar as `ordered_bytes`: one row
  count followed by exactly one little-endian length and the concatenated
  chunks for each row. Row and chunk scheduling metadata never enters the
  frame. The declared count, successful row-declaration count, and
  completed-row count agree exactly; receipt `collection_items` increases by
  exactly that count. Each row's individual and aggregate field bytes and
  complete-frame addition are admitted before its length prefix or producer
  callback. An underwrite, pre-absorption
  overwrite, fallible length-source or producer error, arithmetic overflow,
  resource refusal, or cancellation consumes the encoder. The first sink
  failure is sticky even when producer code ignores the immediate `write`
  result. A higher-ranked row callback and private sink fields prevent safe
  code from retaining the mutable encoder borrow beyond that row.
  Child fields are PARENT-BOUND (bead sj31i.52.10): a `FieldSpec` for
  `Child`/`OrderedChildren` must declare its expected child role and
  complete schema identity via `FieldSpec::child_of` /
  `ordered_children_of` + `ChildSpec::for_identity::<J>()` — an unbound
  child field is a compile-time error, and the encoder refuses any
  identity whose role/domain/name/version/context/field-schema differs from the
  declaration (empty ordered collections included). Field-schema comparison is
  structural and depth-capped; pointer identity is only a fast path/cycle guard
  because associated constants have no stable address. Distinct marker types
  with identical roles and complete descriptors are therefore
  admission-equivalent. `ChildSpec` equality/hash instead retain pointer-tail
  identity to stay total on recursive values; those traits are not the encoder's
  admission relation. The binding is folded RECURSIVELY (depth-capped with a
  deterministic poison tag) into the schema-descriptor preimage under
  `FSSCHEM\x02`, so changing the expected role or descriptor changes the parent
  `SchemaId`. `SchemaId::for_schema` can name an over-depth descriptor through
  that poison tag. `CanonicalEncoder` recursively validates every child
  descriptor, enforces the
  per-descriptor field cap plus a derived aggregate expansion ceiling, polls
  cancellation throughout, and refuses invalid or over-depth bindings before
  publishing an identity.
  Directly deriving the same descriptor under `FSSCHEM\x01` and `FSSCHEM\x02`
  therefore produces DIFFERENT roots. Because `SCHEMA_ID_HASH_DOMAIN`, the
  nominal types, and shape-only parsers remain `.v1`, the crate cannot
  distinguish an externally parsed old root from a current root. No cross-era
  authority or automatic migration is claimed; callers must quarantine old
  roots and perform an explicit external migration until the public domains and
  types rotate together.
  Typed children bind their role, complete schema descriptor identity, and all
  32 digest bytes.
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
  injected checks remains caller/policy responsibility, and the resulting
  `Admitted` state is therefore EXPLICITLY POLICY-RELATIVE
  (`PolicyRelativeAdmitted`), never promotion authority (bead sj31i.52.9).
- Promotion-capable admission is exclusively
  `PromotionTrustRoot::admit_for_promotion`: the domain owner independently
  configures exact verifier and key-policy identities WITH their
  canonical-byte observations; admission re-adjudicates the presented binding
  against those retained observations (foreign identities and
  same-ID/different-bytes presentations are typed refusals retaining both
  observations) and mints the `PromotionWitness` — private fields, no public
  constructor, sealed typestates (`compile_fail`-proven), so a foreign
  permit-everything capability can never fabricate one. The witness binds the
  exact subject receipt/preimage, anchor, verifier/policy observations, and
  root context. `PromotionTrustRoot::charter` derives the current v2 charter
  without allocation: explicit verifier/key-policy role tags, both schema
  domains, both COMPLETE recursively child-bound `SchemaId` values, both
  retained IDs, both observation roots and exact lengths, and the context are
  uniformly u64-length-framed and streamed through `DomainHasher`. A
  compile-time target-width guard plus checked conversion prevents host-length
  truncation while preserving the existing const `configure` API. Configure
  refuses verifier or key-policy schemas beyond the complete-descriptor depth
  before a root exists, so a poison-tagged over-depth `SchemaId` can never back
  charter authority or mint a witness. Every witness carries this current
  charter. Historical v1 charters can be
  reconstructed only as nominal `legacy::PromotionRootCharterV1` replay
  evidence. The dedicated legacy reconstruction function bypasses only the
  current schema-depth admission guard, preserving historical over-depth v1
  configurations without minting a current root. The replay type has no parser,
  conversion/equality bridge, strong-identity implementation, or acceptance
  path into current promotion authority. V1's same-domain/different-schema
  collapse is retained for faithful replay, never migration. `audit()` yields
  bounded namespace + observation-root/length
  metadata only. No-claim: the root's guarantees are relative to its
  configuration authority — fs-blake3 cannot vouch that the configured
  verifier is meaningful. The owner-local v2 identity declaration is present;
  registry/golden regeneration and union proof remain central batch work.
- Same typed ID plus differing caller-supplied byte-root or length observations
  yields a refusal preserving both observations; it is not itself proof of a
  cryptographic collision.

## Error model

Compatibility parsers (`from_hex`, `from_slice`) return `Option` and refuse
malformed input with `None`. Canonical construction returns `CanonicalError`
for invalid limits or schemas, arithmetic overflow, resource limits, field
mismatch or incompleteness, declared-length mismatch, nonfinite floats,
invalid set order or duplicates, and cancellation. No partial identity is
published on any of these paths. `ordered_bytes_stream` instead returns
`OrderedBytesStreamError<E>` so exact length-source and row-producer errors
remain distinct from canonical failures; both variants carry structured
fail-closed diagnostics and return no encoder. Crate-controlled inputs do not
panic; caller-provided iterators, row producers, and cancellation callbacks may
themselves panic.
Plain hashing inputs longer than 2^54 chunks are outside the supported
envelope (vastly beyond any in-tree artifact size).

## Determinism class

The typed ID and canonical-frame root are pure functions of the canonical input
bytes — bit-stable across runs, thread counts, platforms, and ISAs. There is no
floating-point arithmetic: finite values are represented by exact bits.
Successful stream partitioning, larger admissible limits, and non-cancelling
probe schedules do not move identity. Receipt metadata may differ when
admissible limits differ. Time, host state, locale, and I/O are not hash inputs.
For ordered-row streaming, changing chunk boundaries, empty-write placement,
or a successful callback schedule is nonsemantic. Changing row boundaries,
row order, declared lengths, or concatenated row bytes is semantic.

## Cancellation behavior

`Blake3`, `DomainHasher`, `hash_bytes`, and `hash_domain` remain synchronous and
non-cancellable. Canonical construction takes an explicit `CancellationProbe`;
`NeverCancel` is the deliberate opt-out. Encoding polls during initialization,
schema validation and descriptor emission, header sizing and absorption,
streamed chunks, long comparisons, payload absorption, and immediately before
publication. Ordered-row streaming additionally polls before each fallible
length declaration, on every sink write (including empty writes), after every
complete row, after exact source exhaustion, and through long chunk absorption.
`cancellation_poll_bytes` must be positive. Cancellation returns `Cancelled {
absorbed_bytes }`, consumes the encoder, and publishes neither root. No latency
claim applies while caller length-source, row-producer, iterator, or probe code
is blocked between crate-controlled checkpoints.

## Unsafe boundary

None. 100% safe Rust; no capsule.

## Feature flags

None.

## Conformance tests

Unit tests in `src/lib.rs`: official empty-input spec vector, oracle
`abc` vector, hex round-trip and rejection, plain and domain-separated
streaming-vs-one-shot across block/chunk edges, single- and multi-chunk
derive-key oracle vectors, and raw/tagged namespace separation. The historical
multi-chunk / multi-level tree vectors continue to run in fs-ledger's
`ledger_001` conformance suite through the re-exported paths.

`tests/identity.rs` covers independent manual frame/schema parity and mutation
sensitivity; explicit tags and roles; field, schema, role, and child
non-confusability; exact-finite float policy; ordered/set/stream encodings;
hostile bounds and invalid schemas; collision refusal; authority transitions
and refusals; bounded audit records; typed parsing; legacy quarantine; and
compatibility behavior. Cancellation regression cases include
`cancellation_at_every_checkpoint_publishes_no_partial_identity` and
`cancellation_covers_schema_validation_and_long_set_comparisons`. Ordered-row
G0/G3/G4/G5 cases compare whole receipts and an independent manual frame;
exercise empty, chunked, reordered, and differently bounded schedules; retain
typed producer failures; refuse wrong counts, under/overwrites, hostile limits,
ignored sink failures, and arithmetic overflow; and cancel at every observed
row/chunk/publication checkpoint. These tests alone do not establish the
complete G4 tier.

## No-claim boundaries

Digest equality is conditional on BLAKE3 collision resistance; it is not proof
of mathematical identity. Content, schema, and semantic IDs do not prove
origin, authenticity, external authority, or scientific correctness. Parsing
validates shape only. Even `Admitted` retains
`ScientificCorrectnessNotProven`: injected authority capabilities are policy
decisions, not built-in cryptographic verification.
`ordered_bytes_stream` bounds only encoder-owned counters and hash state. It
does not bound memory retained by a caller's length source, indexing scheme, or
row producer, and it requires lengths to be independently obtainable by row
index without first retaining every row payload. Consumer-side sorting,
workspace leases, graph indexes, and producer progress/resume semantics remain
outside this leaf's claim.

No Unicode, unit, JSON, locale, or domain-specific normalization is supplied.
Collision adjudication compares caller-supplied retained observations; it
establishes neither their independence nor authenticity and cannot detect
discarded distinctions or prove collision absence. `IdentityAuditRecord` is
bounded descriptive metadata, not signed or authenticated evidence. Legacy
FNV provenance remains quarantined with no conversion or equality bridge.
The internal `FSSCHEM\x01` to `FSSCHEM\x02` marker transition did not rotate
the public schema-ID or dependent canonical-frame domains. Shape-only parsing
cannot identify which marker produced an externally supplied digest; no
cross-era authority, automatic migration, or completed public-version
crosswalk is claimed until those domains and nominal types are versioned
together.

General-purpose keyed hashing, a public KDF API, and extended (XOF) output are
NOT implemented. `hash_domain` uses the standard derive-key modes only for
public, non-secret identity namespaces. No constant-time claim is made
(content addressing, not secret handling). No SIMD or multithreaded throughput
claim is made.
- Domain strings are caller-owned protocol identifiers. Callers must use
  hardcoded, globally unique, versioned contexts; this crate does not register
  or deduplicate names dynamically.
