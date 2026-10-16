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
- `CanonicalEncoder`, `CanonicalLimits`, `SchemaDescriptorKind`, and
  `CanonicalError` — streaming, fail-closed construction under explicit frame,
  field, per-schema field-count, recursive schema-expansion, collection/chunk,
  and cancellation-poll budgets. `max_field_bytes` must admit the longest
  domain, name, context, or field name across the complete recursively bound
  schema descriptor. A one-under cap returns a distinct typed construction
  error naming the binding term, exact floor, and supplied cap; payload
  overruns remain `LimitExceeded { kind: FieldBytes, .. }`.
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
  consuming, explicitly policy-relative transitions driven by separately
  injected generic verifier and admission capabilities. These typestates are
  diagnostic/policy results and never promotion authority.
- `PromotionCapabilityDescriptor`, `OwnerPromotionVerifier`,
  `OwnerPromotionAdmitter`, `PromotionDecisionScope`,
  `PromotionDecisionRequest`, `PromotionAdmissionRequest`, and the four
  nominal request/verification/admission/committed decision IDs — the v3
  owner-executed decision protocol. Descriptors bind exact retained
  implementation/configuration observations and protocol versions, but are
  self-asserted metadata rather than proof of execution.
- `PromotionTrustRoot<V, P>` — a `Copy`, configuration-only compatibility
  profile. It can derive current and legacy charters but owns no executable
  capability and cannot mint a promotion witness.
- `PromotionTrustRoot<V, P, OwnerPromotionCapabilities<RV, RA>>` — the
  non-copyable, live owner root returned by `configure_owner_executed`. It owns
  non-replaceable verifier/admitter values, an epoch and a private instance seal;
  burns sequence state before calls; catches and poisons on capability panic;
  and is the only producer of `PromotionWitness` replay evidence.
- `PromotionWitness` and `OwnerBoundPromotion` — respectively raw,
  cloneable decision/replay evidence and the nominal authority obtained only
  by binding that evidence back to the exact live root and capability types.
  `crosswalk_witness` performs the sole in-process reconstruction: the target
  owner root re-executes both capabilities over an exact predecessor-bound
  scope and mints a new target-root witness. `PromotionAuditRecord` is the
  fixed-size, cloneable-but-non-`Copy`, payload-free replay receipt carrying
  every request axis, descriptor, statement root, approved disposition, and
  transcript ID needed by a separately versioned checker. The two materialized
  request types are likewise non-`Copy` and are exposed to capabilities only
  through immutable borrows, preventing accidental copies of the complete
  authority frame.
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
- Schema admission computes the global descriptor-string floor across the
  complete depth-capped child-binding graph in the same bounded,
  cancellation-polled traversal that validates its structure. Cancellation,
  malformed descriptors, and schema-expansion budget refusals remain
  authoritative during that traversal. For a structurally admitted schema, a
  `max_field_bytes` value equal to the floor admits; a value below it refuses
  before header or caller payload absorption through
  `MaxFieldBytesBelowSchemaFloor`. Once the descriptor is admitted,
  `LimitExceeded { kind: FieldBytes, .. }` is reserved for actual encoded
  field, row, or collection-item payload limits.
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
- Configuration alone cannot promote. The const three-argument
  `PromotionTrustRoot::configure` preserves a `Copy` policy profile and exact
  v1/v2 replay, but `admit_for_promotion` always refuses after optional binding
  diagnostics. Presenting the OWNER verifier/policy IDs and byte observations
  through a foreign generic `PermitAll` is therefore still only
  policy-relative admission.
- Root context is non-empty and capped at `MAX_PROMOTION_CONTEXT_BYTES` UTF-8
  bytes before any charter/request is derived, keeping the retained request and
  audit envelope bounded. The explicit legacy-v2 replay function retains v2's
  historical non-empty/depth validation and may reproduce a larger old context;
  its nominal result cannot enter current promotion APIs.
- Promotion-capable execution is exclusively
  `configure_owner_executed(...).decide_for_promotion(...)`. The resulting
  root physically stores distinct verifier and admission capability values;
  descriptors are read from those stored values, not passed as parallel
  expected IDs. The root re-adjudicates verifier/policy IDs and observations,
  requires an exact epoch and strictly increasing sequence, derives a complete
  request identity over subject role/schema/ID/preimage/length, external
  anchor, verifier/policy roles, complete schemas, IDs and byte observations,
  v3 charter, root context, attempt/context, epoch/sequence and optional
  crosswalk predecessor, then executes its stored verifier and admission
  capability in that order. Admission receives the immutable base request plus
  the exact stored verifier descriptor, returned verifier statement ID, and
  root-derived verification transcript. Each stage returns only a
  statement/reason artifact; the root alone derives stage transcripts and the
  final decision ID.
- Decision sequences are BURNED and execution is marked in-flight before the
  first capability call. Refusal and cancellation return to ready without
  restoring the sequence. A verifier/admitter panic is caught, publishes no
  witness, and permanently poisons that live root. No partially verified,
  partially admitted, cancelled or faulted decision sequence can be resumed or
  minted.
- `PromotionWitness` has private fields and no constructor, but private fields
  alone are not treated as authority. It is cloneable raw replay evidence and
  carries the complete request/stage/final decision identities, both exact
  capability descriptors, scope, subject/anchor/policy bindings, v3 charter,
  and a private `Arc` root-instance seal. A promotion consumer must obtain the
  nominal `OwnerBoundPromotion<..., RV, RA>` from the exact root's
  `bind_witness`; matching public charters/IDs/observations/descriptors from a
  foreign permit-all root refuse as `ForeignRootInstance`. An independently
  reconstructed equivalent root also refuses direct binding even when its v3
  charter matches. `audit()` retains the exact subject/anchor, roles, complete
  schemas, policy IDs/observations, descriptors, statement roots, scope and
  explicitly approved final disposition plus request/stage/final identities as
  bounded values; it retains no payload preimages and remains replay evidence
  rather than bearer authority.
- The only in-process reconstruction is `crosswalk_witness`: the source must
  already be live-root-bound, the new scope must name its exact committed
  decision as predecessor, and the target root re-executes BOTH of its owner
  capabilities before minting a new target-instance witness. A portable or
  cross-process crosswalk still needs an independently authenticated external
  anchor/signature/attestation and is not implemented by this leaf crate.
- Current `PromotionRootCharter` v3 binds the full prior v2 configuration plus
  capability mode, decision epoch, and both implementation/configuration
  observations and protocol versions. V3 deliberately excludes the private
  live-instance seal and mutable replay state, so equal charters are portable
  configuration equality, NOT bearer authority. Historical v2 is reconstructed
  only as nominal `legacy::PromotionRootCharterV2`; v1 remains nominal
  `PromotionRootCharterV1` with its historical same-domain/schema collapse.
  Neither legacy type has a parser, current conversion/equality bridge,
  strong-identity implementation or promotion acceptance path. Owner-local v3
  charter plus request, verification, admission, and committed-decision
  declarations are present; registry/golden regeneration and union proof remain
  central batch work.
- Same typed ID plus differing caller-supplied byte-root or length observations
  yields a refusal preserving both observations; it is not itself proof of a
  cryptographic collision.

## Error model

Compatibility parsers (`from_hex`, `from_slice`) return `Option` and refuse
malformed input with `None`. Canonical construction returns `CanonicalError`
for invalid limits or schemas, a schema-specific `max_field_bytes` value below
the recursive descriptor floor, arithmetic overflow, resource limits, field
mismatch or incompleteness, declared-length mismatch, nonfinite floats, invalid
set order or duplicates, and cancellation. The descriptor-floor refusal
retains `SchemaDescriptorKind`, the exact required floor, and the supplied cap;
it is never reported as a caller-payload overrun. No partial identity is
published on any of these paths. `ordered_bytes_stream` instead returns
`OrderedBytesStreamError<E>` so exact length-source and row-producer errors
remain distinct from canonical failures; both variants carry structured
fail-closed diagnostics and return no encoder. Crate-controlled inputs do not
panic; caller-provided iterators, row producers, and cancellation callbacks may
themselves panic.
Plain hashing inputs longer than 2^54 chunks are outside the supported
envelope (vastly beyond any in-tree artifact size).

Promotion configuration refuses empty context, over-depth schemas, zero
decision epoch, and zero capability-protocol versions. Structural binding
refuses foreign verifier/policy IDs and same-ID/different-byte observations.
Owner execution additionally refuses missing capabilities, the historical
unscoped path, wrong epochs, stale/replayed sequences, re-entrancy, capability
refusal/cancellation/panic, poisoned roots, foreign live-root instances,
non-recomputing witness transcripts, and wrong crosswalk predecessors. Exact
fixed-size capability reason/cancellation artifact IDs and the failing stage
survive the refusal; panic payloads do not. An unwinding panic is caught at this
leaf boundary with `catch_unwind`, but process aborts and foreign-code side
effects are outside the contract.

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

V3 root charters and promotion request/stage/final decision IDs are pure,
cross-ISA functions of their complete declared byte inputs. Repeating one
successful decision scope on a freshly reconstructed equivalent root with
deterministic capabilities yields the same portable decision IDs, but not the
same private live-root seal. Sequence state and the seal are deliberately not
portable hash inputs. Capability behavior itself inherits the determinism
contract of the stored implementation; this crate neither introspects code nor
upgrades a nondeterministic capability merely because its self-declared
descriptor is stable.

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

The promotion protocol has no hidden global cancellation token. A stored owner
capability reports its exact cancellation as
`PromotionCapabilityVerdict::Cancelled { reason }`; the root burns the
sequence before the call, returns a stage-specific typed refusal, publishes no
witness, and remains usable only at a higher sequence. If capability code
blocks without returning, this leaf cannot bound latency. Panic is distinct
from cancellation and permanently poisons the root.

This dependency-free leaf invokes each stored stage at most once per attempt,
but it does not yet admit or enforce a `Cx`, deadline, memory lease, or internal
work budget for code inside those calls. Adding a versioned external budget
authority and bounded consumption receipt is required before the parent Bead's
resource-accounting acceptance item can close; a self-reported descriptor is
not substituted for enforcement.

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
hostile bounds and invalid schemas; four-way plus recursive-child
schema-descriptor floor admission at one-under/exact/one-over boundaries;
typed separation of descriptor admission from payload `FieldBytes` refusal;
identity stability across admitted floors; collision refusal; authority
transitions and refusals; bounded audit records; typed parsing; legacy
quarantine; and compatibility behavior. Cancellation regression cases include
`cancellation_at_every_checkpoint_publishes_no_partial_identity` and
`cancellation_covers_schema_validation_and_long_set_comparisons`. Ordered-row
G0/G3/G4/G5 cases compare whole receipts and an independent manual frame;
exercise empty, chunked, reordered, and differently bounded schedules; retain
typed producer failures; refuse wrong counts, under/overwrites, hostile limits,
ignored sink failures, and arithmetic overflow; and cancel at every observed
row/chunk/publication checkpoint. These tests alone do not establish the
complete G4 tier.

`tests/authority_promotion.rs` is the v3 G0/G3/G4 authority battery. It proves
the exact OWNER-ID/OWNER-observation `PermitAll` bypass remains
policy-relative; configuration-only/unscoped downgrade refuses; owner execution
binds subject/anchor/context/policy/capability/epoch/scope; capability mode and
every descriptor axis move the v3 charter; same-charter foreign executable
types and swapped descriptors cannot bind to the owner live root; wrong
subject/anchor/context/epoch, same-ID/different-bytes, stale/replay and wrong
crosswalk predecessor refuse; cancellation/refusal burn before call; verifier
and admission panic publish no witness and poison; and equivalent
reconstruction succeeds only by a predecessor-bound target-owner crosswalk.
Independent buffered references retain exact v1/v2 replay plus both
configuration-only and owner-executed v3 grammars. Subject, anchor, attempt,
decision context, sequence, and returned stage statements independently move
the canonical request/stage/final decision identities. Central
all-target/doctest/clippy, identity registry regeneration, consumer migration
and DSR are still required before the Bead may close.

## No-claim boundaries

Digest equality is conditional on BLAKE3 collision resistance; it is not proof
of mathematical identity. Content, schema, and semantic IDs do not prove
origin, authenticity, external authority, or scientific correctness. Parsing
validates shape only. Even `Admitted` retains
`ScientificCorrectnessNotProven`: injected authority capabilities are policy
decisions, not built-in cryptographic verification.

Owner execution proves only that the exact capability VALUES stored in this
live safe-Rust root returned approval for the exact request/transcripts. It
does not prove those capabilities are scientifically meaningful, free of
interior mutability, deterministic, uncompromised, or the code described by
their self-asserted descriptors. Descriptor IDs/observations, v3 charter and
decision hashes authenticate nothing by themselves. The private `Arc` seal is
an in-process nominal capability boundary, not a cryptographic secret,
signature, hardware attestation, process identity, durable ledger root or
cross-machine bearer token. Raw `PromotionWitness` is replay evidence only;
only `OwnerBoundPromotion` is the in-process consumption authority.

The owner-executed crosswalk is in-process and target-policy-relative. It does
not authenticate a serialized witness after process restart. Portable/offline,
cross-process or cross-machine promotion requires an independently anchored
signature/attestation and revocation/replay store outside this crate. Sequence
tracking is bounded to one live mutable root and is not a durable global replay
ledger. The attempt ID is a transcript/correlation axis, not an independently
retained uniqueness set: a caller may deliberately retry the same attempt ID
at a higher sequence. `bind_witness` is likewise a reusable live-root proof,
not a one-shot consumption store; consumers that require single use must
durably adjudicate `decision_id`. `catch_unwind` contains Rust unwinding but
cannot contain aborts, signals, memory corruption, or side effects already
performed inside foreign capability code.

Only an approved two-stage completion publishes canonical witness/audit
evidence in this slice. Refusal, cancellation, and panic return typed bounded
diagnostics and burn replay state, but do not yet publish a portable canonical
outcome receipt binding budget consumption and final disposition. That durable
negative-outcome ledger is integration work, not an inferred success claim.
Statement/reason fields retain fixed-size `ContentId` roots, not referenced
payload lengths or bytes; payload-level retention needs its own independently
bounded observation record.

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
