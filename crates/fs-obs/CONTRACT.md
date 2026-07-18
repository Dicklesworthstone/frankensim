# CONTRACT: fs-obs

## Purpose and layer
Structured observability: the ONE event schema for kernels, solvers, test
suites, and (once fs-ledger lands) the ledger `events` table. Layer: UTIL.

## Public types and semantics
- `Event { session, scope, seq, severity, kind, wall_ns }` — envelope + typed
  payload; `to_jsonl()` canonical single-line display transport;
  `content_identity()` is the version-9 typed, exact-bit replay encoding and
  `content_hash()` is its FNV-1a root, EXCLUDING wall-clock.
- `EventIdentityReceipt { declared identity version, canonical bytes, root }`
  retains the complete event-identity proof. `content_identity_receipt()`
  captures it, `from_retained_parts()` reconstructs untrusted stored parts,
  and `Event::admit_content_identity()` accepts only a supported version, a
  root derived from the retained bytes, and bytes exactly equal to the
  admitted event's current canonical identity.
- `EventKind` v1 registry: solver_residual, tile_complete, cancellation,
  budget_delta, gradient_check, conformance_case, benchmark_result,
  storm_assertion, manifest_selection, stratum_expansion, dsr_run,
  campaign_run, submission_decision, work_identity_binding, lease_cursor,
  attach_detach, journey_phase, scope_tile_progress, heartbeat,
  observation_gap, oracle_comparison, tolerance_derivation,
  claim_adjudication, capability_domain_decision, lifecycle_transition,
  artifact_lifecycle, visualization_transform, diagnostic_repair,
  containment_node, containment_gap, race_record, degradation_event,
  import_receipt, certificate_verdict, custom (pre-serialized JSON escape
  hatch whose UTF-8 is identity-bearing opaque bytes, not alleged canonical
  JSON).
- `Severity` (trace/info/warn/error), opaque `SamplingCadence`
  (`never` or `every(NonZeroU64)`), and opaque `EmissionGate` (independent
  trace/info cadences; warn/error always admitted). `Emitter` owns the
  per-scope monotone emitted-event sequence; `emit_gated` takes an explicit
  logical opportunity ordinal and a lazy builder for both payload and optional
  caller clock. `validate_line` is the strict structural validator;
  `lint_failure_record` is failure-records-must-reproduce lint v1.
- `fnv1a64`, `SCHEMA_VERSION`.
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
- `process::{ProcessCapture, ProcessCapturePolicy, ProcessFrame, ProcessGap}`
  — deterministic, I/O-free process-stream admission with critical,
  diagnostic, and telemetry loss classes; bounded frame/gap queues;
  cancellation-aware backpressure; committed-artifact spill pointers; exact
  drop/range/policy accounting; canonical typed `process_frame`/`process_gap`
  event projection; and final-receipt reconciliation.
- `privacy::{LabeledField, FieldPolicy, ShareRequest, ShareManifest}` —
  field-level sensitivity, license, export, and retention policy; bounded
  deterministic reveal/redact/omit decisions; explicit correlation-token
  safety; privilege downgrade; and replay-completeness effects.

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
- Gated sampling is zero-based and keyed by a caller-supplied logical
  opportunity ordinal: `every(N)` admits `0, N, 2N, ...`. Trace and info use
  independent cadences; warn and error cannot be suppressed. A rejected
  opportunity invokes neither the payload/clock builder nor `Emitter::emit`,
  so only returned `Some(Event)` values advance the ordinary emitted-event
  sequence.
- Process capture never drops a critical frame. Queue/inline pressure returns
  the untouched frame for drain/spill/retry; cancellation or sink failure
  returns incomplete evidence, unchecked integrity, and demoted promotion.
  Diagnostic/telemetry loss is consumed only after a quantified gap is
  retained; a full gap ledger applies backpressure.
- `ProcessFrame::into_event` preserves the frame's severity in the canonical
  event envelope and hex-encodes opaque inline bytes only for display; typed
  event identity binds the original bytes. `ProcessGap::into_event` preserves
  the full ordinal range, u128 counts, reason, policy version, and optional
  committed artifact pointer. Malformed ranges/counts fail the event lint.
- A `DurableArtifactPointer` can be constructed only through the committed
  pointer constructor. Inline omission without that token remains explicit
  loss; process exit (including code zero) without a final typed receipt
  remains an observation gap.
- `xtask check-casual-print` scans every core-library Rust source and rejects
  `print!`, `println!`, `eprint!`, `eprintln!`, and `dbg!` outside test and CLI
  boundaries. Its three pre-policy structured emitters are ratcheted by exact
  path plus function owner; listing a file never exempts unrelated functions.
- Share policy never infers sensitivity from payload text. Opaque bytes are
  revealed only when sensitivity, audience, privilege, license, export,
  retention, and request bounds all admit them. Credentials are never
  revealed, including locally; unsafe unsalted/salted correlation of PII or
  secrets refuses the whole manifest before disclosure.
- Share entries are sorted by path and duplicate paths refuse. Intentional
  redaction preserves integrity but makes required replay evidence incomplete
  and promotion demoted; missing license/export authority is unsupported.
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
event-validation boundaries; they reject rather than repair. Process-capture
constructors reject zero bounds, malformed identities/hashes, and zero
ordinals; non-monotone frames are returned untouched. Privacy constructors
reject malformed paths/realms/tokens, zero bounds, duplicate paths, byte/count
budget overflow, and dictionary-testable correlation for protected labels.
`SamplingCadence::every` accepts only `NonZeroU64`; a zero cadence is therefore
not representable. Admission uses divisibility without incrementing the
logical ordinal, so it is defined at `u64::MAX` and introduces no hidden
opportunity-counter overflow. Callers own ordinal epochs and any rollover
policy.

## Determinism class
Deterministic: pure functions; no clocks (callers supply `wall_ns`), no I/O,
no RNG. Gated decisions depend only on severity, cadence, and the explicit
logical ordinal. Producers must derive that ordinal from stable logical work
identity rather than worker, thread, or arrival order.

## Cancellation behavior
All operations O(event size). No Cx required. Process runners project external
cancellation into `CaptureCancellation`; the pure capture state holds no
lifecycle/control lock while applying backpressure.

## Unsafe boundary
None.

## Feature flags
None.

## Conformance tests
The unit and integration batteries cover all 35 kinds' serialize+validate
round-trip; an ordered count-locked mutation of all 154 payload fields; golden
line; envelope/content hash split; every top-level semantic event field;
independent event-identity-version and wire-schema bytes/roots;
same-display/different-bit NaNs; exact opaque Custom bytes; retained receipt
admission and each refusal class; stale identity versions; monotone sequences;
zero-based and maximum-ordinal sampling; independent trace/info cadences;
unconditional warn/error admission; lazy rejected builders; emitted-only
sequence advancement; corruption rejection; failure-record lint;
hostile-string escaping; non-finite tagging; FNV known answers; and exhaustive
replay-frame mutations, including direct domain-prefix binding. The public
conformance case also proves gated events retain the canonical wire and exact
identity contracts while direct `emit` continues the same sequence stream.

The process-capture battery covers critical queue pressure and cancellation,
sink failure, deterministic telemetry sampling, diagnostic drop coalescing,
bounded gap-ledger pressure, durable oversized spill, lossy truncation,
non-UTF-8 payload preservation, per-stream ordering, and final-receipt closure.

The privacy battery covers public and privileged-local projection, credentials
in every privilege state, PII/secret dictionary attacks, keyed opaque
correlation, license/export refusal, expiry/legal hold, deterministic field
ordering, hostile binary values, duplicate paths, budget exhaustion, and
privilege downgrade.

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
  `fs_obs::validate_line` (never a second dialect).
- `emit_gated` guarantees that rejection does not invoke its payload/clock
  builder, clone the emitter envelope, materialize an `Event`, or serialize it.
  It does not assign logical ordinals, run a sink, or promise that unrelated
  caller work is optimized away. No ns/op or end-to-end percentage overhead
  is claimed until an ignored release roofline harness records an explicit
  machine fingerprint and acceptance band.
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
- `process` is the deterministic policy/state layer, not a pipe reader,
  artifact store, async executor, CLI renderer, or source print-macro scanner.
  The repository scanner is `xtask check-casual-print`; runtime owners must
  drain and persist the canonical events produced by the projection APIs and
  independently authenticate committed artifacts. `DurableArtifactPointer`
  records constructor-level admission; it cannot prove an external store
  truthful by itself.
- `check-casual-print` lexes workspace Rust source fail-closed, treats comments
  as trivia, recognizes structurally attached exact `cfg(test)` item scopes,
  rejects protected-macro renames, and ratchets each pre-policy emitter by one
  unique function body plus exact normalized invocation tokens. It makes no
  semantic claim about output hidden inside an external or generated macro
  whose protected spelling is absent from workspace source.
- `privacy` is a labeled-data policy core, not a secret detector, cryptographic
  implementation, access-control service, jurisdiction engine, or artifact
  scanner. External salted/keyed tokens remain caller assertions until an
  admitted cryptographic owner verifies them. Legal hold controls retention;
  it never grants disclosure.
