# CONTRACT: fs-casebook

The shared conformance-case harness (bead huq.5, plan §13.3): named cases
with structured JSON-lines records — the executable half of the CONTRACT.md
discipline. A crate's conformance suite registers cases; running the suite
emits one record per case (case id, inputs digest, tolerance, verdict,
evidence pointers). Replayable cases additionally emit a companion containing
the exact rerun command and complete canonical input frame; disagreements emit
a first-mismatch record suitable for an upstream bug report. A reimplementation
can therefore be held to its predecessor's suite without weakening the stable
v1 case-record schema.

## Purpose and layer

Layer UTIL (test support, like fs-propcheck). Zero runtime dependencies:
policy tooling and any crate's dev-dependencies can use it without pulling
library builds.

## Public types and semantics

- `fnv1a64(bytes) -> u64` — the canonical inputs-digest helper.
- `ToleranceSpec` — `Exact | Ulps(n) | RelativeLe(b) | AbsoluteLe(b) |
  Structural`; rendered stably into records (`exact`, `ulps<=n`,
  `rel<=b`, `abs<=b`, `structural`).
- `CaseOutcome { pass, details, evidence, disagreements }` with
  `pass/fail/with_evidence` constructors — what one executed case reports back;
  `with_disagreement` attaches a typed disagreement and forces the outcome red.
- `Suite::new(name).case(id, inputs_digest, tolerance, run).run()` —
  registration-ordered deterministic execution; returns `SuiteReport`.
- `ReplaySpec::new(command, canonical_inputs)` / `from_hex` — a stable rerun
  selector plus the complete canonical input frame. `case_replayable` derives
  the case digest from these bytes and emits a `ReplayRecord`; callers cannot
  independently supply a drifting digest.
- `ReplayRecord::verify_and_decode` — schema/required-field admission followed
  by canonical-lowercase hex reconstruction and declared-length/FNV
  verification. `verify_and_decode_for` additionally checks the independently
  selected suite and case identity. `ReplayError` localizes version, identity,
  required-field, hex, length, and digest failures.
- `DisagreementRecord::first` — exact-frame comparison returning `None` for
  equality or a record with complete lengths/digests and the first differing
  byte offset. A common-prefix length mismatch records `None` for the side that
  ended first. Its comparison evidence is immutable after construction and is
  exposed through read-only accessors; `json_line()` renders deterministic
  bug-report-ready JSONL.
- `CaseRecord` — the typed per-case record; `json_line()` renders the
  one-line JSON form with deterministic field order and full escaping.
- `SuiteReport { records, replay_records, disagreements }` — `all_passed()` (an
  empty suite is NOT green), `failures()`, unambiguous replay/disagreement
  lookup, and `assert_green()` (panics carrying each failing registration plus
  that exact registration's replay and disagreement rows).

## Invariants

- Execution order is registration order; records preserve it.
- Every record carries the suite name, stable case id, 16-hex-char inputs
  digest, rendered tolerance, verdict, details, and evidence pointers.
- `Suite::case` and `CaseRecord::json_line()` retain their exact v1 behavior and
  bytes. Replay metadata is additive companion JSONL, never an in-place schema
  mutation of legacy case rows.
- Every replayable-case digest is derived from its retained canonical bytes.
  Decoding a replay record must admit its schema and required fields and
  reproduce those bytes, length, and digest or fail closed with a typed
  `ReplayError`; matching it to a selected case also requires exact suite/case
  identity.
- Disagreement localization reports the first unequal byte, or the first
  absent byte when equal prefixes have different lengths. It binds both full
  frames by length and FNV even though it retains only the localized bytes.
- Attaching any disagreement forces its owning `CaseOutcome` red. Private
  registration/discovery ordinals bind every companion to the exact execution
  that produced it without changing stable JSON bytes; the report admits only
  ordered, in-bounds, identity-consistent companions and refuses a green record
  carrying disagreement evidence.
- Duplicate case ids and empty case ids are recorded as structural
  FAILURES at run time (fail closed), never silently accepted. String lookup
  for a duplicate id is intentionally ambiguous and returns no companion;
  registration-indexed failure rendering still retains each duplicate's own
  evidence.
- An empty suite is not green: running nothing proves nothing.
- `json_line()` output is exactly one line; all string content is escaped.
- The case format is data-first: outcomes and reports are values, printing
  is a separable layer, so an IR-speaking front end (post fs-ir-core) wraps
  additively without rewriting suites.

## Error model

Legacy registration remains infallible: static ids and closures; defects
discovered at run time (duplicate/empty ids) become failing records rather
than panics. Replay admission returns typed `ReplayError` values for unsupported
schemas, empty required fields, selected-identity mismatch, or frame-integrity
failure.
`assert_green` is the deliberate merge-gate panic point, used by test mains to
fail the process with all relevant structured rows. `Suite::new` panics on an
empty suite name and `ReplaySpec::new` on an empty command (programmer errors).

## Determinism class

Deterministic: execution order, record field order, digest arithmetic
(FNV-1a 64), canonical lowercase hex, first-byte localization, and rendering
are pure and platform-independent. Float formatting in caller-supplied
`details` strings is the caller's claim, not this crate's.

## Cancellation behavior

None: case execution is synchronous and bounded by the caller's own cases.
Long-running suites belong under the caller's Cx lanes; the harness adds no
blocking, no I/O beyond stdout emission of case/replay/disagreement rows in
`run()`.

## Unsafe boundary

No `unsafe` (workspace forbids it; nothing here needs it).

## Feature flags

None.

## Conformance tests

`tests/casebook.rs` — the demo suite (the copyable example of adding cases:
exact-roundtrip, numeric-tolerance, structural-refusal); the intentionally
failing self-test; fail-closed empty and duplicate-id suites with exact
per-registration companion ownership; exact legacy `CaseRecord` JSON bytes,
escaping, and FNV constants; stable replay companion rows; schema,
required-field, selected-identity, canonical-frame, length, and digest tamper
refusal; identical-frame no-disagreement; seeded first-byte corruption
localization; common-prefix length-boundary localization; immutable
disagreement evidence and owning-case identity binding; companion JSON
escaping; report-integrity red controls; and `assert_green` retention of the
bug-report-ready rows.

## No-claim boundaries

- The exact-frame comparator localizes bytes; it does not decide what a
  numerical tolerance means physically or compare records across hosts by
  itself. Cross-ISA evidence remains a Gauntlet/fs-detaudit lane.
- A replay record is a complete deterministic recipe artifact, not a command
  executor or sandbox. The harness does not launch subprocesses.
- Disagreement JSON is bug-report-ready evidence; this crate does not choose
  which upstream is wrong and does not file issues or mutate external systems.
- No certification tiers: converter certification by sheaf axioms is
  fs-conform's scope.
- IR admission/ledger wrapping is fs-ir's scope; this dependency-free harness
  remains usable by non-IR test suites.
- Records are stdout JSON lines, not authenticated ledger rows; ledger
  binding is fs-obs/fs-ledger scope.
