# CONTRACT: fs-casebook

The shared conformance-case harness (bead huq.5, plan §13.3): named cases
with structured JSON-lines records — the executable half of the CONTRACT.md
discipline. A crate's conformance suite registers cases; running the suite
emits one record per case (case id, inputs digest, tolerance, verdict,
evidence pointers) so failures reproduce from their log line alone and a
reimplementation can be held to its predecessor's suite.

## Purpose and layer

Layer UTIL (test support, like fs-propcheck). Zero runtime dependencies:
policy tooling and any crate's dev-dependencies can use it without pulling
library builds.

## Public types and semantics

- `fnv1a64(bytes) -> u64` — the canonical inputs-digest helper.
- `ToleranceSpec` — `Exact | Ulps(n) | RelativeLe(b) | AbsoluteLe(b) |
  Structural`; rendered stably into records (`exact`, `ulps<=n`,
  `rel<=b`, `abs<=b`, `structural`).
- `CaseOutcome { pass, details, evidence }` with `pass/fail/with_evidence`
  constructors — what one executed case reports back.
- `Suite::new(name).case(id, inputs_digest, tolerance, run).run()` —
  registration-ordered deterministic execution; returns `SuiteReport`.
- `CaseRecord` — the typed per-case record; `json_line()` renders the
  one-line JSON form with deterministic field order and full escaping.
- `SuiteReport { records }` — `all_passed()` (an empty suite is NOT
  green), `failures()`, `assert_green()` (panics carrying every failing
  record's JSON line).

## Invariants

- Execution order is registration order; records preserve it.
- Every record carries the suite name, stable case id, 16-hex-char inputs
  digest, rendered tolerance, verdict, details, and evidence pointers.
- Duplicate case ids and empty case ids are recorded as structural
  FAILURES at run time (fail closed), never silently accepted.
- An empty suite is not green: running nothing proves nothing.
- `json_line()` output is exactly one line; all string content is escaped.
- The case format is data-first: outcomes and reports are values, printing
  is a separable layer, so an IR-speaking front end (post fs-ir-core) wraps
  additively without rewriting suites.

## Error model

Infallible by construction at the API surface: registration takes static
ids and closures; defects discovered at run time (duplicate/empty ids)
become failing records rather than panics. `assert_green` is the one
deliberate panic point, used by test mains to fail the process with the
structured failure records in the message. `Suite::new` panics only on an
empty suite name (programmer error).

## Determinism class

Deterministic: execution order, record field order, digest arithmetic
(FNV-1a 64), and rendering are pure and platform-independent. Float
formatting in caller-supplied `details` strings is the caller's claim, not
this crate's.

## Cancellation behavior

None: case execution is synchronous and bounded by the caller's own cases.
Long-running suites belong under the caller's Cx lanes; the harness adds no
blocking, no I/O beyond stdout emission in `run()`.

## Unsafe boundary

No `unsafe` (workspace forbids it; nothing here needs it).

## Feature flags

None.

## Conformance tests

`tests/casebook.rs` — the demo suite (the copyable example of adding
cases: exact-roundtrip, numeric-tolerance, structural-refusal), the
intentionally-failing self-test (structured failure record + non-green
report + `assert_green` panic carrying the record), fail-closed empty and
duplicate-id suites, JSON escaping, and the digest helper's FNV-1a 64
constants.

## No-claim boundaries

- The harness runs and records; it does not decide what a tolerance means
  physically and does not compare records across runs or ISAs itself —
  cross-ISA evidence is a Gauntlet lane executing the same suite on both
  hosts and comparing emitted records.
- No certification tiers: converter certification by sheaf axioms is
  fs-conform's scope.
- No IR front end yet: the v0 format is designed for additive wrapping
  once fs-ir-core's integration language lands; no claim that wrapping
  exists today.
- Records are stdout JSON lines, not authenticated ledger rows; ledger
  binding is fs-obs/fs-ledger scope.
