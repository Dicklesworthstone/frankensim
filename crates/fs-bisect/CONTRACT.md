# CONTRACT: fs-bisect

Physics-VCS bisect (plan addendum, Proposal 10): git-bisect for a wrong number.
Plus the Gauntlet failure-compounding workflow (`compound`, bead 6nb.9):
minimize every counterexample into a permanent regression family.

## Purpose and layer

Layer L6 (version control / orchestration). Runtime dependencies are the
in-tree `fs-blake3` content hash and dependency-free UTIL `fs-propcheck` shrink
trait — no numerical stack. The implementation remains pure control flow over
a caller-supplied `CommitOracle` and data plumbing for `compound`'s capture →
minimize → probe → family → replay pipeline.

## Public types and semantics

- `Verdict { Good, Bad }`; `CommitOracle::evaluate(commit) -> Verdict` —
  commit 0 is the oldest (assumed-good baseline), `len-1` the newest. The
  oracle IS "replay commit k, then evaluate the predicate".
- `bisect(len, &oracle) -> BisectRun` — `O(log n)` binary search assuming a
  monotone predicate; returns `Culprit { index, confirmed: false }`,
  `AllGood`, `AllBad`, or `Empty`. `BisectRun { result, probes }` logs the
  search path.
- `verify_monotone(len, &oracle) -> Option<(usize, usize)>` — `O(n)` scan for a
  non-monotonicity witness (a Bad followed by a later Good).
- `bisect_checked(len, &oracle) -> BisectRun` — verifies monotonicity first;
  a non-monotone predicate yields `NonMonotone { bad, later_good }` instead of
  a mis-localization.
- `bisect_two_tier(len, &low, &full) -> BisectRun` — narrows with a cheap
  `low`-fidelity oracle, then CONFIRMS the culprit at `full` fidelity; if full
  rejects the low candidate it re-searches entirely at full fidelity. The
  culprit is `confirmed = true` (a *verified* localization vs the *estimated*
  single-fidelity one).
- `compound` module (bead 6nb.9): `FailureCase<I>` (id, seed, typed input,
  `InvariantClass`, contract surface, detail); public re-export of
  `fs_propcheck::Shrink` (one deterministic candidate contract across G0 and
  permanent failure families) + `minimize` (`fails == true` means FAILURE;
  greedy first-failing descent to a fixpoint or
  step budget, `converged` flag, typed `NotFailing` refusal on a passing input,
  and hard per-step/aggregate evaluation ceilings including the seed call);
  `probe_neighborhood` over a count- and aggregate-byte-bounded, uniquely
  labeled caller set; internally
  constructed `RegressionFamily<I>` (minimum + failing neighbors, nonempty
  unique tracking refs, `recommended_admission`) only after minimization reaches
  a fixpoint. The accepted-step budget still checks whether the value reached
  at the exact cap is a fixpoint. `Canon` requires a stable globally unique
  codec id/schema, a pure encoding of every predicate-relevant field, writes
  through a bounded `CanonWriter`, and supplies the payload snapshots that
  `content_hash` binds with the complete codec schema using
  domain-separated in-tree BLAKE3. The escaped JSON-lines `manifest` exposes
  codec/schema fields and a hash trailer. `replay` re-canonicalizes every live
  member globally before work, immediately before and after its predicate,
  and again after all callbacks; candidate generation, retained minimization
  witnesses, and cross-neighbor shared state receive the same checks, so
  persistent identity drift is refused. Members
  that authentically stop failing are REPORTED as `now_passing`, never silently
  dropped. `COMPOUND_CANON_VERSION` is the golden-couplings surface const: any
  codec/tag/hash/manifest assembly change bumps it and deliberately re-freezes
  dependents.

## Invariants

- On a monotone sequence with a Good prefix and a Bad suffix, `bisect` returns
  the first Bad index.
- `bisect_two_tier` never returns a full-fidelity-rejected candidate: it
  re-searches. Its culprit is always `confirmed`.
- All functions are pure and deterministic; the probe log records every
  evaluation in order.
- Failure compounding validates identifiers and descriptions before descent,
  rejects empty tracking, duplicate labels/references, reserved custom
  invariant names, the reserved `minimized` neighbor label, incomplete
  minimization, and limit+1 work. Content identity uses the type/schema domain
  and canonical bytes sealed at family construction. Callback bracketing and a
  final replay pass refuse persistent canonical mutation; authoritative replay
  additionally requires the documented pure-and-complete `Canon` contract.

## Error model

No errors/panics on valid indices; degenerate inputs map to explicit result
variants (`Empty`, `AllGood`, `AllBad`, `NonMonotone`). Failure compounding
returns `CompoundError::{NotFailing, InvalidField, LimitExceeded,
DuplicateIdentity, MinimizationIncomplete, ReplayIdentityDrift,
CallbackIdentityDrift}`; it never
silently truncates, repairs, or replays an unauthenticated family.

## Determinism class

Fully deterministic: a bisect is a pure function of `(len, oracle)`; the same
oracle reproduces the same culprit + probe path (sound only if the oracle's own
commit replay is deterministic — the ledger's `at(t)`/ExecMode contract).
Compounding identity is BLAKE3 over the versioned codec/schema domain and
construction-time canonical snapshot. Persistent interior mutation leaves the
identity stable but makes replay fail closed at a pre/post/final check.

## Cancellation behavior

None here; the oracle's own (possibly expensive, low- or full-fidelity)
evaluation runs under the caller's cancellation scope.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/bisect.rs` (Proposal 10, 10 cases): first-bad localization (short +
long, `O(log n)` probe count); empty / all-good / all-bad / singleton
boundaries; `verify_monotone` witness; `bisect_checked` flags non-monotone;
two-tier agreement, re-search on full-fidelity rejection, and endpoint
confirmation; confirmed-flag semantics; determinism.

`tests/compound_battery.rs` (6nb.9 + j3q2 acceptance, 22 cases): a deliberately
broken cross-crate golden modeled on the real powi incident (sequential vs
square-multiply chains; minimizes to the exact k=7 divergence boundary for
base 0.7, neighborhood shows the sharp region edge, family replays live and
goes stale under the fixed implementation); a falsifier hit on a wrong tail
constant (systematic error minimizes to n=1, whole neighborhood fails);
frozen content identity
`6aed2aab4250ca30b657e6e016f628eb5ba33e09c3e49b732156a5058eb9141f`
and complete JSON-lines manifest hash
`43bd1ebddb606eb5a156ca06c642cf03ddd06770223584c5ea45ae031d0cf6b2`
(registered in `golden-couplings.json` against
`fs-bisect:compound-canon=3`); bitwise
minimization determinism + `NotFailing` refusal; canon concatenation-collision
resistance; stable codec-domain separation for identical payload bytes;
per-field content-hash sensitivity; exact work boundaries and limit+1
refusals, including the aggregate seed-inclusive evaluation ceiling; nonempty
tracking/reserved-invariant/canonical-byte gates; incomplete-minimization
refusal; parsed unique-key JSON-lines and escaping; exact-cap fixpoint
recognition; fail-closed predicate/neighbor/replay mutation; and an executable
no-claim fixture for intentionally incomplete codecs.

The battery's `serde`/`serde_json` dev dependencies are an isolated independent
parser oracle for the emitted JSON-lines manifest and duplicate-key rejection;
they are not linked into `fs-bisect` production code.

## No-claim boundaries

- `bisect` ASSUMES monotonicity (documented); use `bisect_checked` when the
  predicate may be non-monotone. Detection is `O(n)`; plain `bisect` stays
  `O(log n)`.
- The colors (estimated for the low-fidelity search, verified for the
  full-fidelity confirmation) are represented here by the `confirmed` flag; the
  caller attaches the `fs-evidence` `Color` when it records the result.
- Commit replay determinism is the ledger's contract; this crate assumes the
  oracle is a faithful replay-plus-predicate.
- `compound` does not write to the ledger and emits no fs-obs events —
  recorded follow-up once the huq.16 observability schema lands; manifests are
  returned to the caller, not persisted here.
- `compound` does not enact admission rules: `recommended_admission` is
  carried prose (as check-powi was born from the powi incident); enacting it
  is the responding agent's task.
- `minimize` finds A minimal failing input under the caller's shared `Shrink`
  order,
  not THE global minimum (greedy, not exhaustive); `converged = false` marks a
  budget-limited descent honestly. It deliberately retains a distinct,
  stricter operational contract than `fs_propcheck::minimize`: the predicate
  uses the opposite polarity (`true` means failure), and fs-bisect additionally
  enforces canonical callback bracketing, identifier validation, and hard
  candidate/evaluation limits.
- Caller callbacks (`fails`, `shrink_candidates`, `neighbors_of`, and `Canon`
  implementations) execute outside FrankenSim's control and must enforce their
  own cancellation and internal allocation budgets. `CanonWriter` bounds bytes
  retained by this crate before allocation, but cannot prevent a codec from
  doing excessive private work. Typed members themselves remain caller-owned
  allocations. Persistent callback mutation is detected, but a callback that
  transiently mutates and restores hidden state remains outside the claim.
- `Canon` completeness and purity are implementor obligations: Rust cannot
  prove that a caller omitted no predicate-relevant field. An incomplete codec
  may replay changed hidden semantics under unchanged bytes and is explicitly
  unauthenticated; authoritative families should use derived/field-by-field
  codecs over immutable replay values.
- Stable codec ids and versions are an implementor contract; this crate binds
  and validates their syntax but cannot prove that two independent owners did
  not choose the same id. A semantic codec change must change its schema
  version.
- The manifest carries codec/schema-qualified canonical bytes but this generic
  crate does not supply a decoder for arbitrary caller-defined `I`; replay uses
  the sealed in-memory typed family. Durable type-specific decoding belongs
  with each regression-family owner.
- The v3 manifest/content fixtures are reproduced in debug and release on the
  aarch64 Apple M4 Pro. Post-v3 x86-64 reproduction is pending and cross-ISA
  equality is not yet claimed.
