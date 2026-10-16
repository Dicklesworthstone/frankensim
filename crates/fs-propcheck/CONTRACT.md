# CONTRACT: fs-propcheck

> The Gauntlet G0 property-testing engine (bead frankensim-4nh8) plus the G3
> metamorphic relation layer (bead frankensim-2uce): seeded generators,
> integrated shrinking, and structured relation evidence.

## Purpose and layer
In-house property-based testing with proptest-class shrinking. Layer: UTIL
(zero runtime dependencies; consumers normally use it as a dev-dependency,
while `fs-bisect` re-exports the same trait for production failure families).

## Public types and semantics
- `Stream` — deterministic splitmix64 test-input source keyed by
  `(suite_seed, case_index)`. NOT a simulation RNG: kernels keep fs-rand's
  Philox discipline; this stream never touches ledgers or physics. Inclusive
  integer ranges use rejection sampling; finite half-open f64 ranges avoid
  subtraction overflow and upper-endpoint rounding; case indices are streamed.
- `Shrink` — the shared convention publicly re-exported by
  `fs-bisect::compound`: candidates
  strictly smaller, most aggressive first, fixed order, empty = fixpoint.
  Provided for i64/u64/f64 (including ordered infinity descent), vectors, and
  coordinated pairs/triples; law suites implement it for their own case types.
- `minimize(seed_failure, property, max_steps) -> Result<ShrinkReport,
  MinimizeError>` — `property == true` means PASS; a passing seed is refused,
  and greedy first-failing-candidate descent checks for a fixpoint at the exact
  accepted-step budget. The wrapper also applies fixed 100k-evaluation and
  4,096-candidate-per-step ceilings.
- `minimize_with_budget(seed_failure, property, MinimizeBudget)` — the same
  descent with caller-selected accepted-step, evaluation, and per-step
  candidate ceilings. Hitting a ceiling retains the best known failing input
  and reports `converged: false`.
- `check(name, suite_seed, cases, generate, property)` — the runner: on
  returned-false or caught-panic failure, shrink, emit one dependency-free
  escaped JSONL row to stdout and a replay artifact, then panic with the
  shrunk input, artifact path, and replay seed. The default process-scoped
  artifact lives under `CARGO_TARGET_DIR` or the OS temp directory;
  `FSIM_PROPCHECK_REPLAY_FILE` selects an exact CI artifact path.
  `FSIM_PROPCHECK_REPLAY=<case_seed>` reruns exactly the failing case; malformed
  values fail closed instead of silently falling back to a full run.
- `metamorphic::RelationCase<Input, Transform>` — one generated base input and
  transform. Its `Shrink` implementation delegates to the coordinated tuple
  lattice, so a violation minimizes both components rather than freezing the
  original transform while shrinking only the input.
- `metamorphic::MetamorphicOperator` and `MetamorphicRelation` — dependency-free
  G3 harness traits. A relation owns a stable id, canonical family, transform,
  output expectation, and declared tolerance semantics. Closure-backed
  constructors cover rigid motion, coherent unit rescaling, refinement
  monotonicity, adjoint/finite-difference agreement, and conversion-path
  independence without making UTIL depend on domain crates.
- `metamorphic::check_relation(...)` — adapts a declared relation to the same
  case stream, replay selector, shrink budgets, and artifact writer as `check`.
  A returned violation's final shrunk failure row records relation id/family,
  transform, both outputs, tolerance, and finite signed admission margin in
  `failure_context`.
- `metamorphic::Tolerance` — exact finite-bit, absolute,
  absolute-or-relative, and non-increasing scalar semantics. Limits are
  nonnegative and finite; non-finite observations and overflowed comparisons
  fail closed. Nonnegative margins pass.

## Invariants
1. Same `(suite_seed, case_index)` reproduces the identical generated
   input on every platform (splitmix64 is integer-exact).
2. Every provided `Shrink` impl strictly decreases its measure (unsigned |x|
   for integers, including `i64::MIN`, and absolute magnitude for floats;
   (len, elementwise |x|) lexicographic for vectors), so greedy
   descent terminates without a step budget in the non-adversarial case;
   singleton vectors never yield themselves. The step, property-evaluation,
   and per-step candidate budgets are backstops; `converged: false` means a
   ceiling prevented the engine from proving a local fixpoint.
3. A failing or panicking `check` never exits silently: escaped JSONL stdout
   row + the same append-only replay row on disk + a panic carrying the
   best-known counterexample and replay seed. Artifact I/O failure is named in
   the final panic instead of being silently ignored.
4. `check_relation` never treats an inapplicable transform or operator panic as
   a skipped case. Consumer generators must produce admitted cases; panics are
   caught as failures, and a zero-case request is refused as vacuous. A
   returned relation violation is reevaluated after shrinking, and only that
   final input/transform/output/tolerance context is written to the receipt.
   Shrinking preserves the initial failure class, so a relation violation
   cannot minimize into an unrelated operator panic. Panic receipts retain the
   generic panic detail but cannot name outputs that were never produced.

## Error model
Test-harness semantics: inverted/non-finite f64 bounds and `vec_of(usize::MAX)`
are caller bugs and panic with the offending request. `MinimizeError::SeedPasses`
refuses fake failure seeds, and `EmptyEvaluationBudget` refuses an envelope
that cannot evaluate one. Malformed replay configuration panics before cases
run. The runner's final failure panic IS the API.

## Determinism class
Deterministic: splitmix/rejection case selection is integer-exact; candidate
order and shrink trajectories are fixed across runs and thread counts. Finite
f64 draws use an explicit IEEE-754 convex interpolation and endpoint clamp; the
engine makes no distribution-density or non-IEEE target claim.

## Cancellation behavior
None: test-time only. Case enumeration is streamed; each shrink descent has
explicit step, property-evaluation, and candidate-count ceilings. Property cost
and caller-defined `shrink_candidates` vector construction remain caller-owned.

## Unsafe boundary
None. `unsafe_code` denied.

## Feature flags
None.

## Conformance tests
`tests/propcheck_battery.rs`, `tests/metamorphic_battery.rs`, plus in-crate
parser/escaping tests: stream
determinism; unbiased integer rejection and extreme/adjacent finite f64 bounds;
explicit caller-boundary refusals; shrink-lattice strict decrease over
`i64::MIN`, singleton vectors, infinities, and coordinated tuples; passing-seed
refusal, exact-budget convergence, and evaluation/candidate work ceilings; the
planted `(a>=100,b>=7)` minimum;
clean-property coverage; returned-false and caught-panic failure drills with
replay-file diagnostics; deterministic artifact path/JSONL writing; malformed
replay parsing and JSON escaping. The G3 battery pins all five canonical
relation labels, fail-closed scalar tolerance boundaries, one declaration
reused by two operators, and a planted violation that jointly shrinks to
`input=100, transform=7` before checking its structured failure context.
Owning-crate G3 adopters exercise six production surfaces without removing
their fixed pins: `fs-la::gemm_f64`, `fs-fft::Fft::forward`,
`fs-geom::fixtures::SphereChart::convert`, both actions of
`fs-solver::CsrOp::general`, `fs-opt::descend_fn`, and
`fs-render::Lambertian::furnace_radiance`. Their relation ids, seeds, case
counts, tolerances, and no-claim boundaries live in the corresponding crate
contracts. Separately, `fs-conform` applies the same engine to its test-local
`Mtx` implementation of the public `Converter` SDK trait; that is SDK harness
coverage and is not counted as a production geometry conversion.

## No-claim boundaries
No statistical coverage claims for f64 generation (uniform-ish plus deliberate
specials, not density-calibrated); no concurrency; no simulation-RNG duties.
Shrinking is deterministic greedy local descent, not a global-minimum proof.
`check` can catch Rust unwind panics, not aborting processes. Replay selection
is process-global, so multi-property replay should target the named test.
Rust invokes the process-global panic hook before `catch_unwind`, so a
panic-based property may print intermediate shrink-probe panics even though the
runner catches and classifies them.
An explicitly shared `FSIM_PROPCHECK_REPLAY_FILE` is append-only but has no
cross-process record-locking claim; CI shards should use distinct paths. The
default process-scoped filename avoids that collision.
The engine count-checks a caller-defined shrink vector after construction; it
cannot preempt arbitrary work or allocation inside that trait method.
The G3 layer alone is not evidence that a production operator obeys a canonical
relation. Six owning-crate adapters now declare concrete production checks, but
their claims become executable evidence only when those crate suites pass; the
engine contract does not substitute for that batch proof. Broader domain
adapters for motors, runtime unit systems, refinement estimators, adjoint
directions, and certified conversion routes remain consumer-owned. A custom
relation comparator receives and receipts its declared tolerance but remains
responsible for applying those semantics honestly;
`Tolerance::evaluate_scalar` is the checked built-in path. The engine does not
infer transform applicability, physical units, certificate composition, or an
independent oracle from closure types.
Generic transform/operator/comparator panics are caught by the shared runner,
but this engine-only slice does not yet attach relation-stage context to a panic
that occurs before both outputs and a margin exist.

## Golden couplings (bead y4pt)
No goldens pinned. Any future pinned shrink-trajectory golden must have a
row in golden-couplings.json per docs/GOLDEN_POLICY.md.
