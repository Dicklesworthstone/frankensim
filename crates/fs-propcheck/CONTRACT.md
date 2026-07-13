# CONTRACT: fs-propcheck

> The Gauntlet G0 property-testing engine (bead frankensim-4nh8): seeded
> generators + integrated shrinking, adopted by law suites across layers.

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
tests/propcheck_battery.rs plus in-crate parser/escaping tests: stream
determinism; unbiased integer rejection and extreme/adjacent finite f64 bounds;
explicit caller-boundary refusals; shrink-lattice strict decrease over
`i64::MIN`, singleton vectors, infinities, and coordinated tuples; passing-seed
refusal, exact-budget convergence, and evaluation/candidate work ceilings; the
planted `(a>=100,b>=7)` minimum;
clean-property coverage; returned-false and caught-panic failure drills with
replay-file diagnostics; deterministic artifact path/JSONL writing; malformed
replay parsing and JSON escaping.

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

## Golden couplings (bead y4pt)
No goldens pinned. Any future pinned shrink-trajectory golden must have a
row in golden-couplings.json per docs/GOLDEN_POLICY.md.
