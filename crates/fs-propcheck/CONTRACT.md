# CONTRACT: fs-propcheck

> The Gauntlet G0 property-testing engine (bead frankensim-4nh8): seeded
> generators + integrated shrinking, adopted by law suites across layers.

## Purpose and layer
In-house property-based testing with proptest-class shrinking. Layer: UTIL
(zero runtime dependencies; dev-dependency consumers only).

## Public types and semantics
- `Stream` — deterministic splitmix64 test-input source keyed by
  `(suite_seed, case_index)`. NOT a simulation RNG: kernels keep fs-rand's
  Philox discipline; this stream never touches ledgers or physics.
- `Shrink` — the fs-bisect::compound convention verbatim: candidates
  strictly smaller, most aggressive first, fixed order, empty = fixpoint.
  Provided for i64/u64/f64/Vec/T-pairs; law suites implement it for their
  own case types.
- `minimize(seed_failure, property, max_steps) -> ShrinkReport` — greedy
  first-failing-candidate descent; identical trajectory for identical
  inputs (deterministic debugging).
- `check(name, suite_seed, cases, generate, property)` — the runner: on
  failure, shrink, emit one JSONL row, panic with the replay seed.
  `FSIM_PROPCHECK_REPLAY=<case_seed>` reruns exactly the failing case.

## Invariants
1. Same `(suite_seed, case_index)` reproduces the identical generated
   input on every platform (splitmix64 is integer-exact).
2. Every provided `Shrink` impl strictly decreases its measure (|x| for
   scalars; (len, elementwise |x|) lexicographic for vectors), so greedy
   descent terminates without a step budget in the non-adversarial case;
   the budget (10k) is a backstop, and `converged: false` is reported
   honestly when hit.
3. A failing `check` never exits silently: JSONL row + panic carrying the
   minimized counterexample and the replay seed.

## Error model
Test-harness semantics: inverted generator bounds are caller bugs and
panic with the offending bounds. The runner's failure panic IS the API.

## Determinism class
Deterministic: bit-stable case streams and shrink trajectories across
runs, thread counts, and ISAs (integer arithmetic only in the engine).

## Cancellation behavior
None: test-time only, no long-running compute paths (a single check's
budget is bounded by cases x property cost, owned by the caller).

## Unsafe boundary
None. `unsafe_code` denied.

## Feature flags
None.

## Conformance tests
tests/propcheck_battery.rs: stream determinism + decorrelation, generator
bounds (10k draws incl. special values), shrink-lattice strict decrease,
the seeded-violation self-test (a planted (a>=100, b>=7) law break must
minimize to exactly (100, 7)), clean-property run (500 cases), and the
full failing-path drill (panic carries kernel + replay seed).

## No-claim boundaries
No statistical coverage claims (generation is uniform-ish, not
distribution-calibrated); no concurrency; no simulation-RNG duties; the
f64 generator draws special values by design and makes no density claim.

## Golden couplings (bead y4pt)
No goldens pinned. Any future pinned shrink-trajectory golden must have a
row in golden-couplings.json per docs/GOLDEN_POLICY.md.
