# CONTRACT: fs-roofline

> Status: ACTIVE (harness v0). Owns the roofline measurement discipline of
> plan §14: measured axes, intensity-derived limits, dispersion-reported
> attainment, fingerprint-keyed ledger rows, staleness detection.

## Purpose and layer

Performance claims as falsifiable targets: benchmark every registered
kernel against its arithmetic-intensity-derived limit on the actual
machine, ledger the result under the machine fingerprint, and alert when
the fingerprint drifts. Layer: L6 (consumes fs-substrate, fs-simd,
fs-exec, fs-session, fs-ledger). Runtime dependencies remain workspace
crates plus `std`.

## Public types and semantics

- `MachineAxes` (`axes::probe()`) — measured axes only: STREAM-triad
  bandwidth (single/all-core, from fs-substrate) and FMA-chain peak FLOPs
  (single/all-core, in-house microbench), plus the fs-substrate topology
  fingerprint. Never spec-sheet numbers.
- `KernelSpec` — identity + intensity model (`bytes_per_elem`,
  `flops_per_elem`), threading axis, optional `target_fraction`.
- `RooflineKernel` — owns its buffers; `run_once` is the timed unit.
- `measure` / `run_registry` — warmup + repetitions →
  [`Attainment`]: median rate, achieved GB/s and GFLOP/s, binding
  `RoofSide`, attainment fraction, relative IQR `dispersion`, and a
  `Verdict` (`WithinBand`/`BelowBand`/`NoTarget`/`EnvironmentInvalid`).
  The invalid verdict carries a reason and is never a performance pass or
  failure.
- `run_admission_error` / `run_is_citable` — the publication boundary:
  require exact pre/post axis agreement, private timed provenance, positive
  work count and sample durations, raw-sample re-derivation, unmodified
  spec/derived fields, and unique kernel identities. Analytic helper rows are
  deliberately non-citable.
- `SECTION_14_1_TARGETS` — the plan's target table as data; `landed`
  flips only when the owning kernel registers here (no silent coverage
  gaps).
- `record_run` — one atomic ledger transaction. Admitted runs publish
  metrics, `benchmark_result` events, and `tune` rows keyed (kernel ×
  `roofline-v2:<kernel-version>` × fingerprint LE bytes). Rejected input
  publishes one rejection event and an Error op, never normal-looking metrics
  or tuning evidence; storage failures roll back the entire write set.
- `staleness` — `MatchingIdentityAgeUnknown` / `FingerprintDrift` /
  `NeverMeasured` per kernel version and fingerprint. The name deliberately
  does not claim freshness because the tune schema has no timestamp.
- `kernels::default_registry` — the stable test/meta registry: fs-simd
  axpy/dot/sum (report-only bands in v0). `SeededSlowKernel` is the separate
  meta-test kernel claiming a band it cannot meet.
- `kernels::production_registry` / `GemmKernel` — the shipped command's
  registry adds real f64 GEMM through `fs_session::gemm_f64_session`.
  One kernel instance owns one tuner and cancellation gate: the first
  warmup closes measure → cache → model → dispatch and later warmups/timed
  reps reuse its validated row. With `--ledger`, the kernel also owns a
  dedicated connection during measurement, adopts a valid row from a prior
  process, or persists the cold sweep before dispatch. Execution is
  sequential through exclusive `&mut RooflineKernel`; no wrapper lock hides
  tune state.
- `roofline` CLI bin — axes line, per-kernel JSONL, §14.1 coverage table,
  optional `--ledger` recording + staleness report.

- `regress` module (plan §14.4, bead fz2.4): the regression layer.
  `gate` — DISPERSION-AWARE bands (k·σ against the rolling baseline,
  never a naive threshold), and a red arrives WITH its diagnosis: the
  phase-share flame-graph diff ranked by growth. `Cusum` — the
  complementary slow-drift detector (slack k, threshold h) over
  expanding-baseline standardized scores. `slower_this_month` — the
  canonical dashboard question as ONE call: (kernel, drop %, guilty
  phase). Calibration is meta-tested: zero false alarms across 20
  kernels × 60 stable nights at the default settings; a 0.3σ/night
  drift invisible to the single-night band trips the CUSUM mid-month.

## Invariants

1. Axes are measured on the machine that runs the kernels, in the same
   process; the compute axis is compiler-achievable FMA throughput
   (conservative where autovectorization misses — the honest direction for
   a denominator). Probe calibration (bead xdgf): timed samples span
   ≥ 5 ms (single microsecond-scale passes sat inside the frequency-
   ramp/scheduler noise floor and wandered tens of percent), and the
   accumulator lane count is REGISTER-FILE-sized per arch (48 on
   aarch64, 64 on x86) — the former 64-lane constant spilled on NEON
   and read the axis ~25% low, which inflated attainments past 1.0.
2. `attainment = measured_rate / min(bandwidth_limit, compute_limit)` with
   limits derived from the spec's intensity model (meta-tested against
   hand calculations).
3. Receipt schema v1 carries bit-exact pre-run axes, intensity spec, target,
   every raw timed sample, median/p25/p75/dispersion, and exact derived-result
   bits. Rounded decimal fields are display-only. A standalone reader can
   rederive the reported rate, roof, and variance bar.
4. Ledger rows are keyed by kernel version and fingerprint. A drifted
   fingerprint or version refuses reuse; matching identity has unknown age.
5. Axes must be finite and positive, have a nonzero logical-CPU count, meet
   the 5 GB/s and 5 GFLOP/s single-thread reference-family floors, and have
   aggregate axes at least half their single-thread counterparts. These
   absolute guards catch the extreme bead-1n61 collapse. A second axis probe
   after the registry must also agree within 25% on every axis; changing
   contention poisons the run.
6. Specs, rates, targets, and dispersion are screened before verdict
   arithmetic. Any non-finite/negative input or attainment above 1.5 makes
   the run invalid. One invalid row poisons every verdict in that registry
   run because the shared axes can no longer certify any sibling result.
7. Citable rows are created only by `measure`: `elements > 0`, every elapsed
   sample finite and positive, stored sample count equal to `reps`, and all
   public outputs bit-exactly rederive from the private receipt and exact axes.

## Error model

Measurement APIs are infallible (they report what they saw, including
zero rates, with invalid evidence normalized to finite JSON plus an
explicit reason). Ledger interaction returns `fs_ledger::LedgerError`
(structured, machine-actionable). The CLI refuses malformed arguments with
a structured JSON error on stderr and a nonzero exit.

## Determinism class

Not deterministic: wall-clock measurement of a shared machine. The
REPORTING is deterministic given the same measured times (order statistics
with deterministic tie-breaking). Seeds are not applicable; repetition
counts and dispersion make the noise visible instead of hidden.

## Cancellation behavior

The session-backed GEMM carries an `fs_exec::CancelGate` into fs-la's
bounded-poll, request → drain → finalize path. The current synchronous
roofline CLI owns but does not externally expose that gate; cancellation
therefore cannot be requested through the command yet. Other registry
kernels remain bounded by `reps × run_once` and have no tile cancellation
surface.

## Unsafe boundary

None. Safe Rust only; SIMD reached through fs-simd's safe façades.

## Feature flags

None. All v0 behavior is `[S]` default-path.

## Conformance tests

`tests/conformance.rs`: registry run + reporting shape (rf-001);
seeded-slow kernel demonstrably below band on real axes (rf-002);
ledgered run with versioned fingerprint-keyed tune rows, lint-clean (rf-003);
identity-match-age-unknown/drift/never-measured plus rejection-without-publication
(rf-004/004b); re-run reproducibility
within stated dispersion allowance (rf-005); CLI smoke incl. §14.1
coverage table and structured refusals (rf-006). Unit tests cover
attainment hand-calculations, order statistics, and axes sanity. The GEMM
registry regression executes the production session call twice and proves
exactly one cold sweep plus two recorded dispatch decisions (warm-row reuse,
not a test-only wrapper).

## No-claim boundaries

- No machine performance claims live in this crate: numbers become
  citable only as ledgered rows with fingerprints (plan §14.1 discipline).
- The compute axis is compiler-achievable FMA throughput, not theoretical
  ISA peak; modest attainment above 1.0 is reported as-is. Attainment above
  1.5 means the shared axis is not credible for gating, whether because it
  was crushed/stale or because a specialized kernel outran the probe by too
  much; the run is retained as invalid evidence and must be re-probed.
- §14.1 family targets other than GEMM remain `landed: false` until their
  kernels register. The GEMM production kernel is registered; the target
  table's `landed` bit remains false until the owning concurrent edit is
  reconciled with its focused execution proof.
- Per-CCD bandwidth axes, P/E-core-class split, frequency-state capture,
  and thermal controls are future scope (v0 measures whole-machine axes).
- Static floors plus pre/post agreement cannot detect a host that is already
  degraded before the first probe and remains equally degraded through the
  second. Fingerprint-specific trusted historical baselines are still required
  before promoting this harness from reporting to a binding nightly gate.
- `RooflineKernel::elements()` and intensity are asserted by the registered
  implementation. Receipt v1 proves what was timed and how the arithmetic was
  derived; it does not prove a custom trait implementation performed the work
  it claimed. Default-registry review is the v1 trust root; implementation
  hashes remain follow-up scope.
- Verdict gating in CI is deliberately absent on shared runners; bands
  bind only on ledgered reference machines (nightly lane, later).

## Fail-closed evidence screening (bead fz2.4.1)

Every public regress entry point screens its floating inputs before
any verdict arithmetic: `gate` returns `GateVerdict::Invalid { reason }`
— never Green — for non-finite or negative attainment or phase
durations anywhere in the history and for unusable specs (non-finite
or non-positive k_sigma, min_baseline < 2); `Cusum::first_alarm`
alarms AT the first non-finite residual (NaN previously reset the
shortfall via `max`, silently suppressing detection) and an invalid
detector spec cannot certify quiet; `standardize` maps history to −∞
from the first non-finite entry so poison never enters the expanding
baseline; `slower_this_month` reports poisoned kernels FIRST with an
infinite drop and the flaw as the "why", never skipping them.
Metamorphic property (tested): rescaling phase durations by a constant
(time-unit change) preserves verdicts and attribution ranking.

## No-claim boundaries (regress)

- This module is the STATISTICS + attribution + gate arithmetic; the
  nightly both-machine CI wiring rides the ci-gauntlet pipeline
  (huq.4's lane), memory-regression tracking rides fs-alloc's
  allocation-site diffs, and FrankenPandas trend dashboards ride
  fs-report — the named seams, each consuming these verdicts.
- Suspect-commit bisection hooks are diff-vs-last-green consumers of
  the same attribution output, not re-implemented here.
