# CONTRACT: fs-roofline

> Status: ACTIVE (harness v0). Owns the roofline measurement discipline of
> plan §14: measured axes, intensity-derived limits, dispersion-reported
> attainment, fingerprint-keyed ledger rows, staleness detection.

## Purpose and layer

Performance claims as falsifiable targets: benchmark every registered
kernel against its arithmetic-intensity-derived limit on the actual
machine, ledger the result under the machine fingerprint, and alert when
the fingerprint drifts. Layer: L6 (consumes fs-substrate, fs-simd,
fs-ledger). Runtime deps: `std` + those three workspace crates.

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
  `Verdict` (`WithinBand`/`BelowBand`/`NoTarget`).
- `SECTION_14_1_TARGETS` — the plan's target table as data; `landed`
  flips only when the owning kernel registers here (no silent coverage
  gaps).
- `record_run` — one ledger op (frozen Five Explicits) + per-kernel
  metrics, `benchmark_result` events, and `tune` rows keyed
  (kernel × `roofline-v1` × fingerprint LE bytes).
- `staleness` — `Fresh` / `FingerprintDrift` / `NeverMeasured` per kernel
  against the current fingerprint.
- `kernels::default_registry` — fs-simd axpy/dot/sum (report-only bands in
  v0) and `SeededSlowKernel` (meta-test kernel claiming a band it cannot
  meet).
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
3. Every attainment row carries dispersion and repetition count; verdicts
   are reporting-only in v0 — no CI gate consumes them on shared runners.
4. Ledger rows are keyed by fingerprint; a drifted fingerprint makes every
   prior number stale, and `staleness` says so.

## Error model

Measurement APIs are infallible (they report what they saw, including
zero rates). Ledger interaction returns `fs_ledger::LedgerError`
(structured, machine-actionable). The CLI refuses malformed arguments with
a structured JSON error on stderr and a nonzero exit.

## Determinism class

Not deterministic: wall-clock measurement of a shared machine. The
REPORTING is deterministic given the same measured times (order statistics
with deterministic tie-breaking). Seeds are not applicable; repetition
counts and dispersion make the noise visible instead of hidden.

## Cancellation behavior

No long-running loops beyond `reps × run_once`; a run is bounded by its
arguments. Tile-level cancellation integration arrives when kernels run
under fs-exec scopes (deferred with fs-exec integration).

## Unsafe boundary

None. Safe Rust only; SIMD reached through fs-simd's safe façades.

## Feature flags

None. All v0 behavior is `[S]` default-path.

## Conformance tests

`tests/conformance.rs`: registry run + reporting shape (rf-001);
seeded-slow kernel demonstrably below band on real axes (rf-002);
ledgered run with fingerprint-keyed tune rows, lint-clean (rf-003);
staleness fresh/drift/never-measured (rf-004); re-run reproducibility
within stated dispersion allowance (rf-005); CLI smoke incl. §14.1
coverage table and structured refusals (rf-006). Unit tests cover
attainment hand-calculations, order statistics, and axes sanity.

## No-claim boundaries

- No machine performance claims live in this crate: numbers become
  citable only as ledgered rows with fingerprints (plan §14.1 discipline).
- The compute axis is compiler-achievable FMA throughput, not theoretical
  ISA peak; kernels beating it (via intrinsics fs-simd tiers) can exceed
  attainment 1.0 — reported as-is, never clamped.
- §14.1 family targets (LBM/GEMM/SpMV/FEEC/batched/FFT/rays) are
  `landed: false` until their kernels register.
- Per-CCD bandwidth axes, P/E-core-class split, frequency-state capture,
  and thermal controls are future scope (v0 measures whole-machine axes).
- Verdict gating in CI is deliberately absent on shared runners; bands
  bind only on ledgered reference machines (nightly lane, later).

## No-claim boundaries (regress)

- This module is the STATISTICS + attribution + gate arithmetic; the
  nightly both-machine CI wiring rides the ci-gauntlet pipeline
  (huq.4's lane), memory-regression tracking rides fs-alloc's
  allocation-site diffs, and FrankenPandas trend dashboards ride
  fs-report — the named seams, each consuming these verdicts.
- Suspect-commit bisection hooks are diff-vs-last-green consumers of
  the same attribution output, not re-implemented here.
