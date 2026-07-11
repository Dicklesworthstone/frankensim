# CONTRACT: fs-roofline

> Status: ACTIVE (harness receipt v3). Owns the roofline measurement discipline of
> plan §14: measured axes, intensity-derived limits, dispersion-reported
> attainment, fingerprint-and-baseline-keyed ledger rows, staleness detection.

## Purpose and layer

Performance claims as falsifiable targets: benchmark every registered
kernel against its arithmetic-intensity-derived limit on the actual
machine, ledger the result under the machine fingerprint plus historical
baseline identity, and alert when either drifts. Layer: L6 (consumes fs-substrate, fs-simd,
fs-blake3, fs-exec, fs-session, fs-ledger). Runtime dependencies remain workspace
crates plus `std`.

## Public types and semantics

- `MachineAxes` (`axes::probe()`) — measured axes only: STREAM-triad
  bandwidth (single/all-core, from fs-substrate) and FMA-chain peak FLOPs
  (single/all-core, in-house microbench), plus the fs-substrate topology
  fingerprint. Never spec-sheet numbers.
- `KernelSpec` — identity + intensity model (`bytes_per_elem`,
  `flops_per_elem`), threading axis, explicit `TargetAxis`, and optional
  `target_fraction`.
- `RooflineKernel` — owns its buffers; `run_once` is the timed unit.
- `measure` / `run_registry` — warmup + repetitions →
  [`Attainment`]: median rate, achieved GB/s and GFLOP/s, binding
  `RoofSide`, binding-roof `attainment`, target-axis `target_attainment`,
  relative IQR `dispersion`, and a
  `Verdict` (`WithinBand`/`BelowBand`/`NoTarget`/`EnvironmentInvalid`).
  The invalid verdict carries a reason and is never a performance pass or
  failure.
- `run_admission_error` / `run_is_citable` — the publication boundary:
  require an explicit `AxisBaselinePolicy` whose selected baseline trusts both
  exact pre/post probes, private timed provenance, positive
  work count and sample durations, raw-sample re-derivation, unmodified
  spec/derived fields, and unique kernel identities. GEMM additionally requires
  at least one warmup and an identical sealed decision/path binding after every
  timed repetition. Analytic helper rows are deliberately non-citable.
- `finalize_registry_tuning` / `RooflineKernel::finalize_tuning` — apply that
  single aggregate admission decision to process-local state. Both outcomes
  consume the kernel's pending marker; rejection also invalidates its local
  tuner decision. Every registry hook runs even if an earlier hook fails;
  failures are returned together in deterministic registry order after cleanup
  drains. The registry length and ordered kernel/version identities must match
  the result set. Success returns a non-cloneable, one-shot
  `FinalizedRegistryRun` bound to the exact axes, baseline receipt, and result
  receipts. Publication does not happen here.
- `SECTION_14_1_TARGETS` — the plan's target table as data; `landed`
  flips only when the owning kernel registers here (no silent coverage
  gaps).
- `record_run` — one atomic ledger transaction. Admitted GEMM runs publish the
  exact sealed session row from their stable binding (fresh or adopted), plus
  metrics, `benchmark_result` events, and
  roofline `tune` rows keyed (kernel ×
  `roofline-v6:<kernel-version>:run=<finalized-receipt>:op=<operation-id>` × fingerprint LE bytes
  × 32-byte baseline hash). Rows are append-only per finalized receipt and bind
  the exact executable-content identity, operation, baseline, repetition count,
  and post-probe axes. Every measured payload is also stored as a
  content-addressed `roofline-benchmark-result` artifact and linked as an `Out`
  edge of that exact operation; the row binds the artifact hash. `record_run`
  requires and consumes the matching one-shot token. Rejected input publishes the exact baseline-admission
  event plus one rejection event and an Error op, never normal-looking metrics
  or tuning evidence; storage failures roll back the entire write set.
- `staleness` / `staleness_at` — `Fresh` / `Expired` / `ClockRollback` /
  `FingerprintDrift` / `BaselineUnavailable` / `BaselineDrift` / `BuildDrift` /
  `CorruptEvidence` / `NeverMeasured` per kernel version, fingerprint, selected
  baseline, and exact current executable. Exact current-key rows are
  semantically revalidated against canonical params, artifact bytes and output
  edge, successful operation receipt, admitted baseline, and executable
  identity. `Fresh` requires the newest current-build op to be no more than 30
  days old. `staleness_at` takes explicit wall nanoseconds for deterministic
  replay; `staleness` supplies the current clock.
- Ordered result manifest (bead gp3.15, ledger row schema v3): every recorded
  run binds a versioned `result_manifest` (ordinal × kernel × version ×
  payload content hash, canonical JSON) into the operation `ir` and folds its
  domain-separated hash into the finalized run receipt
  (`finalized-run.v2`). Staleness revalidation reconstructs the entire receipt
  from the manifest and the rows stored **today** — baseline receipt bytes,
  every payload in manifest order, manifest hash — and compares it to the
  op-bound receipt. Replacing one payload plus its matching artifact/params
  while retaining the old receipt, adding rows beyond the manifest, or
  removing/altering any sibling row classifies as `CorruptEvidence` for every
  row of that run (`tests/manifest_tamper.rs`). Pre-manifest v2 rows cannot
  prove receipt membership and are retired the same way; identical honest
  reruns stay `Fresh`.
- `kernels::default_registry` — the stable test/meta registry: fs-simd
  axpy/dot/sum (report-only bands in v0). `SeededSlowKernel` is the separate
  meta-test kernel claiming a band it cannot meet.
- `kernels::production_registry` / `GemmKernel` — the shipped command's
  registry adds real f64 GEMM through
  `fs_session::gemm_f64_session_with_pool`.
  One kernel instance owns one tuner, reusable TilePool, and cancellation gate: the first
  warmup closes measure → cache → model → dispatch and later warmups/timed
  reps reuse its validated process-local row. With `--ledger`, the kernel owns
  a dedicated read-only cache connection during measurement: it can adopt a
  valid row from a prior process, but a cold sweep remains buffered as a sealed
  `ValidatedGemmTuneRow`. The measured `Attainment` owns the publication clone;
  only `record_run`, after pre/post axes and every receipt pass
  `run_admission_error`, publishes it in the same transaction as the citable
  evidence. Execution is sequential through exclusive `&mut RooflineKernel`;
  no wrapper lock hides tune state.
- `roofline` CLI bin — axes line, per-kernel JSONL, §14.1 coverage table,
  optional `--baseline <jsonl>` plus required `--firmware <identity>`, and
  optional `--ledger` recording + composite staleness report. Omitting the
  baseline is an explicit report-only candidate run (`citable:false`).

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
   hand calculations). This remains the binding-roof report. Verdicts compare
   `target_attainment` to the declared target: GEMM's 75% row is always divided
   by measured compute peak, even when memory bandwidth is the binding roof.
3. Receipt schema v3 carries bit-exact pre-run axes, intensity spec, target
   axis, warmup count, every raw timed sample, median/p25/p75/dispersion, and
   exact derived-result bits. Rounded decimal fields are display-only. A
   standalone reader can rederive the reported rate, roof, target ratio, and
   variance bar.
4. Ledger rows are append-only per finalized run and keyed by kernel version,
   finalized receipt, fingerprint, and admitted baseline. A drifted fingerprint,
   version, baseline, or executable refuses reuse; a malformed exact-key row is
   corrupt rather than fresh. Payload bytes must hash to a retained artifact
   linked as an output of the row's exact successful op. Current-build evidence
   is fresh through 30 days inclusive, expired afterward, and a clock earlier
   than its newest op fails closed as rollback.
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
   GEMM additionally requires `warmup_runs > 0` so a cold autotune sweep cannot
   masquerade as steady-state kernel time.
8. GEMM warmup and timed repetitions never write a newly measured tune row to
   the durable ledger. The registry buffers the sealed row process-locally;
   `record_run` publishes the exact bound row in the evidence transaction,
   including when the measurement adopted it from a different ledger;
   rejection clears both
   the pending marker and local tuner state, and an already valid durable row
   remains readable without a new sweep. Finalization always consumes the
   marker so registry reuse cannot relabel an old row as newly measured.
   Fresh/adopted/local bindings insert into an empty destination, no-op on an
   identical tuple, and fail closed rather than overwrite a conflicting row.
   Replacement belongs to a separate explicit cache-refresh protocol; a cloned
   or delayed benchmark receipt has no overwrite authority. A successful commit
   consumes its fresh marker, while rollback retains it for retry without
   allowing later reuse to replace a newer cache row.
   Registry finalization additionally compares each kernel's current execution
   binding and pending-row identity with the last repetition sealed in its
   `Attainment`; results from an older run cannot preserve or publish a newer
   unfinalized kernel state. A mismatch forces rejected finalization, which
   drains that newer pending marker and invalidates its process-local decision
   rather than allowing unbound state to survive the refusal.
9. Every citable GEMM repetition binds the same exact scoped tune key, shape,
   canonical MC/NC plan, tuned source, operation-specific SIMD tier, build
   identity, derive-key-domain tune-row hash, and deterministic execution-path
   receipt. Receipt JSON embeds the canonical sealed tune-row preimage, so a
   historical benchmark remains independently re-verifiable after the mutable
   dispatch cache advances. The path proves nonempty completed TilePool panel
   traversals with sequential declared run ordinals, so NC/KC panels have
   distinct deterministic stream identities, while excluding nondeterministic
   steal, worker-distribution, and latency samples.
10. Citable admission additionally requires an explicit `AxisBaselinePolicy`
    whose selected baseline returns `Trusted` for both exact pre/post axis
    receipts. `record_run` binds the bit-exact axis records, declared
    environment/day, canonical baseline preimage, domain hash, and verdict in
    the same transaction. `None` is an explicit unbaselined candidate policy,
    never a permissive default; it publishes no metrics or tune rows.
11. The ledger op's versions field contains a domain-separated hash of the
    actual current executable bytes, not a checkout label or ambient CI
    variable. Each roofline row binds that identity and staleness revalidates it
    against both the successful op and the executable currently asking, so two
    rebuilt binaries cannot share provenance merely because their source
    revision string matches.

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
fresh/expiry/clock-rollback plus fingerprint/baseline drift, corrupt evidence,
and never-measured states (rf-004), plus rejection-without-publication (rf-004b);
re-run reproducibility
within stated dispersion allowance (rf-005); CLI smoke incl. §14.1
coverage table and structured refusals (rf-006). Unit tests cover
attainment hand-calculations, order statistics, and axes sanity. The GEMM
registry regression executes the production session call twice and proves
exactly one cold sweep plus two recorded dispatch decisions (warm-row reuse,
not a test-only wrapper). Durable-cache regressions prove all three publication
states: rejected measurement leaves zero rows, admitted evidence atomically
persists the session row plus its roofline row, and a new process adopts the
exact row identity without re-sweeping. Receipt regressions tamper every bound
decision/path field, alter only one repetition, and remove warmup; all refuse
admission. A registry-hook unit test proves the hook receives the complete
pre/post admission decision and drains every hook after a middle failure.
Baseline integration tests prove an unbaselined first run and stable sustained
contention cannot cross `run_is_citable`, while admitted ledger rows carry one
canonical `axis_baseline_admission` event.

## No-claim boundaries

- No machine performance claims live in this crate: numbers become
  citable only as ledgered rows with fingerprints (plan §14.1 discipline).
- The compute axis is compiler-achievable FMA throughput, not theoretical
  ISA peak; modest attainment above 1.0 is reported as-is. Attainment above
  1.5 means the shared axis is not credible for gating, whether because it
  was crushed/stale or because a specialized kernel outran the probe by too
  much; the run is retained as invalid evidence and must be re-probed.
- §14.1 family targets other than GEMM remain `landed: false` until their
  kernels register. GEMM is `landed: true` because the shipped production
  registry now executes the persistent session tune path; this is an
  implementation-coverage claim, not a claim that any machine met the 75%
  target.
- Per-CCD bandwidth axes, P/E-core-class split, frequency-state capture,
  and thermal controls are future scope (v0 measures whole-machine axes).
- Static floors plus pre/post agreement cannot detect a host that is already
  degraded before the first probe and remains equally degraded through the
  second. Every citable API and the CLI now require a separately admitted
  `BaselineAxes` through `AxisBaselinePolicy`; reporting-only lanes may run
  unbaselined but their measurements are candidate evidence, never citable.
- `RooflineKernel::elements()` and intensity are asserted by the registered
  implementation. Receipt v3 proves what was timed and how the arithmetic was
  derived; it does not prove a custom trait implementation performed the work
  it claimed. Default-registry review is the v1 trust root; implementation
  hashes remain follow-up scope. Likewise, public callers currently supply the
  pre/post `MachineAxes`; receipt admission proves agreement between the values
  supplied, not that a caller actually performed the post-run probe. An opaque
  production-run protocol that owns registry selection and both probes remains
  `frankensim-epic-perf-fz2.5`; custom/public registry runs are candidate
  evidence, not independently citable proof.
- Verdict gating in CI is deliberately absent on shared runners; bands
  bind only on ledgered reference machines (nightly lane, later).
- Baseline promotion authority (bead fz2.7): the `authority` module adds
  a `PromotionAuthorityVerifier` capability with typed verdicts
  (authorized / wrong-signature / unknown-key / revoked-key) over a
  `PromotionAttestation` signed on the record's content hash — which
  binds the canonical schema/domain hash, sorted source receipts,
  band/drift policy, machine identity, and promotion time, so editing
  ANY of them (including the free-text `promoted_by`) invalidates the
  signature. `AttestedBaselineStore::admit_verified` refuses
  unattested/forged/edited records, unknown or revoked keys, and
  promotions whose named source receipts are not in the retained set;
  `citable_axis_admission_authorized` (and
  `AxisBaselinePolicy::with_authority`) re-verify before EVERY citable
  decision and bind the verifying key identity into the admission
  receipt. Rotation = authorize the new key and re-promote; revocation
  retroactively demands re-promotion. TRUST CLASS, stated honestly: the
  in-tree `StaticKeyRegistry` is an operator-governed registry with
  domain-separated keyed-hash tags — tamper-EVIDENT, not unforgeable
  (its tag function is public); unforgeable signatures require an
  EXTERNAL verifier implementation injected through the same trait. The
  UNBOUND `AxisBaselinePolicy` tier remains operator-trusted and its
  receipts say `"authority":"operator-trusted"`.

## Trusted historical axis baselines (bead dfh3)

Sustained-contention detection that pre/post agreement cannot provide.
`BaselineAxes` is an opaque schema-v1 record of the machine's quiet axes with
provenance (named operator, justification, promotion day, and sorted unique
source receipt identities), an age policy (≤ 365 days), and a declared
environment identity (topology + OS + arch + firmware string, compared
verbatim). Its bit-exact canonical JSON has a domain-separated BLAKE3 identity.
Trust laws:

1. First-run measurements are CANDIDATE evidence (`Unbaselined`) — a
   probe can never authorize itself as its own baseline.
2. `promote_baseline` is the only public constructor: ≥ 3 floor-plausible,
   same-identity `BaselineCandidate` values with distinct retained source
   receipt hashes that mutually agree within the reprobe drift band, plus a
   named operator and non-blank justification. The
   promoted value is the per-axis maximum over the runs (a too-low
   baseline would inflate later attainment). Updates are re-promotions;
   no in-place mutation API exists.
3. `citable_axis_admission(pre, post, baseline, identity, now_day)`
   composes absolute floors (last-resort sanity, unchanged) + pre/post
   agreement + baseline bands. Only `BaselineVerdict::Trusted` supports
   citable gates. Distinct refusals: `Degraded` (an axis below 0.70 of
   baseline — the 6-GB/s-on-a-100-GB/s-host counterexample), `Suspect`
   (above 1.15× — not the machine the baseline describes), `Stale`
   (past the age policy), `IdentityDrift` (fingerprint/topology/OS/
   arch/firmware mismatch), `ClockRollback`, `InvalidAxes`, `ReprobeFailed`,
   and `InvalidBaseline`. Every verdict serializes as valid JSON; non-finite
   diagnostic ratios become `null`.
4. `BaselineStore` is a strict bounded JSON-lines store, one baseline
   per fingerprint in deterministic order. Admission revalidates every sealed
   invariant; malformed lines, sub-floor/impossible axes, duplicate source
   receipts/fingerprints, oversized records, and non-monotone replacement are
   refusals, not last-write-wins.

Drills (unit tests in `baseline.rs`): quiet-trusted, sustained
contention refused despite pre/post self-agreement, suspiciously-fast
refused, stale-by-age refused (boundary day still trusts), firmware and
fingerprint drift refused, first-run-not-self-authorizing, all six
promotion refusals, future/rollback clocks, valid-JSON refusal receipts, source
receipt uniqueness, and store round-trip + tamper/duplicate refusal.

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
