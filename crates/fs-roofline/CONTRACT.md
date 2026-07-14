# CONTRACT: fs-roofline

> Status: ACTIVE (ledger row v4, production protocol v3). Owns the roofline measurement discipline of
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
- `RooflineKernel` — owns its buffers; fallible `run_once` is the timed unit.
- `measure` / `run_registry` — warmup + repetitions →
  [`Attainment`]: median rate, achieved GB/s and GFLOP/s, binding
  `RoofSide`, binding-roof `attainment`, target-axis `target_attainment`,
  relative IQR `dispersion`, and a
  `Verdict` (`WithinBand`/`BelowBand`/`NoTarget`/`EnvironmentInvalid`).
  The invalid verdict carries a reason and is never a performance pass or
  failure. Both entry points are fallible and reject work before invoking a
  kernel: warmups are bounded to 1000, repetitions to `1..=1000`, registries
  to 256 kernels, and aggregate kernel invocations to 250,000 with checked
  arithmetic. Result and sample buffers use fallible reservations.
- `AxisBaselinePolicy` / `AttestedAxisBaselinePolicy` /
  `AxisAdmissionSnapshot` — the historical-axis trust boundary. A plain
  `AxisBaselinePolicy` can only freeze a `candidate` snapshot and is structurally
  report-only, even when its numerical verdict is `Trusted`.
  `AttestedBaselineStore::policy_for_run` is the only public mint for an owned
  attested policy: it checks exact run identity, membership of every named
  source hash in the protected retention-inventory declaration, and one atomic
  `ConfiguredPromotionAuthority` decision before measurement. The public mint
  obtains its epoch day internally; callers cannot backdate admission. Its
  consumed `decide` call obtains the day again and freezes that decision with
  the exact pre/post probes. An unavailable clock or a policy held across an
  epoch-day boundary yields an unauthorized snapshot rather than extending the
  old decision. Finalization, output, and ledger recording reuse those canonical snapshot
  bytes verbatim.
- `ExternalPerfGateLane` / `ExternalPerfGateReceipt` /
  `record_external_perf_gate_at_path` — the shared positive-recording boundary
  for the FEEC and FFT performance lanes. The configured durable path is named
  once by `EXTERNAL_PERF_GATE_LEDGER_ENV`
  (`FRANKENSIM_ROOFLINE_LEDGER`); an empty path, `:memory:`, or any SQLite
  `file:` URI is outside this deliberately narrow durable-path boundary. The
  gate JSON is bounded to `MAX_EXTERNAL_PERF_GATE_JSON_BYTES`, must be one
  valid object for the typed lane (`feec-gate` or `fft-gate`), has at most 64
  top-level fields and 64 nested containers, and must carry exactly one top-level
  `citation_eligible:true`, `recorded:true`, `report_only:false`, and the exact
  supplied admission receipt. Candidate, stale, denied, malformed,
  mismatched-lane, and report-only inputs are refusals. The Ledger-typed
  `record_external_perf_gate_in_ledger` is the in-memory test seam; FEEC/FFT use
  the path wrapper and need no direct fs-ledger dependency.

  One recorder-owned transaction stores the exact axis-admission bytes as a
  content-addressed `In` artifact and the exact final-gate JSON as a
  content-addressed `Out` artifact, links both to one successful operation,
  appends one lane-qualified diagnostic event projecting that op, both hashes,
  and the exact gate, and returns a structured receipt containing the local
  op/event ids, lane, and both hashes. Artifact write receipts, exact bounded
  byte re-reads, role-qualified edges, and every operation-envelope field must
  all agree before commit. The append-only event API has no general payload
  reader; the recorder therefore checks its positive row id and an exact
  one-row event count increment in-transaction. Those checks prove insertion,
  not the stored event payload; the event is non-authoritative, while the
  content-addressed artifact and op IR are the byte-exact re-read authority. A
  caller-owned transaction is refused. Admission freshness is checked again
  after all writes and immediately before commit, so a UTC-day
  rollover rolls back the complete gate instead of extending yesterday's
  authority. Both the operation start and completion wall timestamps must be
  monotone and map to that same decision day.
- `run_admission_error` / `run_passes_measurement_admission` — the aggregate
  measurement boundary: require a frozen snapshot whose selected baseline
  trusts both exact pre/post probes, private timed provenance, positive
  work count and sample durations, raw-sample re-derivation, unmodified
  spec/derived fields, and unique kernel identities. Passing this predicate is
  necessary but not sufficient for citation: only the sealed production
  protocol supplies production-registry provenance. GEMM additionally requires
  at least one warmup and an identical sealed decision/path binding after every
  timed repetition. Execution binding v4 binds the complete producer-owned
  `fs_session::GemmExecutionReceipt::receipt_identity`, including every logical
  memory-plan field; the embedded `execution_path` JSON is a diagnostic
  projection, while its adjacent child identity is the authoritative complete
  receipt. Analytic helper rows are deliberately non-citable.
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
  `roofline-v8:<kernel-version>:run=<finalized-receipt>:op=<operation-id>` × fingerprint LE bytes
  × 32-byte baseline hash). Rows are append-only per finalized receipt and bind
  the exact executable-content identity, operation, baseline, repetition count,
  and post-probe axes. Every measured payload is also stored as a
  content-addressed `roofline-benchmark-result` artifact and linked as an `Out`
  edge of that exact operation; the row binds the artifact hash. Receipt-backed
  production runs also retain the exact canonical fs-la normal/build graph as
  an `fs-la-depgraph-receipt` artifact, link it as an `In` edge, and bind both
  its content hash and exported-domain digest into the operation IR and every
  row. A development-salt production run is report-only: it records a
  structured rejection and publishes no metrics or tune rows. `record_run`
  requires and consumes the matching one-shot token. Rejected input publishes the exact baseline-admission
  event plus one rejection event and an Error op, never normal-looking metrics
  or tuning evidence; storage failures roll back the entire write set.
  The public custom-registry path uses the disjoint
  `roofline-candidate-v1:<kernel-version>:run=...` namespace; exploratory rows
  therefore cannot satisfy or poison the production `roofline-v8` scan. Retained
  v7 rows coexist as historical data but are outside the v8 production keyspace
  and cannot satisfy a current lookup.
- `staleness` / `staleness_at` — `Fresh` / `Expired` / `ClockRollback` /
  `FingerprintDrift` / `BaselineUnavailable` / `BaselineDrift` / `BuildDrift` /
  `CorruptEvidence` / `NeverMeasured` per kernel version, fingerprint, selected
  baseline, and exact current executable. Exact current-key rows are
  semantically revalidated against canonical params, artifact bytes and output
  edge, successful operation receipt, admitted baseline, and executable
  identity. Operation fields are metadata-preflighted by fs-ledger before
  materialization; measured artifacts are capped at the exact tune-row length
  and dependency receipts at fs-la's producer-owned 1 MiB ceiling, with a
  source-pin test that fails on independent producer drift. A bound
  refusal is corrupt evidence, never an engine failure. Validation also
  requires the canonical `production-v3` stamp,
  well-formed nonce/pre-probe/post-probe fields, exact dependency artifact
  kind/metadata/bytes, an `In` lineage edge, and a digest recomputed with
  `fs_session::GEMM_DEPGRAPH_RECEIPT_DOMAIN`. Historical rows validate against
  their own retained receipt; equality with the receipt compiled into the
  current executable is required only for a row claiming that current build.
  A public `custom-registry` row is outside the production namespace and thus
  classifies as `NeverMeasured`, never `Fresh` or `CorruptEvidence`. `Fresh`
  requires the newest current-build op to be no more than 30 days old. The
  successful operation envelope is exact: session, seed, versions, wall-time
  budget, roofline capability, `ok` outcome, absent diagnostic, and monotone
  timestamps all revalidate.
  `staleness_at` takes explicit wall nanoseconds for deterministic replay;
  `staleness` supplies the current clock.
- Sealed production protocol: `production::ProductionProbe` / `ProductionRun`
  is the only API path to a production citation candidate. `observe()` performs
  the pre-probe and mints a per-run nonce (callers may read the axes for
  baseline selection but never supply them). `run()` accepts only an opaque
  `AttestedAxisBaselinePolicy`, owns production registry selection, timed
  warmup/repetitions, the post-probe (observed strictly after the timed loop),
  one frozen admission snapshot, aggregate admission, and tune finalization.
  `run_report_only()` accepts `AxisBaselinePolicy` and returns a distinct
  `ReportOnlyProductionRun` with no `citation_eligible()` method and no ability
  to publish production metrics or tune rows. Both `record()` methods consume
  the run and retain the exact snapshot; only the attested type stamps the
  operation `ir` with `"protocol":"production-v3"`, the nonce, content hashes of both
  observed axis receipts, and the dependency-receipt artifact/digest. Its
  `citation_eligible()` pre-commit predicate additionally refuses the default
  development salt with a structured report-only reason; CLI output cannot set
  `citable:true` until the atomic record succeeds and every row revalidates as
  `Fresh`. The public `record_run` path stays available for
  harness tests but is stamped `"protocol":"custom-registry"` and is
  explicitly NON-CITABLE — a custom kernel wearing a production name,
  replaying a cloned execution binding, or passing the pre-probe twice can
  never wear the production stamp through this API (`tests/production_seal.rs`,
  in-crate drifted-post/finalizer-failure/probe-ordering battery). Type opacity
  seals the API constructor, not the database: the nonce is a process-unique
  challenge, and callers with general `fs-ledger` mutation authority remain a
  trusted-writer boundary. External package/ledger authentication is required
  before evidence can claim authority against a malicious writer.
  `ProductionRunConfig::validate` runs before registry allocation/timing and
  bounds `n` to `1..=2^24`, warmups to `0..=63`, repetitions to `1..=64`, and
  at most 64 combined warmup and timed invocations per kernel. The generic
  custom-registry measurement API retains its separate 1,000-per-field limits.
  Before constructing any production registry buffer, a
  checked integer model derives every shipped kernel's work from the actual
  production shape: axpy + dot + sum vector traffic and FLOPs, plus the square
  GEMM's `24 * side^2` logical bytes and `2 * side^3` FLOPs where
  `side = max(isqrt(n), 256)`. The complete run is refused above `2^39`
  modeled FLOPs or `2^33` modeled logical bytes. These orthogonal limits retain
  the 32 MiB streaming default while preventing both large-shape amplification
  and thousands of minimum-shape GEMM/receipt invocations; all additions and
  multiplications are checked.
  The CLI also caps promotion probes at 1000 and baseline age at
  36,500 days before allocation or loops. Promotion reads the existing store
  through the same 1 MiB bound as runtime admission, serializes same-store
  writers with an OS lock, and opens each uniquely named same-directory staging
  generation with create-new semantics. It syncs the complete generation,
  verifies open-handle/path identity and single-link ownership, atomically
  replaces the prior identity-checked regular file, and syncs the parent
  directory. It never reopens or truncates stale crash generations;
  symlink/special-file stores and staging paths are refused. This durable
  promotion transaction currently fails closed on non-Unix hosts, where the
  required file identities are unavailable. Existing and staged stores must
  each have exactly one hard link. Atomic replacement assumes a stable, trusted
  parent directory: the identity checks detect path/object drift at each check
  point but cannot prevent a privileged writer from swapping the directory
  entry after the final rename. The executable content identity is
  captured before timing and rehashed immediately before recording; drift
  refuses the transaction.
- Ordered result manifest (bead gp3.15, retained in ledger row schema v4): every recorded
  run binds a versioned `result_manifest` (ordinal × kernel × version ×
  payload content hash, canonical JSON) into the operation `ir` and folds its
  domain-separated hash into the finalized run receipt
  (`finalized-run.v3`). Kernel/version identifiers are bounded to 1..=128
  canonical ASCII bytes and JSON-escaped before any admission decision, so a
  refused custom registry cannot corrupt the retained operation envelope.
  Staleness revalidation reconstructs the entire receipt
  from the manifest and the rows stored **today** — baseline receipt bytes,
  every payload in manifest order, manifest hash — and compares it to the
  op-bound receipt. Replacing one payload plus its matching artifact/params, or
  removing/altering any sibling row, classifies every manifested row in that
  run as `CorruptEvidence`. A forged row added beyond the manifest is itself
  corrupt, while untouched manifested rows remain `Fresh`; an unrelated extra
  row is not allowed to revoke valid evidence. The crate-private
  `production::tests` battery executes these receipt-backed attacks. The
  external battery actively proves public custom-registry rows cannot acquire
  Freshness and retains the now-vacuous historical attacks as ignored source.
  Pre-manifest/dependency-receipt rows cannot prove current membership and are
  retired the same way; identical honest reruns stay `Fresh`.
- `kernels::default_registry` — the stable test/meta registry: fs-simd
  axpy/dot/sum (report-only bands in v0). `SeededSlowKernel` is the separate
  meta-test kernel claiming a band it cannot meet. Every built-in constructor
  and registry is fallible: vector sizes must be in `1..=2^24`, each GEMM
  matrix must contain at most `2^24` elements, and GEMM worker budgets must be
  in `1..=4096`. These checks and checked extent arithmetic run before buffer
  allocation; buffers use fallible reservations. Empty midpoint indexing and
  hostile-size capacity panics therefore cannot enter `run_once` through a
  built-in constructor.
- `kernels::production_registry` / `GemmKernel` — the shipped command's
  registry adds real f64 GEMM through
  `fs_session::gemm_f64_session_with_pool`.
  `production_registry_work` is the crate-private source of truth for the
  registry's four per-kernel intensity rows; sealed configuration admission
  evaluates it before allocation, rather than treating vector length as a GEMM
  work proxy.
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
  optional `--baseline <jsonl>`, `--firmware <identity>`,
  `--authority-policy <tsv>`, `--retained-receipts <lowerhex-lines>`, and
  optional `--ledger` recording + composite staleness report. Plain stores are
  auto-detected and always call `run_report_only`. An attested envelope calls
  `run` only when both strict, bounded configuration files parse and
  `policy_for_run` admits the exact machine. Missing, partial, malformed,
  denied, revoked, tampered, or cross-machine input is measured through
  `run_report_only`, emits `citation_eligible:false`, and retains its candidate
  snapshot in the same `--ledger`; there is no separate ledger setting. A
  configuration refusal longer than the durable diagnostic ceiling is retained
  as a UTF-8-safe prefix plus its original byte length and a domain-separated
  digest, so diagnostic size cannot discard an already completed report-only
  measurement.

- `regress` module (plan §14.4, bead fz2.4): the regression layer.
  `gate` — DISPERSION-AWARE bands (k·σ against the rolling baseline,
  never a naive threshold), evaluated in a shared normalized scale. An exact
  zero-dispersion baseline yields z=0 for equality and signed infinity for a
  change, so tiny finite regressions cannot disappear behind an absolute
  epsilon. A red arrives WITH its diagnosis: the
  phase-share flame-graph diff ranked by growth. `Cusum` — the
  complementary slow-drift detector (slack k, threshold h) over
  expanding-baseline standardized scores. `slower_this_month` — the
  canonical fallible dashboard question as ONE call: (kernel, drop %, guilty
  phase), refusing non-finite or negative thresholds. Its percentage uses a
  common opening/trailing scale and its phase diagnosis compares those exact
  same seven-night windows; sustained regressed middle nights cannot dilute the
  causal phase shift. Public histories,
  kernels, phase identities/counts, per-history observations, dashboard-wide
  aggregate nights/phase observations, and standardized streams have explicit
  deterministic caps. Night histories use
  strictly increasing unique logical indices (gaps allowed); duplicates or
  reversals are invalid evidence. Calibration is meta-tested: zero false alarms across 20
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
   Every retained artifact read reaches bytes only through fs-ledger's shared
   bounded path; an oversized op row or artifact cannot be materialized merely
   to decide that it is corrupt.
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
10. Citable admission requires an `AttestedAxisBaselinePolicy` minted by
    `AttestedBaselineStore::policy_for_run`, followed by a frozen snapshot whose
    selected baseline returns `Trusted` for both exact pre/post axis receipts.
    The snapshot binds the bit-exact axis records, declared environment/day,
    canonical baseline preimage and hash, attestation, sorted retained-source
    set, internally observed mint/decision day, and atomic authority verdict
    plus policy receipt. `AxisBaselinePolicy`
    (including `None`) is explicitly candidate/report-only and cannot publish
    production metrics or tune rows.
11. The ledger op's versions field contains a domain-separated hash of the
    actual current executable bytes, not a checkout label or ambient CI
    variable. Each roofline row binds that identity and staleness revalidates it
    against both the successful op and the executable currently asking, so two
    rebuilt binaries cannot share provenance merely because their source
    revision string matches.
12. An external FEEC/FFT gate can claim `recorded:true` only after the shared
    recorder accepts an authority-admitted snapshot and a lane-matched positive
    recording envelope. The exact snapshot (including attestation key,
    signature, authority-policy receipt, and complete sorted source-receipt
    set) and exact final gate remain separately content-addressed and connected
    to one operation and projected by its lane-qualified diagnostic event. Any
    validation, write, re-read, edge, event receipt/count, operation, commit,
    or same-day revalidation failure is one atomic refusal.

## Error model

Measurement APIs return `Result`: invalid resource envelopes are refused before
any kernel execution. Sample reservations happen before the affected kernel's
warmup, and a registry aborts rather than executing a kernel whose measurement
buffers cannot be reserved. A `run_once` failure is annotated with the kernel
identity, warmup/timed phase, and invocation index; it aborts the measurement
without constructing an attainment row. Any registry execution or lifecycle
failure invokes every kernel's idempotent, non-publishing `abort_tuning` hook;
process-local tune rows and decision bindings from a tokenless partial run
cannot survive into registry reuse. The production GEMM propagates session,
tune-key, and execution-receipt refusals through this path rather than panicking.
Observations such as zero rates remain successful measurements with invalid
evidence normalized to finite JSON plus an explicit reason. Ledger interaction returns `fs_ledger::LedgerError`
(structured, machine-actionable). The external gate recorder uses the same
error vocabulary for invalid lane/payload/admission, non-durable paths, caller
transactions, and exact write/read mismatches. The CLI refuses malformed arguments with
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
attainment hand-calculations, order statistics, axes sanity, exact-cap
dependency evidence accepted by both roofline and fs-plan, and a retained
cap+1 dependency classified as `CorruptEvidence` by roofline while fs-plan
returns the typed pre-materialization refusal. The GEMM
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
contention cannot cross `run_passes_measurement_admission`, while admitted ledger rows carry one
canonical `axis_baseline_admission` event.
Attested-store transport drills pin the self-round-tripping bound: string-aware
balanced-object parsing treats braces/quotes/backslashes inside strings as data,
multibyte UTF-8 fields round-trip byte-identically, raw control bytes and
surrogate or non-hex `\u` escapes refuse at parse, non-canonical spellings
refuse against the reserialized envelope, the exact `MAX_BASELINE_STORE_BYTES`
transport is admitted while parse- and admission-side limit+1 refuse, and every
refusal (including malformed attestation fields at `admit_verified`) leaves the
store unchanged.

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
  second. Every citable API and the CLI therefore require a separately attested
  `BaselineAxes` through `AttestedBaselineStore::policy_for_run`; plain or
  unbaselined measurements are candidate evidence, never citable.
- `RooflineKernel::elements()` and intensity are asserted by the registered
  implementation. The measurement receipt proves what was timed and how the arithmetic was
  derived; it does not prove a custom trait implementation performed the work
  it claimed. The sealed production registry and executable content identity
  are the implementation trust root. Public callers still supply custom
  pre/post `MachineAxes`; those rows are exploratory and never classify Fresh.
  `ProductionProbe` owns both probes and registry selection for citable rows.
- The sealed production preflight bounds the shipped kernels' declared logical
  FLOPs, logical bytes, and invocation count. It is not a wall-clock, energy,
  physical-memory-traffic, or instruction-count certificate; one bounded cold
  tuning sweep and machine-axis probes add control/calibration work outside the
  four registry intensity rows.
- Dependency receipts are operator-observed, structurally validated evidence,
  not Cargo-authenticated invocation proofs. Protocol v3 proves
  that the exact bytes compiled into fs-session were retained, rehashed, linked
  to the operation. The receipt command admits only byte-identical results from
  two complete Cargo-discovery-plus-package-hashing derivations, so detected
  movement fails closed instead of producing a mixed snapshot. This is a
  coherence check over an operator-observed filesystem, not filesystem
  transactional atomicity; undetectable ABA mutation remains outside the claim.
  A current-build row is compared with the current
  executable's exported receipt; a foreign-build row remains valid history
  when its own retained artifact and digest agree. Correspondence between the operator-selected Cargo tree and
  the actual rustc invocation remains an explicit no-claim.
- Verdict gating in CI is deliberately absent on shared runners; bands
  bind only on ledgered reference machines (nightly lane, later).
- The external gate recorder certifies atomic retention and correspondence of
  an authority-admitted axis snapshot, typed lane envelope, and exact gate
  bytes. It does not independently re-run the FEEC/FFT kernel, validate either
  lane's FLOP/traffic model, prove the target threshold, authenticate a
  malicious general ledger writer, or provide post-record revocation and
  freshness. Those measurement semantics remain owned by the lane tests and
  their crate contracts; the recorder deliberately permits a retained citable
  gate whose measured target is false.
- Baseline promotion authority is configured at the entrypoint through
  `ConfiguredPromotionAuthority::from_text`. Its bounded canonical TSV policy
  names exact authorized messages and revoked key IDs; every atomic
  `PromotionAuthorityDecision` binds the verdict and a hash of those exact
  policy bytes. `AttestedBaselineStore::policy_for_run` selects one exact
  machine baseline, verifies the envelope against that policy, requires every
  canonical source hash in the strictly sorted protected retention inventory,
  and mints one owned policy. Missing records, edited envelopes, wrong authorization,
  unknown or revoked keys, missing sources, and identity drift are refusals.
  The resulting admission receipt freezes baseline, attestation, sources,
  authority verdict, policy receipt, and probe pair. This snapshot proves the
  decision made for that run; it does not provide later live revocation or
  post-record freshness. The mint and consumed decision each read the system
  epoch day; clock failure or a day-boundary mismatch freezes an unauthorized
  snapshot. `AxisBaselinePolicy` stays structurally `candidate`
  and never carries positive authority.
- The retained-source input is a protected operator declaration of hash
  membership, not a source-artifact reader. `fs-roofline` does not fetch or
  rehash the named bytes at admission time and makes no independent claim that
  those bytes remain retrievable; artifact-backed availability proof belongs to
  a later retention service boundary.

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
3. `candidate_axis_admission(pre, post, baseline, identity, now_day)`
   composes absolute floors (last-resort sanity, unchanged) + pre/post
   agreement + baseline bands for a report-only numerical assessment.
   `BaselineVerdict::Trusted` is still non-citable until an opaque attested
   policy freezes the same assessment with positive authority. Distinct
   refusals: `Degraded` (an axis below 0.70 of
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

Every public regress entry point screens its floating and collection inputs before
any verdict arithmetic: `gate` returns `GateVerdict::Invalid { reason }`
— never Green — for non-finite or negative attainment or phase
durations anywhere in the history and for unusable specs (non-finite
or non-positive k_sigma, or a min_baseline outside 2..=4095); `Cusum::first_alarm`
alarms AT the first non-finite residual (NaN previously reset the
shortfall via `max`, silently suppressing detection) and an invalid
detector spec cannot certify quiet; fallible `standardize` refuses an oversized
stream and maps an admitted history to −∞
from the first non-finite entry so poison never enters the expanding
baseline; duplicate or reversed logical nights likewise invalidate `gate`;
`slower_this_month` rejects malformed global thresholds and reports poisoned
or non-chronological kernels FIRST with an infinite drop and the flaw as the
"why", never skipping them. Mean/dispersion, relative-drop, and phase-share
normalization use scaled arithmetic, so extreme finite values cannot overflow
into a silent NaN or turn an improvement into an invalid regression. Exact
zero-dispersion changes remain scale-independent rather than using an absolute
attainment floor.
`standardize` applies the same rule to its expanding baseline: equal values
score zero, a nonzero residual against exact zero dispersion saturates at the
sign-appropriate finite f64 bound, and positive improvements therefore cannot
overflow into poison or a false one-sided CUSUM alarm.
Phase medians include zero for nights where a phase is absent; sparse phases
cannot select their baseline only from nights where they happened to appear.
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

## Bounded authenticated staleness-history checkpoints (bead vm3i)

The `checkpoint` module seals one exhaustive staleness verification into a
chained checkpoint so sustained probing stops re-paying O(H·(K+4)) ledger
queries per kernel:

- `checkpoint_staleness_history(ledger, kernel, version, fingerprint,
  baseline)` runs the exhaustive per-row verifier once and appends a
  checkpoint row to the tune table under the reserved kernel
  `roofline-staleness-checkpoint:<kernel>` with a v2 shape prefix (invisible to
  production row queries, machine-keyed identically to the rows it covers).
  Free-form kernel, version, and shape bytes are lowercase-hex framed in the
  canonical JSON body, so JSON delimiters cannot collide. Per covered row it
  retains a domain-separated content hash over the
  row's full stored identity (kernel ‖ shape_class ‖ machine ‖ params ‖
  measured, length-prefixed), the validated build identity, the op-bound
  dependency-receipt digests, the recorded timestamp, and a verdict.
  Checkpoints chain: `digest_i = blake3_domain(chain-domain,
  prev_digest ‖ body)`; ordinals are dense from 0 and the row's params
  restate the expected digest. Insertion is `tune_put_if_absent` plus a
  read-back equality check, so a colliding ordinal can never overwrite
  sealed history. Load, comparison, insertion, and a final strict chain
  re-read share one transaction owned by the sealing call. Sealing refuses an
  already-open caller transaction, because without savepoints a post-insert
  verification failure could otherwise leave a row for the caller to commit;
  every owned error path attempts rollback and reports a rollback failure
  explicitly.
- Both durable digest layers are registered identities rather than implicit
  helper hashes. `fs-roofline:staleness-row-content` v1 binds the exact five
  stored tune-row fields with ordered u64 length prefixes.
  `fs-roofline:staleness-checkpoint-chain` v2 binds the canonical body schema,
  kernel/version/ordinal, previous-digest presence and value, ordered entry
  count and every entry field; it declares the row-content, executable-build,
  and dependency-receipt child identities. Mutation batteries cover every
  semantic field, lowercase-hex framing, domain, and version, while strict
  `load_chain` canonical reserialization is the transport guard. The identity
  closure explicitly binds the content-hash byte/hex/parser primitives.
- Retained v1 shapes under the same reserved kernel remain binding evidence.
  Their frozen v1 parser and digest are verified first. Because the historical
  v1 writer authenticated snapshots but not snapshot transitions, migration
  conservatively normalizes the union of every v1 entry: shapes never
  disappear, only byte-identical valid metadata remains valid, and any corrupt
  verdict, valid rewrite, or temporary omission becomes a canonical permanent
  tombstone. Later v1 snapshots cannot move that tombstone's row hash or
  restore it to valid. A v2 genesis restarts its own ordinal at zero, names the
  exact verified v1 tip as `prev`, and carries that normalized union. Any later
  v1 append, truncation, or mutation breaks the bridge; v1 history cannot be
  laundered by an upgrade.
- `staleness_at_checkpointed(...)` mirrors `staleness_at` exactly on the
  selection lattice (NeverMeasured/FingerprintDrift/BaselineUnavailable/
  BaselineDrift decided identically when no exact-machine seal exists). With a
  verified chain, covered rows are checked by content hash and replayed
  through the build/dependency scan from sealed metadata — two bounded
  read queries total, independent of history depth (tested at 40 retained
  runs). Rows newer than the seal (the delta) run the full exhaustive
  validator. Any anomaly FAILS CLOSED: a missing, gapped, duplicate,
  reparse-failing, digest-mismatching, or params-inconsistent chain routes
  the probe to the exhaustive path, and a covered row that was altered or
  removed classifies CorruptEvidence. Exact-machine history is loaded before a
  no-row lattice return, so deletion of *every* covered production row is
  CorruptEvidence rather than NeverMeasured.
- Tombstones are permanent: a row sealed `corrupt` stays corrupt in every
  later checkpoint (re-sealing inherits tombstones before re-validating),
  and restoring the original bytes never un-corrupts the fast-path verdict
  — an operator who trusts checkpoints keeps seeing the incident.
- Sealing over a chain that fails verification is refused
  (`LedgerError::Invalid`): a broken chain is permanent evidence, never
  extended and never overwritten. Each new snapshot is monotonic: prior valid
  rows cannot disappear or change bytes, prior corrupt rows remain tombstones
  even when their production row is absent or restored, and new rows are
  exhaustively verified.
- No-claim: the chain authenticates against *later mutation of covered
  history*; it is tamper-EVIDENT, not unforgeable. A writer with tune-table
  access can mint a fresh chain for uncovered history; what they cannot do
  is alter rows under an existing seal (or truncate them) without the next
  checkpointed probe classifying CorruptEvidence. Suffix truncation of the
  newest checkpoint itself still requires an independently retained head
  anchor to detect; the in-ledger chain makes no such external-retention claim.
  The shared reserved-kernel read is bounded over the aggregate of v1 and v2
  rows across every machine and version for that production kernel: 1024 rows
  and 16 MiB total, not 1024 rows per individual chain. An unrelated-machine
  or unrelated-version flood can therefore make sealing/probing refuse early;
  an exact machine/prefix query and compaction are future work and must preserve
  tombstone permanence.
