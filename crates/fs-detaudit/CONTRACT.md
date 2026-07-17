# CONTRACT: fs-detaudit

The G5 determinism-audit harness (bead frankensim-epic-gauntlet-6nb.6):
same-ISA bit-identity audits across worker counts and repeats with
first-divergence localization over staged hash traces, the cross-ISA
divergence classification report, and measured ExecMode
fast/deterministic deltas.

## Purpose and layer

Layer UTIL (test/audit support). Zero runtime dependencies; real audit
subjects (fs-exec pools, fs-rand streams, the PV/OED study) enter through
dev-dependency integration tests and through callers registering
closures.

## Public types and semantics

- `fnv1a64(bytes)` — the artifact/content hash helper.
- `StageHash { label, hash }`, `StagedTrace` (`final_hash()`) — one run's
  ordered stage hashes over exact result bits.
- `Subject { name, run: Fn(workers) -> StagedTrace }` — a workload
  claiming ExecMode::Deterministic reproducibility.
- `WorkerMatrix` — `host_default()` gives the bead's canonical
  {1, 2, P, P+2, 2P} (deduplicated); `explicit(counts)` for fixtures.
- `AuditConfig { matrix, repeats }`; `audit(&Subject, &AuditConfig) ->
  AuditReport` — every (worker count, repeat) trace must equal the
  baseline bit-for-bit; `AuditReport::identical()`, `json_lines()`.
- `DivergenceLocator { workers, repeat, first_stage, stage_label,
  baseline_hash, observed_hash }` — the FIRST differing stage.
- Cross-ISA: `IsaLedger`/`LedgerRow { hash, value_bits }`,
  `DivergencePolicy` (artifact → `DivergenceClass::FmaContraction` or
  `LibmUlp { max_ulps }`), `classify_cross_isa(a, b, policy) ->
  CrossIsaReport` with `clean()` and `render_markdown()`.
- `measure_mode_delta(name, repeats, det, fast) -> ModeDeltaReport` —
  wall times, gain ratio, and per-mode hash reproducibility, observed.
- `cross-isa-report` bin: ledger JSONL × 2 + policy TSV → markdown
  report; non-zero exit when not clean.

## Invariants

- The baseline is the first matrix entry, repeat 0; comparison is whole-
  trace bit equality — a matching final hash with a differing interior
  stage still diverges.
- The locator names the first differing stage index and label; length
  mismatches locate at the first absent stage.
- Every divergence is reported; nothing is averaged, retried, or hidden.
- Cross-ISA: artifacts absent from the policy must match bit-for-bit;
  libm-ULP classifications are VERIFIED against value bits when supplied
  (finite, same-sign; otherwise unclassified); envelope violations,
  missing rows, and undeclared divergences are never clean — and
  reduction-shape divergence has no admissible declaration in
  deterministic mode.
- Mode deltas are measured per call; no throughput figure is assumed.

## Error model

The library surface is infallible: audits and classifications return
typed reports; defects become report rows (`Unclassified`, divergence
locators), never silent successes. Constructor panics are programmer
errors only (`WorkerMatrix::explicit` on empty/zero). The report bin
exits 2 on usage/IO errors and 1 on a not-clean report.

## Determinism class

Deterministic: the engine's own arithmetic (FNV-1a, trace comparison,
classification, rendering) is pure and platform-independent.
`measure_mode_delta` wall times are measurements, excluded from any
bit-stability claim. Subjects carry their own determinism claims.

## Cancellation behavior

None: the engine runs subject closures synchronously; long-running
subjects manage their own Cx lanes (the PV subject does exactly that in
its own closure).

## Feature flags

None.

## Unsafe boundary

No `unsafe` (workspace forbids it; nothing here needs it).

## Conformance tests

- `tests/detaudit.rs` — deterministic subject audits clean across an
  explicit matrix; the bead's seeded arrival-order reduction is CAUGHT
  and localized to the reduce stage with repeat-consistent divergence;
  cross-ISA fixtures classify every declared category with a clean
  render; undeclared divergence, violated ULP envelopes, and missing
  rows are never clean; mode-delta reproducibility loss is observed.
- `tests/real_subjects.rs` — fs-exec pooled non-associative reduction
  and fs-rand logical stream bit-identical across the host matrix; the
  PV/OED demo study replays bit-identically; plus the ignored
  `emit_isa_ledger` per-host emitter feeding the report bin.

## No-claim boundaries

- The harness audits what subjects hand it: no claim of nondeterminism
  absence beyond the exercised matrix, repeats, and staged hashes.
- The cross-ISA report is classified-divergence evidence, not a
  cross-ISA equality certificate, and is only as complete as the ledger
  rows fed to it.
- No thread scheduling of its own; worker counts are passed to subjects,
  which own their pools.
- No nightly scheduling: the refresh cadence for the report artifact is
  CI/ops scope; this crate provides the reproducible generator command.
- No throughput claim for ExecMode::Fast: deltas are per-call
  measurements on the caller's workload and host.
