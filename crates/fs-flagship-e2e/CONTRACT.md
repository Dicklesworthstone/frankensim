# fs-flagship-e2e CONTRACT

Flagship e2e suite (bead `frankensim-epic-flagships-mye.5`): staged
smoke/mid/full replay lanes for the flagship crates, cross-flagship
audits, failure drills, forensic logs, and a deterministic lab
notebook artifact. This crate is the suite scaffold while the bead is
in progress; golden constants still need to be frozen before the
smoke battery is a closing proof.

## Purpose and layer

Layer **L6 (HELM)**. This crate composes the existing flagship and
support crates (`fs-vessel`, `fs-ornith`, `fs-frame`, `fs-lbm`,
`fs-race`, `fs-ledger`, `fs-scenario`, `fs-marquee`, `fs-exec`) into
one system-level e2e surface. It is not a new physics solver. Its
claim is orchestration, replay identity, cross-flagship consistency,
and structured failure evidence.

## Public types and semantics

- `Tier` names the staged fidelity lanes: `Smoke` for the fast gate,
  `Mid` for nightly-scale envelopes, and `Full` for weekly or
  on-demand production-scale envelopes.
- `StageArtifact` records a flagship name, tier, metric stream,
  content hash, and wall-clock duration. The hash is computed only
  from deterministic metrics; wall time is logged but excluded from
  identity.
- `content_hash(metrics)` folds metric names and IEEE-754 bit
  patterns through a fixed FNV-64 stream.
- `artifact(flagship, tier, metrics, wall_s)` constructs a
  `StageArtifact` with its content hash already computed.
- `log_row(stage, kind, payload)` emits the suite's structured JSON
  row shape: `stage`, `kind`, and `payload`.
- `notebook(artifacts)` emits the deterministic lab-notebook body
  over stage hashes and metric bit patterns.
- `lbm_core_roll_hash()` runs a canonical D2Q9 roll fixture so vessel
  and ornithoid consumers share one public audit point for the LBM
  core.

## Invariants

1. Content hashes are metric-only. Wall-clock seconds are evidence,
   not identity.
2. Re-running the same deterministic smoke stage must reproduce the
   same metric hash before that hash is eligible to become a golden.
3. Shared machinery changes should surface once in the shared audit,
   not as silent drift across individual flagships.
4. Mid and full stages are wired with `#[ignore]` until their
   cadence and envelopes belong to the perf/CI lanes.
5. Failure drills must produce expected structured outcomes:
   cancellation storms, budget exhaustion, ledger crash recovery, and
   model-form escalation.

## Error model

The crate is a conformance suite, so programmer-contract violations
panic through test assertions. Runtime evidence is emitted as
structured log rows and deterministic artifacts rather than a
recoverable application API.

## Determinism class

Smoke-stage identity is deterministic by construction: fixed seeds,
fixed metric order, fixed hash function, and wall time excluded from
the golden body. Stochastic or long-running future stages must use
envelopes rather than pretending wall-clock or sample-path identity.

## Cancellation behavior

The suite itself is synchronous. Cancellation behavior is tested
through lower-level public surfaces, especially `fs_exec::KillRegistry`
inside the e-race failure drill.

## Unsafe boundary

`unsafe_code = "deny"` through workspace lints. This crate introduces
no unsafe code and no unsafe capsules.

## Feature flags

None. Mid and full fidelity stages are gated by ignored tests rather
than Cargo features.

## Conformance tests

`tests/e2e_battery.rs` defines the suite:

- **fe2e-001** vessel smoke-stage hash replay and mass-drift gate.
- **fe2e-002** ornithoid smoke-stage hash replay.
- **fe2e-003** frame smoke-stage hash replay.
- **fe2e-004** marquee lane status recording; the suite records the
  owning lane status instead of pretending a disabled runner.
- **fe2e-005** shared LBM D2Q9 roll hash for vessel/ornithoid shared
  core behavior.
- **fe2e-006** e-race consistency over identical normalized losses.
- **fe2e-007** failure drills for cancellation storms, budget
  exhaustion, ledger crash recovery, and model-form escalation.
- **fe2e-008** forensic JSON row self-audit and bitwise notebook
  replay.
- `fe2e_mid_stages` and `fe2e_full_stages` are intentionally ignored
  lane placeholders until the perf/CI cadence lands.

Current caveat: smoke golden constants are still placeholders in the
initial scaffold. The crate compiles, but the non-ignored smoke
battery is not yet the closing proof until those constants are
measured, reviewed, and frozen with rationale.

## No-claim boundaries

- No new vessel, ornithoid, frame, or LBM physics claim is made here;
  this crate composes public APIs from those crates.
- No production-scale full-fidelity flagship run is claimed. Mid and
  full lanes are wired as ignored tests with envelope homes.
- No CI authority is claimed. DSR remains the repository automation
  source of truth.
- No evidence package or FrankenScript study driver is emitted yet;
  the lab notebook is an in-crate deterministic artifact body.
- No closed-bead proof is claimed while `mye.5` is still in progress
  and the smoke goldens remain placeholders.
