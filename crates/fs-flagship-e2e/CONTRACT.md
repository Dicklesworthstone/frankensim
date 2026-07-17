# fs-flagship-e2e CONTRACT

Flagship e2e suite (bead `frankensim-epic-flagships-mye.5`): staged
smoke/mid/full replay lanes for the flagship crates, cross-flagship
audits, failure drills, forensic logs, and a deterministic lab
notebook artifact. Current v2 golden constants are FROZEN: vessel
`0x4541_d7f3_2926_1082`, ornith `0xae56_945a_312e_0378`, frame
`0x9c09_b06a_7883_57fc`, and shared LBM core
`0x1539_430c_dae4_7762`. Replay equality is verified before freezing;
bump only with a semantic justification in the owning flagship or a
registered shared core.

The vessel smoke hash was formerly `0xe621_48d4_490c_a887` under the
radix-2 FFT schedule. The mixed radix-4/2 schedule intentionally
changes the floating-point operation order in `fs-cheb`'s DCT path,
which feeds the vessel stability objective. Its independent DCT,
Orr-Sommerfeld, vessel-property, and replay checks remain green; the
new bit identity is recorded here rather than silently accepted. The
metric-level audit found that only `robust_offband` moved, from
`-0.0004364607241673659` to `-0.00043646072421213883` (about
`4.48e-14` absolute); the other five metrics retained their exact bits,
and restoring only the old final-field bits reconstructs the old hash.

The subsequent mixed radix-8/4/2 schedule moved that same downstream
`robust_offband` field from bits `0xbf3c_9a98_8956_ba53` back to
`0xbf3c_9a98_894a_2018`; the other five metric bit patterns remained
unchanged. Substitution reconstructs the radix-4/2 aggregate
`0x4e42_4a53_6a63_ce8b` and the radix-8/4/2 aggregate
`0x4541_d7f3_2926_1082` exactly. The current value reproduces in debug
and release on aarch64. The upstream FFT stage path has four-quadrant
evidence; this downstream vessel aggregate still has an explicit
x86-64 no-claim until that lane is rerun.

The ornith smoke hash was formerly `0xa6fa_6460_e7c7_972f` while
`PairwiseRace` silently clipped differences against an implicit unit span.
The checked-span repair declares the fixture's analytical 1.52 support and
normalizes without clipping. The winner (11), elimination count (11), winner
L/D, certificate, and ROA retain exact bits; only the race evaluation count
moved from 394 to 925. Replacing that one metric with 394 reconstructs the
former hash exactly. The new hash records the intended, statistically valid
power/cost tradeoff rather than accepting an unattributed golden delta.

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
- `content_hash(metrics)` encodes metric names and IEEE-754 bit
  patterns through the versioned, typed `fs_obs::ident` replay format;
  its current root digest is FNV-1a-64 over those canonical bytes.
- `artifact(flagship, tier, metrics, wall_s)` constructs a
  `StageArtifact` with its content hash already computed.
- `log_row(stage, kind, payload)` constructs the utility JSON row shape
  exercised by the forensic self-audit: `stage`, `kind`, and `payload`.
  The first two fields are JSON-escaped; `payload` is a caller-supplied
  complete JSON value. Live suite evidence uses canonical `fs-obs`
  events rather than printing this utility shape directly.
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
   Vessel and ornith smoke companion events include the complete metric-bit
   notebook evidence so a future golden delta can be attributed field by
   field.
3. Shared machinery changes should surface once in the shared audit,
   not as silent drift across individual flagships.
4. Mid and full stages are wired with `#[ignore]` until their
   cadence and envelopes belong to the perf/CI lanes.
5. Failure drills must produce expected structured outcomes:
   cancellation storms, budget exhaustion, ledger crash recovery, and
   model-form escalation.

## Error model

The crate is a conformance suite, so programmer-contract violations
panic through test assertions. Completed aggregates emit canonical
`fs_obs::EventKind::ConformanceCase` records, and forensic evidence
uses validated object-shaped `Custom` companions. Evidence and
deterministic artifacts are not a recoverable application API.

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

`tests/production_scale.rs` is the ignored scale battery for
`frankensim-ei3b`. An explicit profile admits either an M4 or Threadripper host
only when debug assertions are off; the documented command uses Cargo release,
but the compiled crate cannot authoritatively name its Cargo profile or
optimization level, so the evidence calls this
`debug-assertions-off-profile-unattested`. Missing profile input produces a
named refusal without allocating, while malformed or host-mismatched input
fails after emitting the refusal. The scalar-field rung selects 128^3 on M4 and
256^3 on Threadripper. Its stable configuration identity binds the build-mode
boundary, shape, exact payload and fresh-chunk reservation, pool/lease limits,
OS, architecture, bounded model string, logical CPU count, and crate versions.
It explicitly excludes phase clocks and process RSS.

The tranche separately proves two memory properties. First, an exact
fresh-pool lease limit one byte below the preflight reservation refuses before
allocator or payload mutation: the arena and pool counters remain unchanged,
the external mutation sentinel is untouched, and the exact refusal receipt is
retained. Second, the admitted rung performs one arena allocation whose
initialization is also its first touch, sweeps every f64 cell in deterministic
index order, and drops the arena. With free-list retention disabled, the
operation lease returns to zero and pool accounting reports zero live,
reserved, and free bytes. These are logical allocator/lease claims, not process
RSS or worker-owned NUMA placement claims.

Phase rows are report-only, object-shaped `Custom` observations named
`allocate_initialize_first_touch_ns`, `serial_sweep_ns`, and
`arena_drop_reclaim_ns`. They carry the workload configuration root rather than
mislabeling it as an `fs-substrate` machine fingerprint, so cross-run
performance comparison remains refused. Elapsed durations live only in the
`Custom` payload; they are not misreported as the event's Unix-epoch
`wall_ns`. On Linux only, `/proc/self/status`
`VmHWM` is recorded under its exact process-lifetime high-water semantics; it
includes harness startup and cannot be reset, so this tranche does not use it
as an admitted RSS budget gate. On macOS, current RSS is not relabeled as peak
and the true peak claim is named-skipped.

The companion sparse-D3Q19 rung activates a near-centered 25x25x25 cube of
whole 4x4x4 tiles in a 200^3 domain, leaving 12 tile layers before and 13 after
on each axis: exactly 15,625 tiles and 1,000,000 active cells.
It builds a serial reference from ascending coordinate input and a pooled grid
from the reverse input, then requires both to converge to the same canonical
Morton order. After one zero-force BGK step, all 19,000,000 published f64
population values must be finite and bit-for-bit identical between the serial
reference and the real `TilePool` run. Both canonical mass reductions must
remain within the emitted `8 gamma_n` roundoff envelope. The pooled receipt
retains exactly two deterministic reports -- collide then stream -- each with
1,954 completed kernel groups, per-worker completion accounting, an open
cancellation gate, and a placement identity that round-trips through the
current producer-version admission check.

`allocated_state_bytes()` is used only for its exact logical meaning: three
population buffers per active tile, 29,184 bytes per tile and 456,000,000 bytes
per grid. The test holds a 912,000,000-byte `OperationMemoryLease` charge before
constructing the serial and pooled grids, but labels it a shadow preflight, not
allocation authority. `SparseGrid3` stores ordinary `Vec` and `BTreeMap`
allocations, activation has ordinary heap temporaries, the exact-state oracle
copies another 304,000,000 bytes transiently, and `step_pooled` uses the
runner's legacy internal unbounded lease. None of those allocations is charged
to the shadow receipt. The evidence therefore refuses
`sparse-state-memory-lease-authority` and separately marks only the
`shadow-memory-preflight-ledger` as restricted. It also refuses
`structured-sparse-heap-oom-refusal` before allocating: ordinary infallible
vector growth can still abort under host pressure instead of returning a typed
error. A bounded sparse allocation claim requires a lease-backed storage API,
fallible sparse construction, and a leased pooled-sweep entry point.

Sparse phase rows are likewise report-only `Custom` observations. The stable
configuration identity binds the active Morton-key set, dimensions, population
layout, BGK parameters, seeds, serial/pool activation protocols, logical byte
counts, worker count, TilePool placement identity/version, D3Q19 semantics
version, harness versions, and the mass-acceptance policy version, formula,
population count, multiplier, and computed-bound bits. It excludes clocks,
process RSS, unsurfaced allocator metadata, activation temporaries, and
observed pin success. Linux `VmHWM` and the macOS peak-RSS refusal keep the
same semantics as the scalar rung; no quiet-host performance or attributed
RSS-budget claim follows from them.

The eight completed aggregates retain their existing case identities
and emit canonical `ConformanceCase` records with Info/Error severity,
failure-record linting, JSONL validation, and print-before-terminal-
assert ordering. Ten live forensic companions retain the prior identity
pairs by mapping `stage` to the emitter scope and `kind` to the
`Custom` name: `vessel-smoke/artifact`, `ornith-smoke/artifact`,
`frame-smoke/artifact`, `marquee/status`, `erace-audit/race`, the four
`drill/*` outcomes, and `notebook/emitted`. Their object payloads are
validated before printing. The constructed-only `log_row` fixtures in
fe2e-008 remain the escaping and utility-shape self-audit; they are not
live suite rows.

Input-seed provenance follows the fixtures exactly. fe2e-001,
fe2e-004, and fe2e-005 are fixed-input cases and use zero. fe2e-002
uses `0xE2E` for both its generation LCG and screening call. fe2e-003
uses ensemble input seed `90210`; the Cx stream seed `0xF1A6_5A1D` is
recorded separately as execution provenance and is never presented as
input randomness. fe2e-006 uses `0xAB` for both replayed race runs.
The composite fe2e-007 aggregate uses zero while its companions carry
the cancellation seed `0x570`, the shared surrogate/model seed
`0x0771`, and zero for the fixed ledger drill. The composite fe2e-008
aggregate likewise uses zero while its notebook companion records the
fixed vessel input, ornith input `0xE2E`, frame input `90210`, and the
separate frame execution seed. The process ID used only to isolate the
crash-recovery database path is a resource identity, not a seed.

Setup and operation expectations can still abort before an aggregate
is reached; they remain ordinary Rust test diagnostics. Failure-drill
companions are emitted incrementally before the combined verdict. The
ignored MID/FULL lanes remain assertion-only and emit no aggregate.

Current caveat: the smoke battery is the fast replay gate for the
frozen constants above. Mid/full fidelity envelopes remain ignored
until their perf/CI cadence lands.

## No-claim boundaries

- No new vessel, ornithoid, frame, or LBM physics claim is made here;
  this crate composes public APIs from those crates.
- No production-scale full-fidelity flagship run is claimed. Mid and
  full lanes are wired as ignored tests with envelope homes. The ignored
  scalar-field rung is limited to arena/lease admission, first-touch
  initialization, serial sweep, and reclaim accounting. The sparse rung is
  limited to one million whole-tile cells, one serial-versus-pool exact-state
  comparison, logical retained-state accounting, and two observed TilePool
  passes. NUMA placement, per-CCD bandwidth, CCD-shaped reductions, a second
  production-scale pooled worker count, scale cancellation latency, quiet-host
  timing promotion, cross-ISA comparison, attributed total-heap/RSS coverage,
  and an admitted M4 peak-RSS budget remain explicit named skips until their
  retained host evidence exists.
- No CI authority is claimed. DSR remains the repository automation
  source of truth.
- No evidence package or FrankenScript study driver is emitted yet;
  the lab notebook is an in-crate deterministic artifact body.
- No closed-bead proof is claimed for the ignored mid/full fidelity
  lanes until their perf/CI cadence and envelopes land.
