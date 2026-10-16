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
- `lbm_core_roll_hash()` runs a canonical UNIFORM-TAU D2Q9 BGK roll
  fixture (`Grid::uniform(24, 16, 0.6)`, periodic in `x`, full-width
  walls at `y = 0` and `y = ny-1`, 50 plain `step()` calls, 6 probe
  cells) so vessel and ornithoid consumers share one public audit point
  for the collide/stream core. Its authority is exactly that path. It
  does not cover `Rheology`, `ContactModel`/free surface, interior
  rasterized-obstacle bounce-back, non-periodic inlet/outlet columns,
  momentum-exchange or partial-saturation variants, or D3Q19; a change
  confined to those paths will not move the hash. (The `xo2k` migration
  of fs-lbm's rheology `powf` paths to `det::pow` is the repo's own
  counterexample: shared machinery the vessel rides moved and
  `GOLDEN_LBM_CORE` did not.)

## Invariants

1. Content hashes are metric-only. Wall-clock seconds are evidence,
   not identity.
2. Re-running the same deterministic smoke stage must reproduce the
   same metric hash before that hash is eligible to become a golden.
   Vessel and ornith smoke companion events include the complete metric-bit
   notebook evidence so a future golden delta can be attributed field by
   field.
3. Shared-machinery changes surface once in a shared audit, within that
   audit's stated coverage. Two audits exist:
   - the LBM roll hash covers the uniform-tau collide/stream path only
     (see `lbm_core_roll_hash` above);
   - fe2e-006 is the CROSS-CONSUMER e-race audit. It covers (a) race-core
     replay determinism on a fixed normalized loss table, (b) agreement
     between each flagship's public screening wrapper and the same core
     driven under that flagship's own declared convention —
     `fs_ornith::screen_generation` at declared span 1.52 / ceiling 1.5,
     and `fs_vessel::race::screen_lips` at scale 200 with the
     data-derived support `200 x (fixture spread + jitter width)` — and
     (c) ONE shared loss table (the ornithoid's own `-L/D` scores) raced
     under BOTH declared conventions, gated on reaching the same
     SELECTION: same winner, same elimination count.

     What (c) does NOT claim is that the two conventions cost the same.
     They do not: on the shared 12-candidate table the ornithoid's
     convention spends 925 evaluations and the vessel's 859, because the
     ornithoid normalizes onto a fixed ceiling with +/-0.01 jitter while
     the vessel scales by 200 with a data-derived support and +/-5e-5
     jitter. Both counts are emitted; `cross_consumer_evals_claimed_equal`
     is a literal `false` in the row so the distinction cannot be lost by
     re-reading. Equal winner plus equal elimination count pins the whole
     survivor set only when the elimination is total (11 of 12 here);
     that condition is emitted as `cross_consumer_survivor_set_pinned`
     rather than assumed.

     The vessel side became auditable only when its convention moved out
     of `crates/fs-vessel/tests/battery.rs` into the `fs_vessel::race`
     library surface (bead `frankensim-extreal-program-f85xj.2.31`): a
     convention that exists only in a test cannot be driven by any other
     crate, which is exactly how the original case could claim a
     cross-flagship audit while invoking neither flagship.
4. Mid and full stages are wired with `#[ignore]` until their
   cadence and envelopes belong to the perf/CI lanes.
5. Failure drills must produce expected structured outcomes, and every
   field a drill row publishes must be DERIVED FROM STATE THE DRILL
   MOVED — never a literal written into the format string:
   - cancellation storm: the race completes with a surviving winner
     after mid-race kills;
   - budget exhaustion: a real refinement counter funds
     `LBM_REFINE_BUDGET = 1` LBM refinement, then the remaining 5 of 6
     candidates DEGRADE to the surrogate + conformal path; funded and
     degraded counts and the in-band count are all measured, and the
     gate fails if nothing degrades or the funded lane cannot answer;
   - ledger crash recovery: after a transaction is begun and dropped
     without commit, the reopened ledger is READ BACK — the committed
     artifact materializes byte-identical, the artifact written inside
     the uncommitted transaction is absent, `events` holds exactly the
     one committed row, and artifact integrity verifies clean.
     `schema_version().is_ok()` is recorded but is not the gate: it
     stays true both when a recovery loses the committed prefix and when
     it replays the uncommitted transaction;
   - model-form escalation: `fs_surrogate::certify_or_escalate` takes a
     real decision on the FITTED band half-width — escalating at a
     tolerance below it, serving the surrogate at a tolerance above it —
     and the escalated query is then actually served by the funded LBM
     lane. The count of conformal violations is reported as a measured
     number, not as an escalation claim.

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
- **fe2e-005** shared uniform-tau D2Q9 roll hash for the collide/stream
  path both flagships ride (coverage bounded as stated above).
- **fe2e-006** race-core replay determinism on a fixed normalized loss
  table; consumer/core agreement for BOTH public wrappers
  (`fs_ornith::screen_generation` and `fs_vessel::race::screen_lips`)
  each under its own declared convention; and cross-consumer SELECTION
  agreement over one shared loss table. Evaluation counts across
  conventions are reported, never equated.
- **fe2e-007** failure drills for cancellation storms, budget
  exhaustion, ledger crash recovery, and model-form escalation, with
  every emitted outcome field measured.
- **fe2e-008** forensic JSON row self-audit and bitwise notebook
  replay.
- `fe2e_mid_stages` and `fe2e_full_stages` are intentionally ignored
  lane placeholders until the perf/CI cadence lands.

Six further tests are DRILL FALSIFIERS. They emit no aggregate and no
forensic companion; they exist so a green drill cannot be a silent pass.
Each seeds the fault its bead names and asserts the drill's own gate
rejects it:

- `fe2e_006_consumer_core_agreement_is_falsifiable` — a drifted
  normalization ceiling or declared span must break the ORNITHOID's
  consumer/core agreement, while the pre-fix self-replay shape stays
  green under both (bead `frankensim-extreal-program-f85xj.2.31`).
- `fe2e_006_vessel_consumer_core_agreement_is_falsifiable` — the same
  for the VESSEL: drifting the declared support away from the jitter
  width (slack 1e-3, 1e-2, 1e-1) must break agreement, and a validator
  noisier than the declared support must produce a STRUCTURED
  `RaceError::PairwiseInput` refusal rather than a verdict. It also
  records the honest negative: a pure rescale (20x, 2000x) is an exact
  invariance of the pairwise e-process, not a drift, so the audit
  neither can nor should see it (bead `…2.31`).
- `fe2e_006_cross_consumer_selection_is_falsifiable` — a drifted vessel
  declared support (slack 10.0) must break cross-consumer SELECTION
  agreement over the shared table; the test also asserts that the two
  conventions' evaluation counts DIFFER, so the reported wording cannot
  silently harden into an equality claim (bead `…2.31`).
- `fe2e_007_budget_drill_counts_real_degradation` — funding every
  candidate must yield `degraded == 0` and fail the gate; a funded lane
  that cannot answer must fail it too (bead `…2.32`).
- `fe2e_007_ledger_drill_detects_a_broken_recovery` — a recovery that
  discards the committed prefix and one that replays the uncommitted
  transaction must both fail, and both are shown to leave
  `schema_version().is_ok()` true (bead `…2.33`).
- `fe2e_007_escalation_drill_takes_a_real_decision_and_spends_it` — the
  escalated query must reach the funded lane exactly once, and a lane
  that returns no answer must not be reported as served (bead `…2.32`).

The falsifiers stub the high-fidelity lane so the ACCOUNTING can be
falsified without paying for six LBM refinements; fe2e-007 itself spends
the real refinements (one funded, one escalated).

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
It builds a serial reference from ascending coordinate input and two pooled
grids from the reverse input, requiring all three to converge to the same
canonical Morton order. The pooled grids run sequentially: the primary
`TilePool` uses the host logical-CPU count, while the alternate uses half that
count rounded up. The admitted production profiles therefore exercise two
distinct worker counts without retaining both pooled grids at once. After one
zero-force BGK step per grid, all 19,000,000 published f64 population values
must be finite and bit-for-bit identical between the serial reference and each
real `TilePool` run. All three canonical mass reductions must remain within the
emitted `8 gamma_n` roundoff envelope. The pooled receipts retain four
deterministic reports -- collide then stream for each worker count -- with 1,954
completed kernel groups, per-worker completion accounting, open cancellation
gates, and placement identities that round-trip through the current
producer-version admission check.

`allocated_state_bytes()` is used only for its exact logical meaning: three
population buffers per active tile, 29,184 bytes per tile and 456,000,000 bytes
per grid. The test holds a 912,000,000-byte `OperationMemoryLease` charge before
constructing the serial grid and one pooled grid at a time, but labels it a
shadow preflight, not allocation authority. The three sequential grid
constructions account for 1,368,000,000 logical state bytes in total while the
peak retained grid state remains 912,000,000 bytes. `SparseGrid3` stores
ordinary `Vec` and `BTreeMap` allocations, activation has ordinary heap
temporaries, each exact-state oracle copies another 304,000,000 bytes
transiently (608,000,000 bytes total across both comparisons), and
`step_pooled` uses the runner's legacy internal unbounded lease. None of those
allocations is charged to the shadow receipt. The evidence therefore refuses
`sparse-state-memory-lease-authority` and separately marks only the
`shadow-memory-preflight-ledger` as restricted. It also refuses
`structured-sparse-heap-oom-refusal` before allocating: ordinary infallible
vector growth can still abort under host pressure instead of returning a typed
error. A bounded sparse allocation claim requires a lease-backed storage API,
fallible sparse construction, and a leased pooled-sweep entry point.

Sparse phase rows are likewise report-only `Custom` observations. The stable
configuration identity binds the active Morton-key set, dimensions, population
layout, BGK parameters, seeds, serial/pool activation protocols, peak and total
logical byte counts, both worker counts, both TilePool placement identities and
their producer version, two pooled runs, D3Q19 semantics version, harness
versions, and the mass-acceptance policy version, formula, population count,
multiplier, and computed-bound bits. It excludes clocks, process RSS,
unsurfaced allocator metadata, activation temporaries, and observed pin
success. Linux `VmHWM` and the macOS peak-RSS refusal keep the same semantics
as the scalar rung; no quiet-host performance or attributed RSS-budget claim
follows from them.

The eight completed aggregates emit canonical `ConformanceCase` records
with Info/Error severity,
failure-record linting, JSONL validation, and print-before-terminal-
assert ordering. Ten live forensic companions retain the prior identity
pairs by mapping `stage` to the emitter scope and `kind` to the
`Custom` name: `vessel-smoke/artifact`, `ornith-smoke/artifact`,
`frame-smoke/artifact`, `marquee/status`, `erace-audit/race`, the four
`drill/*` outcomes, and `notebook/emitted`. Their object payloads are
validated before printing. The fe2e-006 case identity is
`fe2e-006-erace-cross-consumer` (formerly
`fe2e-006-erace-core-and-ornith-consumer`, and before that
`fe2e-006-erace-audit`). The lineage is the fix history: the original
name described a cross-consumer comparison the case never performed;
the second named exactly what it then covered — the core plus the one
consumer that had a public wrapper; the current name is earned, because
`fs-vessel` now publishes `race::screen_lips` and the case really does
drive both consumers. The row no longer carries the interim
`vessel_wrapper: "none-public"` marker; `vessel_wrapper`,
`vessel_core`, `vessel_declared_span`, `shared_table_vessel` and
`cross_consumer_selection_agree` carry measured values in its place.
The constructed-only `log_row`
fixtures in fe2e-008 remain the escaping and utility-shape self-audit;
they are not live suite rows.

Input-seed provenance follows the fixtures exactly. fe2e-001,
fe2e-004, and fe2e-005 are fixed-input cases and use zero. fe2e-002
uses `0xE2E` for both its generation LCG and screening call. fe2e-003
uses ensemble input seed `90210`; the Cx stream seed `0xF1A6_5A1D` is
recorded separately as execution provenance and is never presented as
input randomness. fe2e-006 uses `0xAB` for both replayed race runs and
records `0xE2E` separately as `ornith_input_seed`, the seed of the
generation its ornithoid consumer/core agreement check screens and of
the shared cross-consumer loss table, plus `0x7E55` as
`vessel_input_seed`, the seed the vessel's screening jitter hashes with
(the same constant its own vsl-005 battery uses); the aggregate case
keeps `0xAB`.
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
  limited to one million whole-tile cells, two serial-anchored exact-state
  comparisons across distinct pooled worker counts, logical retained-state
  accounting, and four observed TilePool passes. NUMA placement, per-CCD
  bandwidth, CCD-shaped reductions, scale cancellation latency, quiet-host
  timing promotion, cross-ISA comparison, attributed total-heap/RSS coverage,
  and an admitted M4 peak-RSS budget remain explicit named skips until their
  retained host evidence exists.
- No CI authority is claimed. DSR remains the repository automation
  source of truth.
- No evidence package or FrankenScript study driver is emitted yet;
  the lab notebook is an in-crate deterministic artifact body.
- No closed-bead proof is claimed for the ignored mid/full fidelity
  lanes until their perf/CI cadence and envelopes land.
