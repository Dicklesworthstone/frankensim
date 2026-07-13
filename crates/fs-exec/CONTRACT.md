# CONTRACT: fs-exec

## Purpose and layer
Two-lane executor (plan §5.2): the latency lane is asupersync's async
scheduling for orchestration; the throughput lane is a work-stealing
fork-join tile pool with weighted quanta, CCD-local-first stealing, and
fixed-shape reductions. Owns the `Cx`/`TileKernel` contract every hot
kernel programs against. Layer: L0. Depends on asupersync, fs-alloc,
fs-blake3, fs-substrate, fs-obs.

## Public types and semantics
- `Cx<'s>` — the per-tile context (plan Appendix B): `checkpoint()` /
  `is_cancel_requested()` poll the run's `CancelGate` (the MANDATORY
  tile-boundary poll), `arena()` is the tile-scoped fs-alloc arena
  (lifetime-bound; escapes are compile errors), `stream_key()` is the
  logical RNG identity, `budget()` carries asupersync's `Budget`
  vocabulary, `mode()` the `ExecMode` provenance. `refuse(TileFailure)` records
  a typed tile failure, requests sibling drain, and returns the existing
  `Cancelled` break marker without converting the refusal into a panic.
- `StreamKey { seed, kernel_id, tile, iteration }` + `key128()` — RNG
  stream identity derived from LOGICAL work identity only, never from the
  worker (Decalogue P2). fs-rand's Philox consumes the 128-bit key.
- `CancelGate` — request → drain → finalize: `request()` is idempotent and
  stamps the first-request time; workers finish their current tile, stop
  claiming, and the run returns a structured outcome. Timestamps feed
  reports only, never results.
- `TileKernel` (`type Out: Reduce; tiles() -> TilePlan; run(tile, &Cx) ->
  ControlFlow<Cancelled, Out>`) and `TilePlan { tiles, kernel }` with the
  FNV-stable `kernel_id()`.
- `Reduce` — fold identity + `merge`, applied over per-tile slots on the
  FIXED-SHAPE pairwise tree: split at the largest power of two below `n`,
  recurse — shape a pure function of the tile count, items visited in
  ascending index order, so `merge` need not be commutative. Implemented
  for `()`, `u64`, `f64`, `Vec<T>`, `reduce::Compensated`.
- `reduce` module (plan §5.4, the P2 machinery): `pairwise_fold`,
  Neumaier `Compensated` partials (Reduce-composable), `det_sum`/`det_dot`
  /`det_norm2` (256-element blocked, unfused products), `det_min`/`det_max`
  (IEEE total order), `det_argmin`/`det_argmax` (ties -> LOWEST index),
  `det_prefix_sum` (compensated sequential scan), and `audit_accumulator`
  — the G5 order-sensitivity audit that catches arrival-order bugs and
  localizes them to the smallest exposing prefix.
- `TilePool` / `PoolConfig { workers, topo, quantum_weights, seed, mode,
  arena }` — `run(&kernel)` / `run_with_gate(&kernel, &gate) -> (Result<Out,
  RunError>, RunReport)`. `run_declared_budgeted` additionally carries an exact
  caller-supplied asupersync `Budget` into every tile `Cx`; legacy run wrappers
  deliberately retain `Budget::INFINITE`. `run_declared_leased_budgeted`
  (bead wf9.16) further binds the run to a shared
  `fs_alloc::OperationMemoryLease`: root metadata (slots, deque headers and
  initial entries, range plans, victim table headers/final entries/construction
  temporary, per-worker atomics, retained pairwise-fold buffers, report
  vectors) is sized with checked arithmetic and reserved fallibly BEFORE
  worker launch (`RunError::MemoryPlanOverflow`, `MemoryRefused`, or
  `MemoryAllocationRefused`; nothing ran and no root metadata remains
  allocated) and every tile arena charges the lease
  while its chunks are held (mid-run refusals surface as `TileFailed` with
  `AllocError::LeaseExhausted` and drain like any typed refusal). Legacy
  wrappers pass an unbounded lease; the caller reads `lease.receipt()` for
  canonical accounting of the observed admission trace. Identical fixed plans
  with identical arena-cache state have deterministic cumulative demand;
  cache history and near-limit concurrent refusal/peak fields are
  state/schedule-sensitive. The enforced envelope: tracked executor root
  metadata, lease-bound arena chunks, and — since bead wf9.16.1 — output
  payload storage, because the leased entries bound
  `K::Out: LeaseAdmittedOut` (heap-bearing custom outputs fail to compile
  there; see the no-claims section for the full contract). Thread stacks,
  allocator bookkeeping, panic strings, and heap a kernel allocates and
  drops entirely within its own tile body remain explicit no-claims.
  Workers are scoped per run; each owns one contiguous `TileRun`
  (bead wf9.16.2: two u64s, zero post-launch allocation) seeded from the
  contiguous, weight-proportional `weighted_ranges` plan, stealing the
  BACK HALF of a victim's run in `victim_order` (same-CCD ring first).
  `PoolConfig::for_host` / `TilePool::for_host` select the recorded host-probe
  topology; `workers()` and `placement_identity()` expose the normalized
  execution identity that tune producers must bind. The latter has a readable
  topology/mode/pinning-intent prefix plus a derive-key BLAKE3 suffix over
  normalized workers, weights, arena policy, the `ArenaPool`'s recorded
  hugepage decision/outcome, and exact requested pin groups. Pin success is not
  claimed by this identity. Its already-shipped schema is now explicit as
  `TILEPOOL_PLACEMENT_IDENTITY_VERSION = 2` under
  `org.frankensim.fs-exec.tilepool-placement.v2`; this declaration does not
  rotate any v2 byte or root. `admit_retained_placement_identity` accepts only
  that exact producer version and a byte-for-byte identity recomputed from the
  normalized pool.
- `RunError { Cancelled, TilePanicked, TileFailed, WorkerSpawn, MemoryRefused,
  MemoryPlanOverflow, MemoryAllocationRefused, ReductionPanicked, Incomplete }`
  — structured, teaching
  outcomes with tile provenance. `RunReport` — the caller-declared `RunId`
  that keyed every tile stream, steal counts, cross-CCD steal counts,
  cancel-latency samples, `cancel_latency_p99_ns()`, canonical `to_json()`.
  Outcome-class precedence is `WorkerSpawn` > `TilePanicked` > `TileFailed` >
  `Cancelled`, preserving panic containment when a sibling also refuses.
  Within the selected tile-failure class, deterministic mode reports the
  lowest observed logical tile (and its message/failure), never mutex-arrival
  order. `TileFailure::Allocation` retains the original `fs_alloc::AllocError`
  as its error source.
- `LatencyLane` — thin configured handle on the asupersync runtime
  (`block_on`, `runtime()`); no fs-exec scheduling policy of its own.
- `victim_order(worker, workers, topo)` / `weighted_ranges(tiles, weights)`
  — pure, deterministic; these functions ARE what workers use, so fixture
  verification verifies runtime behavior.
- `Racer` / `RacerConfig` / `RaceBranch` / `RaceRun` / `BranchReport` /
  `NoWinner` — speculative races (plan §5.2 behavior 1). Victory rule:
  Deterministic = lowest-index accepted result, recomputed from OUTCOMES
  (timing moves when kills land, never who wins); Fast = first accepted
  arrival, recorded. Canonical JSON reports escape every dynamic branch name.
  Early kills: branch j dies once any i<j is accepted;
  a parent gate (kill-handle) cancels the whole tree via a bounded-stride
  watcher. Liveness caveat documented: below-leader branches must
  terminate on their own budgets before the decision seals.
- `solver` module (behavior 2): `SolverState` (in-house little-endian
  codec, floats as raw bits, self-contained bytes — no pointers, content-
  hash references — so "migrate" can someday mean another machine),
  `ResumableSolver::step` (bounded pause granularity), `drive` (pause IS
  the cancellation path), `fork` (round-trips through bytes, proving
  serializability at fork time), `StepVerdict`, `SolverProgress`.
- `Tuner` / `TuneRow` / `TuningDecision` / `TuneSource` / `ScheduleKind` /
  `TuneEvidence` / `TuneObservation` / `WallTimeSummary` / `WorkUnit` /
  `ThroughputUnit` / `TuneError` — the autotuner (plan §5.5).
  `calibrate(&probe)` first
  requires the probe's stable fingerprint to match the tuner's machine
  key, then measures a real stencil-edge sweep through the real pool (argmin with
  lowest-index tie law), reduction cost, steal cost, and selects the
  schedule kind from measured per-core bandwidth; rows are keyed kernel ×
  shape-class × MACHINE FINGERPRINT with typed evidence and a refresh
  counter (recalibration idempotent). Tune-row identity v2 is the exact,
  domain-separated canonical JSON transport; it carries evidence schema v2.
  Evidence v2 preserves every wall-time sample plus explicit observation and
  per-observation sample counts with revalidated min/max summaries, distinguishes
  completed-tile and steal counters by `WorkUnit`, and records throughput as
  a checked nearest integral milli-unit with `ThroughputUnit` (the resulting
  integer is persisted exactly). Its optional
  `candidate_separation_ppm` is emitted only by the explicit ranked-wall-time
  constructor and is only the descriptive fastest-to-runner-up gap; unranked,
  singleton, or mixed-unit rows have `null`, never a fabricated statistical
  confidence. Persistence is a strict JSON-lines file store;
  migrating it to fs-ledger requires retaining these typed fields rather than
  relabeling opaque integers. Foreign-fingerprint rows are stale and ignored
  on load. The loader accepts only tune-row identity version 2 and evidence
  version 2 in their exact declared domains and canonical writer grammar,
  requires count prefixes to match the exact arrays, and re-derives summaries
  and separation from the exact observations,
  recognized units, full-width canonical integers, strictly positive wall-time
  samples (internally measured sub-nanosecond elapsed values are represented by
  the 1 ns floor), and positive integral refresh counters; suffixes and
  alternate numeric spellings are corruption.
  Parsing and generated-row emission are bounded before growth (16 MiB store,
  1 MiB canonical row, 64 KiB string, 4096 observations per row, 4096 samples
  per timing observation, and 4096 wall-time samples in aggregate per row).
  Every locally generated row must be a canonical writer-to-parser fixed point
  before preparation, commit, insertion, or persistence. Duplicate kernel ×
  shape-class rows for the selected fingerprint are corruption rather than
  last-write-wins. `TuneRow` keeps its state private and exposes a validating
  constructor plus read-only accessors; the public `to_canonical_json` path
  revalidates the complete row before emitting exact transport bytes, so
  callers cannot construct or stamp an unadmitted row with the current
  domain/version.
  Decisions (`tile_edge_for`, `schedule`) are RECORDED as tuning-decision v1
  exact canonical JSON with a domain, identity version, kernel, canonical
  parameters, and source. The strict parser rejects field reordering, stale
  domains/versions, alternate parameter spellings, and trailing bytes. Studies
  pin decisions through typed helpers or this validating replay API and replay
  uses recorded plans, never re-tuned ones (replay fidelity). The process-local
  diagnostic history is a deterministic bounded window: at most 4096 entries
  and 1 MiB of canonical serialized decision bytes (including domain and field
  framing), with oldest-prefix batch eviction. `decision_history()` exposes the
  evicted count and `is_complete()`; a window with an evicted prefix MUST NOT be
  presented as a complete replay record. Production dispatch receipts belong
  in the Design Ledger. General decision kernel identities are nonblank and
  bounded to the canonical 64 KiB tune-string domain before cloning.
  Cold-start defaults: 8-cube tiles, bandwidth-rich schedule.
- `GemmBlockPlan` / `GemmExecutionIdentity` / `GemmTuneKey` /
  `PreparedGemmRow` / `PreparedGemmDecision` / `GEMM_KERNEL_PREFIX` — the
  MC/NC blocking lane for the parallel-GEMM consumer (bead yqug). Plans live
  on a bounded lattice
  (`mc` a multiple of 8 in [8, 1024]; `nc_cap` a multiple of 128 in
  [128, 8192]); only the canonical `mc=X,nc-cap=Y` spelling parses (pins
  fail closed otherwise). A canonical scoped key binds the producer's
  bit-semantics version, shape class, requested thread count, normalized
  maximum thread budget (not the candidate-dependent spawned-worker count),
  explicit memory limit (`u64::MAX` is the canonical unbounded class), exact
  probe dimensions, resolved ISA tier, placement policy, and
  implementation identity, plus a required producer-supplied build/codegen
  identity. The scoped-key identity is v4 in domain
  `org.frankensim.fs-exec.gemm-tune-key.v4`; it additionally frames the exact
  three-element probe count. v3 keys lack that explicit count, v2 keys lack the
  memory seam, and older keys also lack the build seam; none are accepted as
  current GEMM keys. Row lookup, pin lookup, ledger
  lookup, and the recorded decision all use that SAME scoped key, so neither a
  neighboring shape nor a different execution or build configuration can
  reuse the row or pin.
  GEMM evidence must be explicitly RANKED wall-time candidates whose labels
  are canonical plans; the selected plan must equal the minimum-time
  candidate with insertion-order tie-breaking. Cache adoption requires the
  expected `GemmTuneKey` and binds the embedded key, shape, machine, params,
  ranked evidence, and evidence argmin. Generic `adopt_row_json` refuses GEMM
  rows so a caller cannot bypass expected-key binding. `prepare_gemm_row` and
  `prepare_adopt_gemm_row_json` validate without mutation;
  `commit_gemm_row` installs only if tuner state is unchanged, permitting a
  session to persist the canonical params and row first. Likewise,
  `prepare_gemm_decision` resolves without logging and
  `commit_gemm_decision` records only after successful dispatch.
- `KillRegistry` (behavior 3, Bet 8): candidate id -> `Arc<CancelGate>`;
  `kill` (idempotent; unknown id is a non-event), `kill_where` (batch
  elimination, ascending order), `registered_gate` (fetch without silently
  creating), `kill_registered` (structured missing-gate error), `release`.
  Everything a candidate evaluates — pool runs, races, drives — shares its
  explicitly registered gate.

## Invariants
1. Completeness: a non-cancelled, non-panicked run executes every tile in
   `0..plan.tiles` exactly once (exec-001).
2. Fixed-shape reduction: the pairwise tree's shape depends only on the
   tile count, so results are bit-identical across worker counts, steal
   schedules, and repeats — proven with non-associative floats, a
   non-commutative concatenation, and compensated artifact hashes across
   {1,2,P,2P} workers (exec-002/008/009, G5).
3. Stream keys are pure functions of `(seed, kernel_id, tile, iteration)`;
   shuffling worker counts changes nothing (exec-003).
4. Cancellation is request → drain → finalize: after `CancelGate::request`,
   workers claim no new tiles, in-flight tiles finish (or observe the gate
   at their own poll points), arenas reclaim to quiescence, and the run
   returns `RunError::Cancelled` with completed/total counts (exec-004/005).
5. Panic containment: a panicking tile is caught with tile provenance,
   siblings drain via the gate, the pool remains usable, and the process
   NEVER aborts. User-defined `Reduce::merge` panics are caught separately as
   `ReductionPanicked`; OS worker-creation failures cancel and drain every
   already-started worker before returning `WorkerSpawn` (exec-005 and unit
   battery).
6. Steal order is CCD-local-first under the fixture topologies; initial
   quanta are weight-proportional within one tile (exec-006).
7. Per-tile arenas come from one `ArenaPool` (chunk-recycled); the pool's
   quiescence oracle is the leak check after every run.
8. Race losers are FULLY drained before `race` returns (scope join), their
   arenas reclaimed (quiescence oracle); the winner (index and bits) is
   identical across timing jitter in Deterministic mode (exec-010).
9. Pause -> serialize -> deserialize -> resume reproduces the
   uninterrupted solver trajectory bit-exactly at any pause depth
   (exec-011, chaotic-map witness); forks are independent and
   serialization-proven at fork time.
10. A registry kill drains the candidate's whole tree at its next poll
    points with arenas quiescent (exec-012, latency ledgered).
11. Tune rows always carry the machine fingerprint; loads drop foreign
    rows and reject non-canonical, dimensionally ambiguous, or out-of-domain
    rows; wall-time summaries and candidate separation are derived from exact,
    strictly positive samples, never trusted as independent claims;
    calibration refuses a probe whose fingerprint differs before any row
    mutation;
    duration narrowing is checked, non-finite/negative/unrepresentable
    throughput is refused before any row mutation, parser and generated-writer
    allocations are bounded, generated rows reparse identically before they
    can enter or leave tuner state, duplicate selected-machine keys are refused;
    recalibration replaces same-key rows with refresh incremented; pinned
    decisions are typed or canonical-validated and reproduce identically on
    ANY machine, calibrated or not; their in-memory diagnostic window is
    count- and byte-bounded, evicts a deterministic oldest prefix, and exposes
    incompleteness instead of silently claiming replay coverage (exec-013 and
    tuner unit battery).
12. GEMM rows and pins are scoped to the complete execution identity; imports
    match the requested scoped key and machine, params are canonical bounded
    plans, selected plans equal the ranked-evidence argmin, parameter families
    cannot shadow one another, and decisions record the exact key used by
    lookup and replay. Row and decision installation are explicit commits so
    failed persistence or cancelled execution cannot fabricate local state or
    a successful-dispatch receipt.
13. Leased runs refuse before launch when tile counts or root-byte arithmetic
    are unrepresentable, when the lease cannot admit the checked root charge,
    or when fallible root-vector reservation fails. The root charge includes
    the `pairwise_fold` split-buffer peak (`2n-1` output elements for `n>0`)
    and one concurrent victim-order partition in addition to final tables;
    it releases on every return/unwind. This invariant applies only to the
    tracked envelope named above, not arbitrary kernel-owned heap.

## Tile-pool placement identity (v2)

The placement payload remains the historical v2 typed-binary sequence: u64
normalized workers; topology CCD/core dimensions; mode tag; declared weight
count plus ordered u32 weights; normalized arena chunk, maximum, optional
limit, free-list, and hugepage-policy fields; the declared byte count plus
exact canonical hugepage-decision JSON; then declared pin-group and per-group
CPU counts plus ordered u32 CPU ids. BLAKE3 derive-key hashing uses the exact
versioned domain above, and the readable key prefix repeats the schema version,
topology, mode, and requested pinning intent.

Every normalized `PoolConfig` placement field and every recorded
`HugepageDecision` field is semantic. The study seed is deliberately excluded:
it keys logical random streams, not hardware placement or tune-row population,
and a dedicated nonmovement test locks that boundary. Requested pin groups are
semantic; observed pin success remains a timing fact and no-claim. Mutation
tests move domain, prefix stem, version, every scalar and enum field, ordered
weights/CPU ids, and each declared count prefix independently without relying
on adding or removing elements. The schema depends on
`fs-alloc:hugepage-decision`, including its strict JSON writer/parser
fixed-point surface; the placement encoder therefore never hashes unchecked
display JSON. Its internal `PlacementCounts` transport fields are also
exhaustively classified as derived, so adding a new count frame cannot bypass
the owner gate.

Retained placement rows fail closed: only producer v2 and an exact recomputed
identity are admitted. Stale/future versions and any identity mismatch require
an explicit migration or recalibration. The re-exported `TilePool` exposes the
current values as `PLACEMENT_IDENTITY_VERSION` and
`PLACEMENT_IDENTITY_DOMAIN` associated constants.

## Snapshot envelope (bead wf9.8.2, v1)

Solver snapshots travel inside a canonical envelope — magic
`FSEXSNAP`, envelope version, stable state TYPE id, payload SCHEMA
version, caller-ledgered provenance, payload length, and payload
checksum — validated in full BEFORE the payload decoder runs. Every
`SolverState` declares `TYPE_ID` (never reused, never changed) and
`SCHEMA_VERSION` (bumped on any layout change); cross-type bytes,
unknown envelope versions, stale schemas, bit flips, truncations, and
appended bytes are each a distinct structured `EnvelopeError`, never a
plausible-but-wrong decode. Schema incompatibility is an explicit
refusal (write a migration when an old version must stay readable).
Length-prefixed vector decoders also refuse wire lengths that do not fit the
reader's `usize` or whose byte extent overflows, before allocation; a 64-bit
length can never truncate into a plausible 32-bit element count. If a valid
envelope carries unconsumed schema bytes, the payload refusal reports the
decoder's exact cursor and remaining-byte count.
`seal(provenance)`/`unseal` carry the run/ledger identity;
`to_bytes`/`from_bytes` are the unattributed convenience over the same
envelope, and `fork` round-trips enveloped bytes. Pause → seal →
unseal → resume remains bit-exact (conformance-tested).

## Stream identity is declared, never scheduled (bead wf9.7.1)

`TilePool` holds NO stream-identity state: the former pool-global
iteration counter made a kernel's RNG keys depend on how many
unrelated runs the pool had executed and on concurrent invocation
order. `run`/`run_with_gate` use the fixed implicit `RunId(0)`;
re-running a kernel under a NEW logical identity (generation, trial,
restart) goes through `run_declared(kernel, gate, RunId(g))`, where
the id comes from the caller's ledger. Keys derive solely from
(study seed, kernel id, tile, declared run) — bit-identical across
pool reuse, concurrency, arrival order, and worker count, and
reconstructible from ledger fields alone (conformance-tested). The
checked width-refusing bridge into fs-rand's key type lives in
fs-rand (`StreamKey::from_exec_parts`, bridge v1) — layering forbids
the reverse dependency.

## Race drain totality (bead wf9.8.1)

`race_with_gate` is PANIC-TOTAL and hang-free: empty races are refused
with a structured `NoWinner` (empty reports + teaching message) before
any thread spawns; the branch body AND the acceptance predicate run
inside one unwind guard, and the terminal-slot epilogue (slot write +
watcher release) is panic-free by construction (poison-tolerant locks
throughout), so the parent watcher is always released — an accept
panic is a `Panicked` outcome, never a hung scope. `KillRegistry`
locks are poison-tolerant, and `kill_registered` returns a structured
`UnregisteredKill` instead of an ignorable `false` for candidates that
must be wired. `registered_gate` lets tournament admission refuse an absent
caller-owned gate without manufacturing a dummy registration; the tournament
holds the fetched `Arc`, so a concurrent registry release cannot disconnect
an admitted evaluation tree. The G4 storm test drives races under registry-owned
gates with external kills: every kill lands registered, every race
returns, arenas end quiescent.

## Error model
All fallible APIs return structured values (`RunError`, `LaneError`) with
teaching `Display` text. Kernel panics become `RunError::TilePanicked`;
reduction panics become `RunError::ReductionPanicked`; fallible scoped thread
creation becomes `RunError::WorkerSpawn` rather than an unwind;
executor-internal invariant violations become `RunError::Incomplete`
(reported, not panicked). The only intentional panics are lock-poisoning
`expect`s (reachable only after a panic already contained elsewhere) and
kernel-authored asserts, which are contained per invariant 5.
Typed kernel refusals become `RunError::TileFailed` after every sibling and
scope arena has drained; they remain distinct from cancellation and panic.
Prelaunch memory failures remain distinct as `MemoryPlanOverflow` (checked
dimension/byte arithmetic), `MemoryRefused` (the shared lease, retaining used
and limit bytes), or `MemoryAllocationRefused` (fallible root backing-store
reservation after lease admission, with the lease charge rolled back).

## Determinism class
Deterministic (P2): results and stream keys are bit-stable across runs,
worker counts, and steal schedules on the same ISA, by construction
(slot-per-tile + length-keyed pairwise tree + logical keys). Tie-breaking
law: argmin/argmax ties resolve to the lowest logical index; float
comparisons use IEEE total order. Cross-ISA: identical shapes reduce
divergence to scalar-arithmetic classes (FMA contraction, libm ULP) owned
by fs-math and reported by the G5 cross-ISA report once the second-ISA
runner lands. `ExecMode::Fast` currently shares the same reduction shape
and exists as recorded provenance for the future relaxation. Timing values (steal counts, latencies) are measurements
quarantined in `RunReport`/events, never in results.

## Cancellation behavior
The throughput lane polls the gate at every tile boundary and requires
kernels to poll `cx.checkpoint()` at bounded strides inside long tiles;
drain semantics per invariant 4. The latency lane inherits asupersync's
region state machine (request → drain → finalize) unmodified. Cancel
latency is MEASURED per run (histogram in `RunReport`, ledgered via events);
see no-claims for the 200 µs target's status.

## Unsafe boundary
None. The pool is safe Rust (scoped threads, mutex deques, atomics);
`catch_unwind` is safe containment. Lock-free deques, if they ever land,
arrive as a registered capsule with a SAFETY.md.

## Feature flags
None. Everything here is `[S]` solid-tier.

## Conformance tests
tests/conformance.rs, cases exec-001..exec-008 (JSON-line verdicts; seeded
cases carry seeds): completeness/arena hygiene, G5 bit-identity across
worker counts, stream-key worker-independence, external-cancel drain with
ledgered latency histogram, the 300-run G4 storm with panic injection,
steal-order/quanta fixtures, latency-lane responsiveness under saturation,
reduction-shape invariance, the exec-009 G5 audit (compensated
artifact hashes bit-stable across {1,2,P,2P} workers; seeded arrival-order
bug caught with prefix localization), exec-010 (deterministic race victory
under jitter + loser drain), exec-011 (bit-exact checkpoint/resume/fork on
a chaotic trajectory), exec-012 (kill-handle drains a deep tree,
latency ledgered), and exec-013 (calibrate -> persist -> consume round
trip; fingerprint keying and mismatch refusal; strict row domains;
idempotent recalibration; canonical pinned replay).
The `reduce` unit suite also runs bead 4nh8's G0 property harness: 512
shrink-armed lengths in `0..=4096`, biased around powers of two, compare the
complete `pairwise_fold` syntax tree against an independently stated
`next_power_of_two(n) / 2` recursion (seed `0xE008_0001`). Existing fixed and
G5 reduction pins remain unchanged.
tests/constellation_smoke.rs pins the
asupersync Budget vocabulary. In-module unit suites cover the gate, keys,
Reduce laws, partitioning, victim orders, self-cancellation, and pool
survival after panics, exact finite-budget propagation, simultaneous typed
allocation refusals, and mixed panic/refusal precedence. GEMM tuner unit drills cover hostile embedded cache
keys, invalid params, unranked evidence, selection/argmin disagreement,
identity-dimension isolation, exact-key replay, parameter-family collisions,
and explicit row/decision commit semantics. The pool unit suite also locks the
complete placement-v2 mutation matrix, independent count-prefix movement,
seed nonmovement, and fail-closed retained-version admission.

## No-claim boundaries
- NO 200 µs cancel-latency CLAIM yet: the reference-hardware p99 gate
  belongs to the roofline/perf harness with release builds and machine
  fingerprints; today the histogram is measured and ledgered per run, with
  a generous sanity envelope in CI (exec-004).
- NO lock-free deque claim: the v1 deques are mutex-based with the correct
  stealing PROTOCOL; Chase–Lev arrives only with roofline evidence that
  justifies its unsafe capsule.
- Workers are scoped per run by default (spawn cost ~tens of µs amortized
  over a kernel run). The parked-crew lane (`with_parked_crew`,
  `with_parked_crew_local` — bead tkr7) is the OPT-IN answer where that
  cost dominates (measured: N-D FFT axis passes, bead 27d3): a crew parks
  once inside its owner's scope and runs dispatch to it over a condvar,
  driving the SAME worker protocol — bitwise-identical results across
  all three lanes by construction. The crew's job handoff is a registered
  unsafe capsule (`crew.rs`, SAFETY.md beside it): a lifetime-erased job
  borrow revoked by a completion latch that `dispatch` blocks on.
- The throughput lane's default (`run`/`run_declared*`) still spawns
  joined `std::thread::scope` workers. Callers holding a live asupersync
  task use `run_scoped` (bead lx0e) or a task-scoped parked crew, whose
  workers are scoped CPU children of the calling task via
  `Cx::scoped_cpu` — task cancellation and budget exhaustion drain runs
  at tile boundaries. The plan's stronger claim that EVERY throughput
  lifetime is an asupersync child scope is still NOT made: untasked call
  sites keep the std lane deliberately (forcing a runtime at every
  perf-lane/tool call site would invert the layering).
- NO achieved thread-pinning/NUMA-binding claim: `victim_order` steers
  locality and supported hosts may attempt the requested affinity through the
  audited fs-substrate capsule, but v1 workers ignore the syscall result. The
  placement key therefore identifies pinning intent, not observed CPU/CCD
  placement; typed requested-versus-observed receipts remain
  `frankensim-3iq7`. P/E quantum WEIGHTS are plumbed but their values await
  the autotuner.
- `run_declared_budgeted` propagates a finite asupersync budget unchanged, but
  generic enforcement of its poll/deadline/cost dimensions is NOT claimed:
  kernels must consume the dimensions they understand. Legacy run wrappers
  still supply `Budget::INFINITE`.
- `run_declared_leased_budgeted` / `run_scoped` enforce one shared lease over
  the checked, tracked root-metadata formula, all chunks acquired by their
  leased arenas, AND output payload storage (bead wf9.16.1): both leased
  entries bound `K::Out: LeaseAdmittedOut`, the sealed inductive contract
  (scalars, fixed arrays/tuples of admitted types, `fs_alloc::LeasedVec`
  of admitted elements, and the `Concat` fold wrapper) under which every
  owned byte is inline — visible to the slot/fold accounting — or
  lease-admitted before allocation. A heap-bearing custom output
  (`Vec`/`String`/opaque payloads invisible to `size_of`) FAILS TO COMPILE
  on the leased entries (`compile_fail` battery in `admit.rs`); `Vec<T>`
  keeps its `Reduce` impl for the legacy unleased wrappers only. Kernels
  admit payloads through `Cx::lease()`; a merge-time re-admission refusal
  in `Concat` surfaces as the documented `ReductionPanicked` containment
  with every charge released by the unwind. Still NOT claimed: thread
  stacks, allocator metadata/capacity rounding, bounded panic-diagnostic
  strings, and heap a kernel allocates and drops entirely within its own
  tile body (invisible to any output/arena boundary). Work stealing is
  ALLOCATION-STABLE (bead wf9.16.2): each worker owns one contiguous
  `TileRun` (two u64s in its cache-padded slot, admitted pre-launch — the
  root formula has NO per-tile deque-entries term), and a half-steal is a
  midpoint split with the exact `ceil(len/2)`-from-the-back arithmetic of
  the previous `VecDeque::split_off` protocol, so the tile→worker transfer
  is preserved verbatim while allocating nothing after launch. The
  contiguous-run invariant that makes this exact (seed contiguous, pop
  only the front, adopt only a victim's back half) is structural; proven
  by the steal-storm battery plus a counting-global-allocator differential
  (`tests/steal_stability.rs`): materially different steal traffic at
  identical run shape shows a ZERO allocation-count delta. None of the
  no-claim allocations may be described as covered merely because their
  container header or initial entry storage is charged.
- The latency lane's ≤100 ms conversational guarantee is HELM's gate;
  exec-007 measures and ledgers turnaround without claiming it.
- `ExecMode::Fast`'s 5–15% relaxed-reduction throughput claim is NOT made:
  Fast currently shares the deterministic tree; the relaxation (and its
  measured delta) waits for the roofline harness.
- Race kill-propagation from a parent gate polls at a 50 µs stride
  (measurement-class latency, ledgered); sub-poll-interval propagation
  needs the perf harness like every other latency claim.
- Ledger SPILL of solver checkpoints (revolve-style schedules, artifact
  rows) is fs-ad/fs-ledger territory; this crate owns the snapshot bytes
  and their bit-exactness only.
- NO "calibrated is faster" assertion in CI: debug-profile timing is
  noise; the improvement is DOCUMENTED via the ledgered calibration
  report, and the perf harness owns throughput verdicts (same doctrine as
  every other latency/perf claim here).
- NO statistical-confidence claim for tune rows: evidence v2 records exact
  observations, wall-time extrema, and (when meaningful) a descriptive
  candidate-separation ratio. It does not estimate repeatability,
  uncertainty, or a probability that the selected candidate is optimal.
- Stencil prefetch-distance calibration rows arrive when that kernel registers
  its microbench. GEMM has a typed row and dispatch lane, while reference-ISA
  performance admission and broader tropical tune-next analytics (Bet 12)
  remain with the Gauntlet and fs-plan respectively.
- Per-core-class (P/E) bandwidth calibration inherits fs-substrate's
  pinning no-claim; the schedule decision uses aggregate per-core numbers.
- Deterministic hash-map wrappers are not shipped: the contract's rule is
  "no HashMap iteration order in results" (BTreeMap or index-keyed slots
  in hot paths); an enforcement lint belongs to CI tooling.
