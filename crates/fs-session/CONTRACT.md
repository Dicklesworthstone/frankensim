# fs-session — CONTRACT

Sessions, capability tokens, and the resource GOVERNOR (plan §11.3):
budgets are ENFORCED, not advisory — plus the agent-proofing trio:
idempotency keys, `estimate()` dry runs, and errors as guidance.

Ambition tags: tokens/governor/idempotency/estimate [F per the bead
label; shipped surface tested to [S] discipline].

## Purpose and layer

Layer **L6** (HELM). Runtime deps: `std`, fs-ir (admission bridge +
study parsing), fs-exec (CancelGate/SolverState/TilePool), fs-la (production
GEMM), fs-ledger (persistence/content hash), fs-blake3 (domain-separated
receipt identity), fs-plan (cost models), fs-obs, fs-qty.
Consumers: the P2 marquee demo, the HELM e2e suite (gp3.11).

## Public types and semantics

- `CapabilityToken { session, ops globs, core_s, mem_bytes, wall_s,
  cores, ledger_scope }` — the explicit grant every IR program executes
  under; `to_admission()` bridges into fs-ir static admission (one token,
  checked statically at admission AND continuously by the governor). fs-ir's
  current memory-ask vocabulary is an `f64` planning projection; exact byte
  enforcement remains the governor's integer-grant responsibility.
- `Governor` — `Send + Sync`; hot paths are mutex-guarded in-memory
  state. `open_session` rejects non-finite or negative floating grants
  before registration and rejects an already-open `SessionId` without
  replacing its immutable token or inheriting its meters, gate, pause state,
  or idempotency receipts. Token and optional gate registration are one atomic
  critical section. `charge(session, Charge)` rejects non-finite,
  negative, or overflowing deltas before mutating meters, then meters
  core-seconds / peak memory / wall and returns `Enforcement`: `Ok` →
  `Throttled` (at the grant) → `Paused` (past 1.2× the grant, with a
  teaching resume hint). Memory admission compares integer bytes exactly,
  including above f64's 53-bit integer precision, and evaluates the hard
  threshold as `used * 5 > granted * 6`; diagnostic f64 fields do not drive the
  verdict. The governor NEVER silently kills.
- `submit_once(session, idem_key, work)` — exactly-once execution:
  the first caller in that session runs and is charged; concurrent/repeat callers block
  on a condvar and receive `Duplicate` with the SAME receipt and NO
  charge. Caller panics and invalid returned charges are contained as a
  terminal `Failed { receipt, what }`; all waiters receive that same
  failure, no charge is committed, and retry requires an explicit new
  key. A private-field `SubmissionReceipt` is a domain-separated BLAKE3
  identity of the owning session, exact key, and terminal charge plus
  enforcement decision or failure. `Executed` and every later `Duplicate`
  expose the same `Enforcement`, so a throttled/paused charge is never hidden
  behind a generic success.
  `idempotency_key(agent_key, program)` length-frames both inputs under a
  separate BLAKE3 domain; blank or oversized supplied keys are refused.
- `apply_memory_pressure(session, level, gate)` — the DECLARED
  degradation ladder (`LADDER`: spill coldest arenas → coarsen
  adaptively → pause-serialize-resume) fires steps `1..=level` in order;
  the pause step requests the session's `CancelGate` so the solver
  checkpoints at its next tile boundary (P7). Every event carries
  attribution and a deterministic ordinal.
- `estimate(study, cost_models, cores)` — the dry run: p10/p50/p90 wall
  from fs-plan quantile models over `:dof`/`:size` features, declared
  memory ask, energy (p50 × cores × 45 W/core), and an HONEST
  `unmodeled_ops` coverage list (never silent gaps). The result is fallible:
  cores and derived wall/energy values must stay finite and non-negative, and a
  memory Count must scale to an exact positive whole-byte `u64`; zero, negative,
  fractional-byte, wrong-unit, and overflowed asks are structured refusals.
  Operation discovery includes namespaced verbs and any undotted verb present
  in the supplied model registry; registry-backed work cannot silently vanish
  from a dry-run estimate because of its spelling.
  Memory is read only from the recognized study's direct `(budget (mem ...))`
  clause: duplicate or malformed memory clauses fail closed, and a body call
  named `mem` has no authority to grant a budget.
  `unmodeled_ops` means no model exists; a present model that refuses its input,
  emits invalid numbers, or reverses its quantiles is an error, not silently
  relabeled as a coverage gap. An explicitly declared `:dof`, `:size`, or
  `:modes` feature must have one numeric value, and duplicate size features are
  refused instead of silently falling back to the unit-size default.
- `CalibrationReport` — estimate-vs-actual rows, ratio quantiles, and a
  content-addressed ledger artifact (`estimate-calibration`): the cost
  models' own report card. Row ingestion rejects negative/non-finite values and
  non-finite ratios before mutation, so its canonical JSON cannot be poisoned
  by `NaN`/infinity spellings.
- `Guidance { code, diagnosis, fixes }` — errors as teaching:
  `from_finding` lifts fs-ir admission findings (the canonical §11.3
  `BudgetInfeasible` fixture) with their cost-model-ranked fixes intact.
- `flush_to_ledger(&Ledger)` — changed consumption snapshots plus new
  degradation and terminal idempotency receipts persisted exactly once as an
  atomic `session.*` event batch. An unchanged repeated call is a no-op;
  failed persistence leaves every generation-aware cursor dirty for retry.
  The call refuses an already-open ledger transaction because it cannot know
  whether the owner will commit or roll back. Explicitly single-threaded:
  fsqlite connections are `!Send` by design. The first successful non-empty
  flush binds one governor to one owning ledger path (and exact handle for
  independent `:memory:` ledgers); a different sink is refused before writes.
- `gemm_f64_session_with_pool(tuner, cache_policy, pool, gate, m, n, k, α,
  a, b, β, c)`
  — the production GEMM autotune loop (bead yqug): measure → cache →
  model → cancellation-correct dispatch through one caller-owned, reusable
  `TilePool`. `gemm_f64_session(..., threads, ...)` is the compatibility
  wrapper that constructs an unpinned host pool. The scoped key binds fs-la's bit
  semantics version, power-of-two shape class, requested/normalized thread
  budget, exact capped probe dims (M/K ≤ 512, N ≤ 2048), resolved SIMD tier,
  the executing pool's canonical topology/mode/weights/arena/pin-groups
  identity, implementation version, and generated compiler/profile/codegen
  build fingerprint, plus `GEMM_TUNER_SCHEMA_VERSION`, which must bump whenever
  the producer lattice, probe/sample policy, ranking, or plan mapping changes;
  the ledger key also binds the machine fingerprint.
  Shape/overflow and
  pre-requested-cancel checks happen before any tune mutation; one-thread,
  small-M, and no-product calls bypass tuning. Pinned plans skip measurement;
  else an exact tuner/ledger row is used; else an up-to-4×2 lattice is
  deduplicated by the `(mc,nc)` values fs-la will actually execute and sampled
  three times. Every output word from every repeat is compared by `f64::to_bits`
  (signed zero and NaN payloads included); drift fails closed. The declared
  model is argmin of minimum wall time with lattice-order tie breaking, not a
  confidence claim. `GemmTuneCache` makes durable access explicit:
  `Disabled`, `ReadOnly`, or `ReadWrite`. Read-only callers may adopt an
  existing validated row but cannot publish during speculative work. A newly
  measured row is sealed as `ValidatedGemmTuneRow` inside `GemmDispatch`; its
  private fields cannot be forged or altered. `receipt_json()` is its canonical
  kernel/shape/machine/params/measured preimage; a public globally unique
  derive-key domain hashes those exact bytes, and
  `publish_to_ledger` participates in an already-open wider transaction.
  Cache adoption returns the same sealed identity on its first dispatch so
  downstream evidence can bind adopted and freshly measured rows uniformly.
  `publish_if_absent_or_identical` lets evidence populate an independent ledger
  but refuses to overwrite a different row already stored under that key.
  `replace_cache_row` is the distinct mutable-cache operation: only a sealed,
  remeasured row can replace stale or malformed dispatch state, and exact
  read-back is required. Replacement is never authority inherited by a delayed
  or cloned benchmark receipt. Read-write mode uses that repair path and
  durably writes the same sealed row before committing it locally. A decision is recorded
  only after fs-la drains and successfully commits the staged output.
  `GemmDispatch.run` carries the final compute counts and every real per-panel
  fs-exec `RunReport`; `execution_receipt()` projects kernel, mode, deterministic
  declared panel ordinal, and completed/total counts into
  `GemmExecutionReceipt`, explicitly excluding steal, latency, and
  worker-distribution measurements from replay identity.
  `GemmDispatch.kernel` is the exact replay key; replay pins the recorded key
  and params rather than reconstructing a weaker base key.

## Invariants

1. **Enforcement is structured**: every over-grant outcome is `Throttled`
   or `Paused` with resource, used, granted, and a resume hint — no kill
   path exists in the API.
2. **Exactly-once within the owning session**: for any `(session,
   idempotency-key)` pair, `work` runs at most once; all callers in that scope
   observe the same content-derived receipt and consumption is charged exactly
   once (16-thread race-tested). The same caller key in another session is
   independent and produces a different receipt, so one tenant cannot suppress
   another tenant's work.
3. **The ladder order is the contract**: spill before coarsen before
   pause, always; pause requests cancellation, and `SolverState`
   snapshots round-trip losslessly (pause-serialize-resume equality).
4. **Estimates state their coverage**: unmodeled ops are listed, their
   wall is excluded, nothing is silently assumed.
5. **Meters are exact under storm**: concurrent charges accumulate
   without loss (32-way storm test asserts exact totals).
6. **Every idempotency key terminates**: success or caller panic transitions
   `Pending` exactly once, wakes every waiter, and carries one shared receipt;
   failed work never charges and same-key retry never executes implicitly.
   Successful receipts bind bit-exact charge fields and the resulting
   enforcement verdict; duplicates replay that verdict without recharging.
7. **Invalid resources fail closed**: NaN, infinities, negative values, and
   accumulated floating-point overflow are rejected before any token or meter
   mutation. Landing exactly on a grant returns `Throttled`.
8. **GEMM tuning cannot create phantom success**: malformed shapes and
   no-op/serial routes cannot create rows or decisions; cancellation or bit
   drift during measurement cannot create a row; cache failure, foreign
   execution identity, and params/body disagreement cannot install a row.
   Read-only mode performs no ledger writes and exports only an unforgeable
   validated row. In read-write mode, ledger persistence precedes local row
   commit; successful compute precedes decision commit. Cancellation during
   final dispatch may retain its already validated measured row, but records no
   successful decision and leaves caller `C` bitwise intact.
9. **Session identity is immutable**: opening an existing `SessionId` is a
   structured `SessionAlreadyOpen` refusal. The original capability, meter,
   cancellation gate, pause state, and terminal idempotency generations remain
   unchanged.
10. **Flush is append-once and retryable**: one terminal submission generation,
    degradation event, or distinct meter snapshot is appended at most once by
    a governor. The whole dirty set commits atomically; refusal or failure
    advances no cursor, while successful unchanged repeats append zero rows.

## Error model

`SessionError`: `UnknownSession`, `SessionAlreadyOpen`, `InvalidResource`,
`Submission`, `Persistence`. `GemmTuneError`: cancellation with
completed/total bounded tile counts, structured TilePool failure with tile
provenance, typed tuner refusal, ledger refusal, or exact-bit drift with
candidate and repeat. Refusals that teach travel as `Guidance` values with
ranked fixes.
A caller-work panic is data, not an unwind across the governor API:
`SubmitOutcome::Failed` records its receipt and diagnosis.

## Determinism class

Governor state transitions are deterministic given the call order;
event ordinals are logical (no wall clocks in ledgered payloads).
Concurrency outcomes (who wins a race) are scheduling-dependent by
nature — the INVARIANTS above are what is guaranteed. GEMM numerical bits are
independent of the selected MC/NC plan; the wall-clock winner is inherently
environment-dependent and therefore travels as scoped evidence plus an exact
replayable decision.

## Cancellation behavior

The governor is itself a cancellation SOURCE (pause step → CancelGate).
Its own operations are short, bounded critical sections. GEMM sweep and final
dispatch use the same caller-owned pool and poll the same gate inside bounded
packing/microtile work. fs-la stages `C`, stops claiming M-band tiles, drains
all workers and Cx arenas, and commits only after the final poll; cancellation
returns compute plus TilePool progress and leaves caller `C` bit-for-bit
unchanged. The pool's worker lifetime is not yet an asupersync child scope; the
precise no-claim and follow-up live in fs-exec's L0 contract.

## Unsafe boundary

Zero `unsafe`.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs` (JSON verdicts, suite `fs-session/conformance`):
ss-001 token→admission bridge end-to-end; ss-002 Ok→Throttled→Paused
with exact meters and structured unknown-session errors; ss-003
16-thread idempotency race (one execution, one charge, shared receipt,
independent keys); ss-004 estimate p10/p50/p90 + energy + declared mem +
honest coverage, calibration ratio quantiles, ledgered artifact
round-trip; ss-004b invalid estimator/calibration numeric domains fail closed
without poisoning JSON; ss-005 ladder order + gate request + toy-SolverState
snapshot equality + attributed ordinal-ordered events; ss-006 the
canonical BudgetInfeasible finding surfacing as ranked `Guidance`;
ss-007 32-way adversarial-grant storm with exact meters and structured
outcomes only; ss-008 seeded caller panic with eight concurrent duplicates,
bounded completion, one shared terminal failure receipt, and zero charge;
ss-009 NaN/infinite/negative grant and charge refusal with no-mutation checks;
ss-010 the exact-grant throttle boundary and atomic accumulated-overflow
refusal; ss-002b duplicate session registration preserving the original token,
meters, gate, and terminal idempotency state; ss-003d atomic incremental ledger
flush, unchanged-call no-op behavior, and dirty-cursor retry after transaction
refusal.

`tests/gemm_tune.rs` (bead yqug drills): cold start sweeps once and
matches serial bits; ledger warm start seeds a fresh session without
re-measuring; stale (foreign-fingerprint) and invalid cache rows are
refused and re-measured; a requested gate cancels the sweep with no row
and untouched output; pre-requested warm and pinned paths leave output and
decision logs untouched; non-canonical/off-lattice/mis-keyed pins are
structured refusals; replay uses the actual recorded scoped decision and
reproduces the live bits; execution identity separates thread budgets and
exact probes even inside one shape bucket; serial/no-op/small and invalid-shape
paths cannot mutate tuning state; every lattice plan matches serial bits.
An n=640 producer test executes both the two-panel NC=512 route and the
single-panel wider route, proving the NC axis reaches fs-la rather than only
changing evidence labels. Caller-pool drills prove the same pool executes
measurement and final dispatch, the receipt contains real TilePool traversal,
legacy std-thread placement keys are refused, and pinned/unpinned placement
policies cannot share tune rows.
In-module injected Gauntlet tests force exact signed-zero/NaN-payload drift in
each repeat, candidate collapse, between-repeat cancellation, cache-write
failure/retry, params/body disagreement, wrong-probe adoption, and cancelled
dispatch without a success decision. The oracle lane
(`--ignored`, release) asserts the live choice is within the declared
25% tolerance of the exhaustive best-of-3 oracle at the real problem
size and reports its machine — measured 1.000/1.000/1.062 on
macos-aarch64 under ambient load; the second-ISA (x86) counterpart is
armed and runs when an x86 host picks it up.

## No-claim boundaries

- **The governor meters what it is TOLD** (`Charge` deltas from the
  executor); OS-level resource sampling and per-thread accounting are
  the executor/observability beads' territory.
- **Degradation steps are orchestration events**: actual arena spilling
  and adaptive coarsening are fs-alloc/solver behaviors triggered by
  these events, not implemented here. Pause IS wired (CancelGate +
  SolverState protocol).
- **Energy is a declared-constant model** (45 W/core), not measured
  power telemetry; the calibration channel is where reality lands.
- **Idempotency persistence is flush-based**: in-process registry +
  session-bound ledgered success/failure receipts; cross-process replay reconstruction
  belongs to the HELM e2e/crash-recovery bead (gp3.11).
- **A governor flushes to one owning ledger sink**: the in-memory cursors prevent
  duplicate appends to that sink, and a later different sink is refused before
  writes rather than receiving a partial history. Cross-ledger replication is
  an event-log concern above this API.
- **Two-lane executor integration** (interactive vs batch lanes with
  core quotas) is deferred to gp3.11's study-scale batteries; the
  enforcement/idempotency/estimate surfaces here are what it composes.
- A mutex self-deadlock in the calibration renderer was found by the
  conformance run and fixed (single lock scope) — reentrancy is a
  documented non-assumption throughout.
- GEMM minimum-wall-time ranking is a deterministic selection rule over the
  recorded samples, not statistical confidence. The x86 oracle lane remains
  armed rather than claimed as measured until it runs on the reference host.
