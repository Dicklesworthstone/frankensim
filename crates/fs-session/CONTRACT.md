# fs-session — CONTRACT

Sessions, capability tokens, and the resource GOVERNOR (plan §11.3):
budgets are ENFORCED, not advisory — plus the agent-proofing trio:
idempotency keys, `estimate()` dry runs, and errors as guidance.

Ambition tags: tokens/governor/idempotency/estimate [F per the bead
label; shipped surface tested to [S] discipline].

## Purpose and layer

Layer **L6** (HELM). Runtime deps: `std`, fs-ir (admission bridge +
study parsing), fs-exec (CancelGate/SolverState), fs-ledger
(persistence), fs-plan (cost models), fs-obs (content hashing), fs-qty.
Consumers: the P2 marquee demo, the HELM e2e suite (gp3.11).

## Public types and semantics

- `CapabilityToken { session, ops globs, core_s, mem_bytes, wall_s,
  cores, ledger_scope }` — the explicit grant every IR program executes
  under; `to_admission()` bridges into fs-ir static admission (one token,
  checked statically at admission AND continuously by the governor).
- `Governor` — `Send + Sync`; hot paths are mutex-guarded in-memory
  state. `open_session` rejects non-finite or negative floating grants
  before registration. `charge(session, Charge)` rejects non-finite,
  negative, or overflowing deltas before mutating meters, then meters
  core-seconds / peak memory / wall and returns `Enforcement`: `Ok` →
  `Throttled` (at the grant) → `Paused` (past 1.2× the grant, with a
  teaching resume hint). The governor NEVER silently kills.
- `submit_once(session, idem_key, work)` — exactly-once execution:
  the first caller runs and is charged; concurrent/repeat callers block
  on a condvar and receive `Duplicate` with the SAME receipt and NO
  charge. Caller panics and invalid returned charges are contained as a
  terminal `Failed { receipt, what }`; all waiters receive that same
  failure, no charge is committed, and retry requires an explicit new
  key. `idempotency_key(agent_key, program)` = agent key + FNV content
  hash.
- `apply_memory_pressure(session, level, gate)` — the DECLARED
  degradation ladder (`LADDER`: spill coldest arenas → coarsen
  adaptively → pause-serialize-resume) fires steps `1..=level` in order;
  the pause step requests the session's `CancelGate` so the solver
  checkpoints at its next tile boundary (P7). Every event carries
  attribution and a deterministic ordinal.
- `estimate(study, cost_models, cores)` — the dry run: p10/p50/p90 wall
  from fs-plan quantile models over `:dof`/`:size` features, declared
  memory ask, energy (p50 × cores × 45 W/core), and an HONEST
  `unmodeled_ops` coverage list (never silent gaps).
- `CalibrationReport` — estimate-vs-actual rows, ratio quantiles, and a
  content-addressed ledger artifact (`estimate-calibration`): the cost
  models' own report card.
- `Guidance { code, diagnosis, fixes }` — errors as teaching:
  `from_finding` lifts fs-ir admission findings (the canonical §11.3
  `BudgetInfeasible` fixture) with their cost-model-ranked fixes intact.
- `flush_to_ledger(&Ledger)` — consumption, degradation, and idempotency
  receipts persisted as `session.*` events. Explicitly single-threaded:
  fsqlite connections are `!Send` by design.
- `gemm_f64_session(tuner, ledger?, gate, threads, m, n, k, α, a, b, β, c)`
  — the production GEMM autotune loop (bead yqug): measure → cache →
  model → cancellation-correct dispatch. The scoped key binds fs-la's bit
  semantics version, power-of-two shape class, requested/normalized thread
  budget, exact capped probe dims (M/K ≤ 512, N ≤ 2048), resolved SIMD tier,
  unpinned work-stealing placement policy, and implementation version; the
  ledger key also binds the machine fingerprint. Shape/overflow and
  pre-requested-cancel checks happen before any tune mutation; one-thread,
  small-M, and no-product calls bypass tuning. Pinned plans skip measurement;
  else an exact tuner/ledger row is used; else an up-to-4×2 lattice is
  deduplicated by the `(mc,nc)` values fs-la will actually execute and sampled
  three times. Every output word from every repeat is compared by `f64::to_bits`
  (signed zero and NaN payloads included); drift fails closed. The declared
  model is argmin of minimum wall time with lattice-order tie breaking, not a
  confidence claim. A row is prepared, durably written with a separate params
  field that must agree with its body, then committed locally. A decision is
  recorded only after fs-la drains and successfully commits the staged output.
  `GemmDispatch.kernel` is the exact replay key; replay pins the recorded key
  and params rather than reconstructing a weaker base key.

## Invariants

1. **Enforcement is structured**: every over-grant outcome is `Throttled`
   or `Paused` with resource, used, granted, and a resume hint — no kill
   path exists in the API.
2. **Exactly-once**: for any idempotency key, `work` runs at most once
   globally; all callers observe the same receipt; consumption is charged
   exactly once (16-thread race-tested).
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
7. **Invalid resources fail closed**: NaN, infinities, negative values, and
   accumulated floating-point overflow are rejected before any token or meter
   mutation. Landing exactly on a grant returns `Throttled`.
8. **GEMM tuning cannot create phantom success**: malformed shapes and
   no-op/serial routes cannot create rows or decisions; cancellation or bit
   drift during measurement cannot create a row; cache failure, foreign
   execution identity, and params/body disagreement cannot install a row.
   Ledger persistence precedes local row commit; successful compute precedes
   decision commit. Cancellation during final dispatch may retain its already
   validated measured row, but records no successful decision and leaves caller
   `C` bitwise intact.

## Error model

`SessionError`: `UnknownSession`, `InvalidResource`, `Submission`,
`Persistence`. `GemmTuneError`: cancellation with completed/total bounded tile
counts, typed tuner refusal, ledger refusal, or exact-bit drift with candidate
and repeat. Refusals that teach travel as `Guidance` values with ranked fixes.
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
dispatch poll the same gate inside bounded packing/microtile work. fs-la stages
`C`, stops claiming work, drains all scoped workers, and commits only after the
final poll; cancellation returns progress counts and leaves caller `C`
bit-for-bit unchanged.

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
round-trip; ss-005 ladder order + gate request + toy-SolverState
snapshot equality + attributed ordinal-ordered events; ss-006 the
canonical BudgetInfeasible finding surfacing as ranked `Guidance`;
ss-007 32-way adversarial-grant storm with exact meters and structured
outcomes only; ss-008 seeded caller panic with eight concurrent duplicates,
bounded completion, one shared terminal failure receipt, and zero charge;
ss-009 NaN/infinite/negative grant and charge refusal with no-mutation checks;
ss-010 the exact-grant throttle boundary and atomic accumulated-overflow
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
changing evidence labels.
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
  ledgered success/failure receipts; cross-process replay reconstruction
  belongs to the HELM e2e/crash-recovery bead (gp3.11).
- **Two-lane executor integration** (interactive vs batch lanes with
  core quotas) is deferred to gp3.11's study-scale batteries; the
  enforcement/idempotency/estimate surfaces here are what it composes.
- A mutex self-deadlock in the calibration renderer was found by the
  conformance run and fixed (single lock scope) — reentrancy is a
  documented non-assumption throughout.
- GEMM minimum-wall-time ranking is a deterministic selection rule over the
  recorded samples, not statistical confidence. The x86 oracle lane remains
  armed rather than claimed as measured until it runs on the reference host.
